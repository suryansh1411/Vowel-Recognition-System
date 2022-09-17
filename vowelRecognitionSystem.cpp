#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

#define PI 3.14159265

// Pre_processing => handles microphone artifact, dc correction and normalisation
long double pre_processing(vector<int>& inData, vector<long double>& processedData)
{
	long double _avg=0.0;             //dc shift = Sum(data)/Size(data);
	long double mx=0.0;               //store maximum absolute value

	int siz=0;

	for(int i=3200;i<inData.size();i++)  //neglect 32 frames (0.2 seconds)
	{
		siz++;
		_avg += ((long double)inData[i]-_avg)/(long double)siz;  //calculate dc shift

		//get absolute maximum value
		if(inData[i]<0.0)
		{
			if(mx<inData[i]*-1.0)
			{
				mx=((long double)inData[i]*-1.0);
			}
		}
		else
		{
			if(mx<inData[i])
			{
				mx=(long double)inData[i];
			}
		}
	}

	for(int i=3200;i<inData.size();i++)
	{
		long double datapoint=( (long double)inData[i]- _avg)*((long double)5000.0/mx);     //dc correction and normalisation
		processedData.push_back(datapoint);
	}

	return (5000.0/mx);
}

void get_frames(vector<long double>& processedData, vector<vector<long double>>& frames)
{
    long double mx=0.0;               //store maximum absolute value
    int mx_index=320;

    for(int i=320;i<processedData.size()-320;i++)
    {
        if(processedData[i]<0.0)
        {
            if(mx<processedData[i]*-1.0)
            {
                mx=((long double)processedData[i]*-1.0);
                mx_index=i;
            }
        }
        else
        {
            if(mx<processedData[i])
            {
                mx=(long double)processedData[i];
                mx_index=i;
            }
        }
    }

    for(int i=0;i<5;i++)
    {
        int starting_index=mx_index-80*(i);

        for(int j=0;j<320;j++)
        {
            frames[i][j]=processedData[starting_index+j];
        }
    }
}

void apply_window(vector<long double>& frame, vector<long double>& window)
{
    for(int i=0;i<frame.size();i++)
    {
        frame[i]*=window[i];
    }
}

long double compute_correlation(vector<long double>& data, int lag)
{
    long double correlation=0.0;
    for(int index=lag; index<data.size(); index++)
    {
        correlation+=data[index]*data[index-lag];
    }
    return correlation;
}

void compute_levensionDurbins(vector<long double>& correlation_coefficients, vector<long double>& lpc_coefficients)
{
    long double e=correlation_coefficients[0];
    long double k;
    vector<long double> aux;

    for(int i=1; i<=12; i++)
    {

        k=correlation_coefficients[i];
        for(int j=1; j<i; j++)
        {
            k-= lpc_coefficients[j-1]*correlation_coefficients[i-j];
        }
        k/=e;
    
        for(int j=1; j<i; j++)
        {
            aux[j-1]=lpc_coefficients[j-1]-k*lpc_coefficients[i-j-1];
        }    
        aux.push_back(k);

        lpc_coefficients = aux;
        e = (1-k*k)*e;

    }
}

void compute_cepstral_coefficients(vector<long double>& lpc_coefficients, vector<long double>& cepstral_coefficients)
{
    for(int i=0;i<12;i++)
    {
        cepstral_coefficients.push_back(lpc_coefficients[i]);
        for(int j=0; j<i; j++)
        {
            cepstral_coefficients[i] += ((double)(j+1)/ (double)(i+1)) * cepstral_coefficients[j] * lpc_coefficients[i-j-1];
        }
    }
}



int main()
{
    //Training
    vector<char> vowels = {'a', 'e', 'i', 'o', 'u'};

    fstream fd_window;
    fd_window.open("Hamming_window.txt", ios::in);
    vector<long double> hamming_window;
    string line;
    while(getline(fd_window, line))
    {
        long double val=stof(line);
        hamming_window.push_back(val);
    }

    vector<long double> raised_sine_window;
    for(int i=0; i<12; i++)
    {
        long double aux=sin((PI*(long double)(i+1))/12.0);
        raised_sine_window.push_back((1.0+6.0*aux));
    }


    fstream fd_out;
    fd_out.open("coefficients.txt", ios::out);

    vector<vector<long double>> vowel_coefficients;

    for(int i=0; i<5; i++)
    {
        vector<long double> avg_coefficients(12, 0.0);
        int size=0;

        for(int j=1; j<=20; j+=2)
        {

            string inputFileName="Data/190101089_";
            inputFileName.push_back(vowels[i]);
            inputFileName.push_back('_');
            inputFileName=inputFileName+to_string(j);
            inputFileName=inputFileName+".txt";

            //Open input File
            fstream fd_input;                         //file descriptor
            fd_input.open(inputFileName, ios::in);

            //Read input File
            vector<int> inData;
            while(getline(fd_input, line))
            {
                int val=stof(line);
                inData.push_back(val);
            }   

            //dc correction and normalisation
            vector<long double> processedData;
            pre_processing(inData, processedData);

            //get Frames
            vector<vector<long double>> frames(5, vector<long double> (320));
            get_frames(processedData, frames);

            for(int v=0;v<5;v++)
            {
                apply_window(frames[v], hamming_window);
                //get correlation
                vector<long double> correlation_coefficients;
                for(int j=0;j<=12;j++)
                {
                    correlation_coefficients.push_back(compute_correlation(frames[v], j));
                }

                //get lcp coefficients
                vector<long double> lpc_coefficients;
                compute_levensionDurbins(correlation_coefficients, lpc_coefficients);

                //get cepstral coefficients
                vector<long double> cepstral_coefficients;
                compute_cepstral_coefficients(lpc_coefficients, cepstral_coefficients);

                apply_window(cepstral_coefficients, raised_sine_window);

                size++;
                for(int k=0; k<cepstral_coefficients.size();k++)
                {
                    avg_coefficients[k]+=(cepstral_coefficients[k]-avg_coefficients[k])/size;
                }
                // fd_out<<endl;
            }
        }

        for(int i=0;i<12;i++)
        {
            fd_out<<avg_coefficients[i]<<" ";
        }
        fd_out<<endl;

        vowel_coefficients.push_back(avg_coefficients);
    }

  
    vector<long double> tokhuras_weights = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    //testing
    for(int i=0;i<5;i++)
    {
        for(int j=2;j<=20;j+=2)
        {
            string inputFileName="Data/190101089_";
            inputFileName.push_back(vowels[i]);
            inputFileName.push_back('_');
            inputFileName=inputFileName+to_string(j);
            inputFileName=inputFileName+".txt";

            //Open input File
            fstream fd_input;                         //file descriptor
            fd_input.open(inputFileName, ios::in);

            //Read input File
            vector<int> inData;
            while(getline(fd_input, line))
            {
                int val=stof(line);
                inData.push_back(val);
            }   

            //dc correction and normalisation
            vector<long double> processedData;
            pre_processing(inData, processedData);

            //get Frames
            vector<vector<long double>> frames(5, vector<long double> (320));
            get_frames(processedData, frames);

            vector<long double> avg_coefficients(12, 0.0);
            int size=0;
            
            for(int v=0;v<5;v++)
            {
                apply_window(frames[v], hamming_window);
                //get correlation
                vector<long double> correlation_coefficients;
                for(int j=0;j<=12;j++)
                {
                    correlation_coefficients.push_back(compute_correlation(frames[v], j));
                }

                //get lcp coefficients
                vector<long double> lpc_coefficients;
                compute_levensionDurbins(correlation_coefficients, lpc_coefficients);

                //get cepstral coefficients
                vector<long double> cepstral_coefficients;
                compute_cepstral_coefficients(lpc_coefficients, cepstral_coefficients);

                apply_window(cepstral_coefficients, raised_sine_window);

                size++;
                for(int k=0; k<cepstral_coefficients.size();k++)
                {
                    avg_coefficients[k]+=(cepstral_coefficients[k]-avg_coefficients[k])/size;
                }
            }
            
            long double mn_distance=INT_MAX;
            int mn=-1;
            for(int p=0; p<5; p++)
            {
                vector<long double> tokhuras_distance(12);
                for(int c=0;c<12;c++)
                {
                    // cout<<"TI"<<endl;
                    // cout<<(avg_coefficients[c]-vowel_coefficients[p][c])<<endl;
                    tokhuras_distance[c]=(avg_coefficients[c]-vowel_coefficients[p][c])*(avg_coefficients[c]-vowel_coefficients[p][c]);
                }

                apply_window(tokhuras_distance, tokhuras_weights);

                long double distance=0.0;
                for(int d=0;d<12;d++)
                {
                    distance+=tokhuras_distance[d];
                }

                if(distance < mn_distance)
                {
                    mn=p;
                    mn_distance=distance;
                }
            }

            cout<<vowels[mn]<<endl;
        }
    }    

}