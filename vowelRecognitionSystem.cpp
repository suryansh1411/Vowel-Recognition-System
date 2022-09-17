#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

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
    vector<long double> frame;

    //Open input File
    fstream fd_input;                         //file descriptor
    fd_input.open("test.txt", ios::in);
    string line;

    //Read input File
    while(getline(fd_input, line))
    {
        long double val=stof(line);
        frame.push_back(val);
    }   

    fstream fd_window;
    fd_window.open("Hamming_window.txt", ios::in);
    vector<long double> hamming_window;

    while(getline(fd_window, line))
    {
        long double val=stof(line);
        hamming_window.push_back(val);
    }
    apply_window(frame, hamming_window);

    vector<long double> correlation_coefficients;
    for(int i=0;i<=12;i++)
    {
        correlation_coefficients.push_back(compute_correlation(frame, i));
        cout<<fixed<<correlation_coefficients[i]<<" ";
    }
    cout<<endl;

    vector<long double> lpc_coefficients;
    compute_levensionDurbins(correlation_coefficients, lpc_coefficients);
    for(int i=0;i<12;i++) cout<<lpc_coefficients[i]<<" ";
    cout<<endl;

    vector<long double> cepstral_coefficients;
    compute_cepstral_coefficients(lpc_coefficients, cepstral_coefficients);
    for(int i=0;i<12;i++) cout<<cepstral_coefficients[i]<<" ";
    cout<<endl;
    

}