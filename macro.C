#include <string>
#include <array>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TF1.h>
#include <TGraph.h>
#include <TObject.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TH1D.h>

std::array<double,100> unStringCSV(string xString){
    std::array<double,100> x{0};
    char c; string val; int j = 0; double xval;
    for(int ii = 0; ii <= xString.length(); ii ++){
        c = xString[ii];
        if(c == ' ') continue;
        if((c == ',') || (ii > xString.length()-1)){
            x[j] = stod(val); val = "";
            j++;
        }
        else val = val + c;
    }
    return x;
}

string extract(string line, char del1 = '[', char del2 = ']'){
    string result =line.substr(line.find(del1)+1,line.find(del2,line.find(del1)+1)-line.find(del1)-1);
    return result;
}

string cut(string line, char del1 = '[', char del2 = ']'){
    line.erase(0,line.find(del2,line.find(del1)+1)+1);
    return line;
}

void macro(string configFileName = "settings.example"){
    ifstream configFile;
    configFile.open(configFileName, ios::in);
    if(configFile.is_open()) cout<<"Reading configuration file.\n";
    string line;  while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the first interesting line

    TString finName = extract(line); // We obtain the imput filename
    std::cout<<"Opened input root file: "<<finName<<"\n"; TFile *fin = new TFile(finName);  // We open the root filename

    getline(configFile,line); const int N = stoi(extract(line)); // we get N

    while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the next interesting line
    double_t mean[N]; double_t Emean[N]; double_t sigma[N]; double_t Esigma[N];
    double_t Res[N]; double_t ERes[N];

    TString histName; TString histTitle; string set; double_t xmin; double_t xmid; double_t xmax;
    int i = 0;
    while (i<N){
        getline(configFile,line);
        if (line[0] != '>') break; // If we finished reading those values, we exit.
        histName = extract(line, '"','"'); line = cut(line,'"','"');
        histTitle = extract(line, '"','"'); line = cut(line,'"','"');
        
        TH1D *hist = (TH1D*)fin->Get(histName);  hist->SetTitle(histTitle);  // We fetch the histogram from the file
        std::cout<<"     Reading histogram "<<histTitle<<"\n";
        while (true){
            if(line.find("[") == -1) break;
            set = extract(line); line = cut(line);
            //std::cout<<"      >Performing fit on peak: "<<set<<"\n";
            char opt = set[1];
            if (set[0] == 'G'){
                xmin = stod(extract(set,',',',')); set = cut(set,',',','); xmax = stod(extract(set,',',','));
                TF1 *funfit = fitGEL(hist, xmin, xmax, opt);
                    mean[i] = funfit->GetParameter(0); Emean[i] = funfit->GetParError(0);
                    sigma[i] = funfit->GetParameter(1); Esigma[i] = funfit->GetParError(1);
                    i++;

            } else if(set[0] == 'D'){
                xmin = stod(extract(set,',',',')); set = cut(set,',',','); xmid = stod(extract(set,',',',')); set = cut(set,',',','); xmax = stod(extract(set,',',','));
                
                TF1 *funfit = fitGEL2(hist, xmin, xmid, xmax, opt);
                    mean[i] = funfit->GetParameter(0); Emean[i] = funfit->GetParError(0);
                    sigma[i] = funfit->GetParameter(1); Esigma[i] = funfit->GetParError(1);
                    i++;

                    mean[i] = funfit->GetParameter(7); Emean[i] = funfit->GetParError(7);
                    sigma[i] = funfit->GetParameter(8); Esigma[i] = funfit->GetParError(8);
                    i++;
            }
        }
        TCanvas *c1 = new TCanvas();
        hist->Draw();
    }

    while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the next interesting line
    TString xAxisName = extract(line);

    double x[N]; double Ex[N]; std::array<double,100> arr;

    while(getline(configFile,line)) if (line[0] == '>') break; line = extract(line);
    arr = unStringCSV(line); for(int i = 0; i <N; i++) x[i] = arr[i];
    
    while(getline(configFile,line)) if (line[0] == '>') break; line = extract(line);
    arr = unStringCSV(line); for(int i = 0; i <N; i++) Ex[i] = arr[i];


    for(i = 0; i < N; i++){
        Res[i]=235.5 * sigma[i]/mean[i];
        ERes[i] = Res[i] * (Esigma[i]/sigma[i] + Emean[i]/mean[i]);
        std::cout<<"E value: "<<x[i]<<"+-"<<Ex[i]<<"\n";
    }

    TGraphErrors *grMean = new TGraphErrors(N,x, mean,Ex,Emean); grMean->SetTitle("Mean peak position depending on " + xAxisName);
    grMean->GetXaxis()->SetTitle(xAxisName);     grMean->GetYaxis()->SetTitle("ADC channel"); 
    TCanvas *cGrMean = new TCanvas(); grMean->SetMarkerStyle(2); grMean->Draw("AP"); // We draw the energy-ADC channel graph

    TGraphErrors *grRes = new TGraphErrors(N,x, Res,Ex,ERes); grRes->SetTitle("Resolution");
    grRes->GetXaxis()->SetTitle(xAxisName);     grRes->GetYaxis()->SetTitle("Resolution"); 
    TCanvas *cGrRes = new TCanvas(); grRes->SetMarkerStyle(2); grRes->Draw("AP"); // We draw the energy-ADC channel graph
}

