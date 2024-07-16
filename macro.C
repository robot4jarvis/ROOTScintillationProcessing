#include <string>
#include <array>
#include "functions.h"
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
                std::cout<<"|"<<xmin<<"|"<<xmax<<"|\n";

            } else if(set[0] == 'D'){
                xmin = stod(extract(set,',',',')); set = "," + cut(set,',',','); 
                xmid = stod(extract(set,',',',')); set = cut(set,',',','); xmax = stod(extract(set,',',','));
                
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
    }

    TGraphErrors *grMean = new TGraphErrors(N,x, mean,Ex,Emean); grMean->SetTitle("Mean peak position depending on " + xAxisName);
    grMean->GetXaxis()->SetTitle(xAxisName);     grMean->GetYaxis()->SetTitle("ADC channel"); 
    TCanvas *cGrMean = new TCanvas(); grMean->SetMarkerStyle(2); grMean->Draw("AP"); // We draw the energy-ADC channel graph


    TGraphErrors *grRes = new TGraphErrors(N,x, Res,Ex,ERes); grRes->SetTitle("Resolution");
    grRes->GetXaxis()->SetTitle(xAxisName);     grRes->GetYaxis()->SetTitle("Resolution"); 
    TCanvas *cGrRes = new TCanvas(); grRes->SetMarkerStyle(2); grRes->Draw("AP"); // We draw the energy-ADC channel graph

    while(getline(configFile,line)) if (line[0] == '>') break; string opt1 = extract(line);
    while(getline(configFile,line)) if (line[0] == '>') break; string opt2 = extract(line);
    
    
    while(getline(configFile,line)) if (line[0] == '>') break; TString foutName = extract(line);
    while(getline(configFile,line)) if (line[0] == '>') break; TString option = extract(line);

    TFile *fout = new TFile(foutName, option);

    
    if (opt1 == 'y'){
        TF1 *funCal = new TF1("funCal",pol2,0.0,16384.0,3);
        funCal->FixParameter(0,0);
        if (opt2 == 'y') {
            funCal->FixParameter(2,0);
        }
        grMean->Fit(funCal);
        fout->WriteObject(cGrMean,"MeanGraph");
        fout->WriteObject(cGrRes,"Energy Resolution");
        fout->WriteObject(funCal,"CalibrationFunction");
        fout->Close();
        while(getline(configFile,line)) if (line[0] == '>') break; line = extract(line);
        if (line == 'y'){
            calibrate(finName, foutName, funCal);
        }

    } else {
        fout->WriteObject(cGrMean,"MeanGraph");
        fout->WriteObject(cGrRes,"Energy Resolution");
        fout->Close();
    }
    

}

