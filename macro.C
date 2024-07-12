#include <iostream>
#include <fstream>
#include <string>
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
using namespace std;

double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;

    // The function is:
    // fval = gauss(x, p0,p1) * p2 + p3*x + p4 * Erfx (x-p0)*p5 + p6
}

// This function, given a histname (that must be inside a Histograms.root file), a xmin and a xmax (ROI), fits a gaussian function (+background) and returns the fitted function
TF1 *fitGEL(TH1D *hist, double xmin, double xmax){

    // we open the file with the Histograms and retreive the second one
    //TFile *f = new TFile("Histograms.root");
    //TH1D *hist = (TH1D*)f->Get(histName);

    // we draw it (just to see that it worked)
    hist->GetXaxis()->SetRange(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL = new TF1("funGEL", GEL, xmin, xmax, 7); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL->SetParameters((xmax+xmin)/2, (xmax-xmin)/2, ymax/2, 0,0,1,0); // We set initial parameters (just a guess)
        funGEL->SetParNames("Centroid","Sigma","PeakVal", "BackgroundSlope","ErcHeight","ErcMult","Baseline");   // We give the parameters a proper name (lineal part ax + b)
        funGEL->SetParLimits(0,xmin, xmax);
        funGEL->SetParLimits(1,0,(xmax-xmin));
        funGEL->SetParLimits(2,0,1.5*ymax);
        funGEL->SetParLimits(4,0,1.5*ymax);
        funGEL->FixParameter(4,0);   funGEL->FixParameter(5,0); 


    hist->Fit("funGEL","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL;
}

void macro(string configFileName = "settings.example"){
        ifstream configFile;
    configFile.open(configFileName, ios::in);
    if(configFile.is_open()) cout<<"Reading configuration file.\n";
    string line;  while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the first interesting line


    TString finName = line.substr(line.find("[")+2,line.find("]")-line.find("[")-3); // We obtain the imput filename
    TFile *fin = new TFile(finName); std::cout<<"Opened input root file: "<<finName<<"\n"; // We open the root filename

    getline(configFile,line); const int N = stoi(line.substr(line.find("[")+1,line.find("]")-line.find("[")-1)); // we get N

    while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the next interesting line
    double_t x[N]; double_t Ex[N]; double_t mean[N]; double_t Emean[N]; double_t res[N]; double_t Eres[N];

    


}
