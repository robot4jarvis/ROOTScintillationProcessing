#define N 10
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


double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;

    // The function is:
    // fval = gauss(x, p0,p1) * p2 + p3*x + p4 * Erfx (x-p0)*p5 + p6
}
double GEL2(double *x, double *par) { // double gaussian + background. The parameters of the extra gaussian are 7,8 and 9 (10 in total)
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false)*par[2]  + TMath::Gaus(xval, par[7],par[8],false) * par[9] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;
}
double pol2(double *x, double *par) {  // function y = a*x*x + b*x + c
    double fval; double xval; xval = *x;
    fval = xval *xval* par[0] + xval* par[1] + par[2];
    return fval;
}
TF1 *fitGEL(TH1D *hist, double xmin, double xmax, char erc = '1'){

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
        if(erc == '0') {funGEL->FixParameter(4,0);};
    hist->Fit("funGEL","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL;
}
TF1 *fitGEL2(TH1D *hist, double xmin, double xmid, double xmax, char erc = '1'){

    // we open the file with the Histograms and retreive the second one
    //TFile *f = new TFile("Histograms.root");
    //TH1D *hist = (TH1D*)f->Get(histName);

    // we draw it (just to see that it worked)
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL2 = new TF1("funGEL2", GEL2, xmin, xmax, 10); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL2->SetParameters((xmid+xmin)/2, (xmid-xmin)/2, ymax, 0,0,1,0,(xmax+xmid)/2, (xmax-xmid)/2,ymax); // We set initial parameters (just a guess)
        funGEL2->SetParNames("Centroid1","Sigma1","PeakVal1", "BackgroundSlope","ErcHeight","ErcMult","Baseline","Centroid2","Sigma2","PeakVal2");   // We give the parameters a proper name (lineal part ax + b)
        funGEL2->SetParLimits(0,xmin, xmid);            funGEL2->SetParLimits(7,xmid, xmax);
        funGEL2->SetParLimits(1,0,(xmax-xmin));         funGEL2->SetParLimits(8,0,(xmax-xmin));
        funGEL2->SetParLimits(2,0,1.5*ymax);            funGEL2->SetParLimits(9,0,1.5*ymax);
        funGEL2->SetParLimits(4,0,1.5*ymax);

        if(erc == '0') {funGEL2->FixParameter(4,0);}
    hist->Fit("funGEL2","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL2;
}

void macro(string configFileName = "settings.example"){
    ifstream configFile;
    configFile.open(configFileName, ios::in);
    if(configFile.is_open()) cout<<"Reading configuration file.\n";
    string line;  while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the first interesting line

    getline(configFile,line); const int N = stoi(line.substr(line.find("[")+1,line.find("]")-line.find("[")-1)); // we get N

    while(getline(configFile,line)) if (line[0] == '>') break; // We skip until the next interesting line
    double_t x[N]; double_t Ex[N]; double_t mean[N]; double_t Emean[N]; double_t sigma[N]; double_t Esigma[N];
    double_t Res[N]; double_t ERes[N];

    TString histName; TString histTitle; string set; double_t xmin; double_t xmid; double_t xmax;
    int i = 0;
    while (i<N){
        getline(configFile,line);
        if (line[0] != '>') break; // If we finished reading those values, we exit.
        histName = line.substr(line.find('"')+1,line.find('"',line.find('"')+1)-line.find('"')-1); line.erase(0,line.find('"',line.find('"')+1)+1);
        histTitle = line.substr(line.find('"')+1,line.find('"',line.find('"')+1)-line.find('"')-1); line.erase(0,line.find('"',line.find('"')+1)+1);
        
        TH1D *hist = (TH1D*)fin->Get(histName);  hist->SetTitle(histTitle);  // We fetch the histogram from the file
        std::cout<<"     Reading histogram "<<histTitle<<"\n";
        while (true){
            if(line.find("[") == -1) break;
            set = line.substr(line.find("[")+1,line.find("]")-line.find("[")-1); line.erase(0,line.find("]")+1);
            //std::cout<<"      >Performing fit on peak: "<<set<<"\n";
            char opt = set[1];
            if (set[0] == 'G'){
                xmin = stod(set.substr(set.find(",")+1,set.find(",",set.find(",")+1))); set.erase(0,set.find(",",set.find(",")+1));
                xmax = stod(set.substr(set.find(",")+1,set.find(",",set.find(",")+1)));
                TF1 *funfit = fitGEL(hist, xmin, xmax);
                    mean[i] = funfit->GetParameter(0); Emean[i] = funfit->GetParError(0);
                    sigma[i] = funfit->GetParameter(1); Esigma[i] = funfit->GetParError(1);
                    i++;

            } else if(set[0] == 'D'){
                xmin = stod(set.substr(set.find(",")+1,set.find(",",set.find(",")+1))); set.erase(0,set.find(",",set.find(",")+1));
                xmid = stod(set.substr(set.find(",")+1,set.find(",",set.find(",")+1))); set.erase(0,set.find(",",set.find(",")+1));
                xmax = stod(set.substr(set.find(",")+1,set.find(",",set.find(",")+1))); set.erase(0,set.find(",",set.find(",")+1));
                
                TF1 *funfit = fitGEL2(hist, xmin, xmid, xmax);
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
    TString xAxisName = line.substr(line.find("[")+2,line.find("]")-line.find("[")-3);


    while(getline(configFile,line)) if (line[0] == '>') break; string xString = line.substr(line.find("[")+1,line.find("]")-line.find("[")-2);
    string val; char c; int j = 0;
    for(int ii = 0; ii < xString.length(); ii ++){
        //reads x string
        c = xString[ii];
        if(c == ' ') continue;
        if((c == ',') || (ii == xString.length()-1)){
            x[j] = stod(val);
            j++;
            val = "";
        }
        else val = val + c;
    }
    

    while(getline(configFile,line)) if (line[0] == '>') break; xString = line.substr(line.find("[")+1,line.find("]")-line.find("[")-2);
    val = ""; j = 0;
    for(int ii = 0; ii < xString.length(); ii ++){
        //reads x string
        c = xString[ii];
        if(c == ' ') continue;
        if((c == ',')||(ii == xString.length()-1)){
            Ex[j] = stod(val);
            j++;
            val = "";
        }
        else val = val + c;
    }

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
}

