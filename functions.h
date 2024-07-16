
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
double fBACK(double *x, double *par) {  // Only the non-gaussian part of the previous function
    double fval; double xval; xval = *x;
    fval =  par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];
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

    TF1 *funBACK = new TF1("funBACK", fBACK,xmin,xmax,7);  // Only the Background
        funBACK->SetParameters(funGEL->GetParameters());
        funBACK->SetLineColor(kGreen);
    funBACK->Draw("SAME");

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

    TF1 *funBACK = new TF1("funBACK", fBACK,xmin,xmax,7);  // Only the Background
        funBACK->SetParameters(funGEL2->GetParameters());
        funBACK->SetLineColor(kGreen);
        funBACK->Draw("SAME");
    return funGEL2;
}
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
void calibrate(TString inputfile, TString outputfile, TF1 *calFun){// Given an imput file with histograms, an output file (which will be overwritten), and a function, recalibrates ALL histograms in imput and puts them in the output.
    // TString inputfile, TString outputfile, TF1 calFun
    // We open the corresponding files
    TFile *fin= new TFile(inputfile);
    TFile *fout= new TFile(outputfile, "UPDATE");

    for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects
        TKey *key = static_cast<TKey*>(k);
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;  // Botam aquesta passa si no és histograma
        TH1 *hist = key->ReadObject<TH1>();

        // We now have the hist "hist"
        TString oldTitle = hist->GetName();
        TString newTitle = oldTitle + "Cal";
        int N = hist->GetNbinsX();
        TH1D *newhist = new TH1D(newTitle,newTitle,N,0,calFun->Eval(N));

        double xlow[N+1];
        for(int ii = 0; ii < N; ii++){
            xlow[ii] = calFun->Eval(ii-0.5);
            newhist->SetBinContent(ii,hist->GetBinContent(ii));
        }
        xlow[N]=calFun->Eval(N+0.5);
        newhist->GetXaxis()->Set(N,xlow);

        fout->WriteObject(newhist, newTitle);
   }
   fin->Close(); fout->Close();
}

class Peak{
private:
    /* data */
public:
    Peak(TH1D *hist, double xmin, double xmid, double xmax, char erc = '1', int numb = 0);
    Peak(TH1D *hist, double xmin, double xmax, char erc = '1');
    TF1 PeakFunction;
    TF1 BackFunction;
    double Mean; double EMean;
    double Sigma, double ESigma;
    double Amplitude; double EAmplitude;
    double FWHM; double EFWHM;
    double Res; double ERes;
    double Area; double EArea;
}

Peak::Peak(TH1D *hist, double xmin, double xmid, double xmax, char erc = '1', int numb = 0){
    TF1 fitDouble = fitGEL2(hist, xmin, xmid, xmax, erc);
    PeakFunction = new TF1("fun",GEL, xmin, xmax, 7);
    BackFunction = new TF1("backFun", BACK,xmin,xmax,7);
    double par[10] = fitDouble->GetParameters();
    if(numb = 0){
        PeakFunction->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
        Mean = fitDouble->GetParameter(0); EMean = fitDouble->GetParError(0);
        Sigma = fitDouble->GetParameter(1); ESigma = fitDouble->GetParError(1);
        Amplitude = fitDouble->GetParameter(2); EAmplitude = fitDouble->GetParError(2);
    }else {if(numb = 1){
        PeakFunction->SetParameters(par[7],par[8],par[9],par[3],par[4],par[5],par[6]);
        Mean = fitDouble->GetParameter(7); EMean = fitDouble->GetParError(7);
        Sigma = fitDouble->GetParameter(8); ESigma = fitDouble->GetParError(8);
        Amplitude = fitDouble->GetParameter(9); EAmplitude = fitDouble->GetParError(9);
    }}
    BackFunction->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
    FWGM = Sigma*2.355; EFWHM = ESigma * 2.355; 
    ERes = Res * (ESigma/Sigma + EMean/Mean);
    Area = Amplitude * Sigma * TMath::sqrt(2*TMath::Pi()); EArea = Area * (ESigma/Sigma + EAmplitude/Amplitude);
}

Peak::Peak(TH1D *hist, double xmin, double xmax, char erc = '1'){
    PeakFunction = fitGEL(hist, xmin, xmax, erc);
    BackFunction = new TF1("backFun", BACK,xmin,xmax,7);
    Mean = PeakFunction->GetParameter(0); EMean = PeakFunction->GetParError(0);
    Sigma = PeakFunction->GetParameter(1); ESigma = PeakFunction->GetParError(1);
    Amplitude = PeakFunction->GetParameter(2); EAmplitude = PeakFunction->GetParError(2);
    BackFunction->SetParameters(PeakFunction->GetParameters());
    FWGM = Sigma*2.355; EFWHM = ESigma * 2.355; 
    ERes = Res * (ESigma/Sigma + EMean/Mean);
    Area = Amplitude * Sigma * TMath::sqrt(2*TMath::Pi()); EArea = Area * (ESigma/Sigma + EAmplitude/Amplitude);
}

