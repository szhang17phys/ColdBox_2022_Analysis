#include "MYCODES.h"


class  MyFunctionObject{
public:

    Int_t n_peaks;
    // use constructor to customize your function object
    Double_t operator()(Double_t *x, Double_t *par) {
        Double_t f;
        Double_t xx = x[0];
        f  = par[0]*exp(-0.5*TMath::Power((xx-par[1])/par[2],2)); // first argument
        f = f+par[3]*exp(-0.5*TMath::Power((xx-par[4])/par[5],2));
        f = f+par[6]*exp(-0.5*TMath::Power((xx-par[7])/(TMath::Power((2),0.5)*par[5]),2));
        for(Int_t i = 1; i<n_peaks; i++){
            f = f+ abs(par[i+7])*exp(-0.5*TMath::Power((xx-(par[4]+(i+1)*(par[7]-par[4])))/(TMath::Power((i+2),0.5)*par[5]),2));
        }
        return f;
    }
};


































class Calibration

{
    
    
    
public:
    
  
  Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
  Int_t nbits = 14; // DIGITIZER bits
  Int_t linhasEvento = 9000;
  
  Int_t channel = 1;
  
  // _______________ Parameters for giveMeSphe_MeanPulseVersion _______________/
  Double_t time_sample = 6000; // in ns
  
  Double_t timeLow = 3650; // integration limit
  Double_t timeHigh = 4200;
  
  Int_t fn = 0;
  
  
  // _______________ Parameters for fit_sphe_wave _______________/
  
  string rootFile;
  
  Int_t n_peaks = 7;
  Double_t peak0 = 0.1;
  Double_t mean0 = -900;
  Double_t sigma0 = 500;
  
  Double_t peak1 = 0.004;
  Double_t mean1 = 5000;
  Double_t sigma1 = 400;
  
  Double_t startpoint = 0.001;
  
  Double_t xmin = -10000;
  Double_t xmax = 40000;

  Double_t deltaplus=1.4;
  Double_t deltaminus=1.2;
  
  Int_t rebin = 4;
  
  Bool_t fixZero = false;
  
  
  Double_t sphe_charge = 0; // wave0
  Double_t sphe_charge2 = 0; // wave0
  
  
  // This was add here to try fitting dark noise data with already found values
  Bool_t darknoise = false;
  
  Bool_t is_poisson_test = false; // if running tests of poisson statistics
  
  
  
  
  
  // _______________ Parameters for smoothing  _______________/
  Int_t smooth_factor = 35;


// ____________________________________________________________________________________________________ //
void fit_sphe_wave(string name){
   
    makeSphe(name);
}




// ____________________________________________________________________________________________________ //
string startingPump(){
    string f = "gaus(0) + gaus(3)";
    Int_t aux = 0;
    for(Int_t i = 0; i<n_peaks; i++){
        f = f + " + gaus(" + to_string(i+6+aux) + ")";
        aux = aux+2;
    }
    return f;
//         TF1 *func = new TF1("func","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",-2000,5000);

}

// ____________________________________________________________________________________________________ //
void getMyParameters(Double_t peaks[],Double_t stdpeaks[],TF1 *func){
    for(Int_t i=0; i<n_peaks; i++){
        peaks[i] = (i+2)*(func->GetParameter(4));
        stdpeaks[i] = TMath::Power((i+2),0.5)*func->GetParameter(5);
//         cout << peaks[i] << endl;
    }
}

// ____________________________________________________________________________________________________ //
void makeSphe(string histogram){
    
    string histogram_tempo = histogram+"_stat";
    TFile *f1 = new TFile(rootFile.c_str(),"READ");
    TH1D *hcharge = (TH1D*)f1->Get(histogram.c_str());
    
    hcharge->Rebin(rebin);
    
    ofstream out;
    out.open("sphe.txt", ios::app);

    // ____________________________ Start of sphe fit ____________________________ //
    hcharge->GetYaxis()->SetTitle("Normalized count");
    hcharge->GetYaxis()->SetTitleOffset(1.0);
    hcharge->GetXaxis()->SetTitle("Charge (ADC*nsec)");
    
    
    TCanvas *c1 = new TCanvas("c1","Carga");
    // c1->SetLogy();
    gPad->SetGrid(1,1);
    gPad->SetTicks(1,1);
    gStyle->SetOptFit();
    
    //First function, will almost fit freely
    TF1 *func = new TF1("func",startingPump().c_str(),xmin,xmax);
    TF1 *fu[2+n_peaks];
    string funame;
    for(Int_t i = 0; i<(2+n_peaks); i++){
        funame = "fu_"+to_string(i);
        fu[i] = new TF1(funame.c_str(),"gaus(0)",xmin,xmax);
    }
    
    Double_t peaks[n_peaks];
    Double_t stdpeaks[n_peaks];
    
    if(peak0==0){
        fixZero = true;
    }
    
    func->SetParameters(peak0,mean0,sigma0,peak1,mean1,sigma1); // this values can change    
    
    Int_t aux = 0;
    for(Int_t i = 0; i<n_peaks; i++){
        func->SetParameter((i+6+aux),startpoint);
        aux++;
        func->SetParameter((i+6+aux),(i+2)*mean1);
        aux++;
        func->SetParameter((i+6+aux),(i+2)*sigma1);
        startpoint = startpoint/5;
    }
    func->SetParName(4,"#mu");
    func->SetParName(5,"#sigma");
    
    func->SetNpx(1000);
    
    if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
    }
    
    if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
    }
    
    getMyParameters(peaks,stdpeaks,func);
    
    
    Double_t scale = 1/(hcharge->Integral());
    hcharge->Scale(scale);
    hcharge->Draw("hist");
    hcharge->Fit("func","R0Q");
    
    // recovery parameters
    getMyParameters(peaks,stdpeaks,func);
    aux=0;
    // use as fixed
    for(Int_t i = 0; i<n_peaks; i++){
        
        func->FixParameter((i+7+aux),peaks[i]-(i+1)*func->GetParameter(1));
        aux++;
        func->FixParameter((i+7+aux),stdpeaks[i]);
        aux++;
    }
        
    if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
    }
    
    if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
    }
    hcharge->Fit("func","R0Q");
    
    // set a new function, now fixed for real
    
    MyFunctionObject MyFunc;
    MyFunc.n_peaks = n_peaks;
    
    TF1 *lastOne = new TF1("lastOne",MyFunc,xmin,xmax,7+n_peaks);
    
    
    lastOne->SetParameter(0,func->GetParameter(0));
    lastOne->SetParameter(1,func->GetParameter(1));
    lastOne->SetParameter(2,func->GetParameter(2));
    
    lastOne->SetParameter(3,func->GetParameter(3));
    lastOne->SetParameter(4,func->GetParameter(4));
    lastOne->SetParameter(5,func->GetParameter(5));
    
    lastOne->SetParameter(6,func->GetParameter(6));
    lastOne->SetParameter(7,func->GetParameter(7));
    aux = 0;
    for(Int_t i = 1; i<n_peaks; i++){
        lastOne->SetParameter((i+8-1),func->GetParameter(i+8+aux));
        lastOne->SetParName((i+8-1),Form("A_{%d}",i+2));
        aux=aux+2;
    }

    
    lastOne->SetParName(0,"A_{baseline}");
    lastOne->SetParName(1,"#mu_{baseline}");
    lastOne->SetParName(2,"#sigma_{baseline}");
    lastOne->SetParName(3,"A_{1}");
    lastOne->SetParName(4,"#mu_{1}");
    lastOne->SetParName(5,"#sigma_{1}");
    lastOne->SetParName(6,"A_{2}");
    lastOne->SetParName(7,"#mu_{2}");
    // lastOne->SetParName(8,"#sigma_{2}");
    
    if(darknoise){
        lastOne->FixParameter(4,sphe_charge);
        lastOne->FixParameter(7,sphe_charge2);
    }
    if(fixZero){
        lastOne->FixParameter(0,0);
        lastOne->FixParameter(1,0);
        lastOne->FixParameter(2,1);
    }
    
    hcharge->Fit("lastOne","R");
    
    Double_t stddev = lastOne->GetParameter(2);
    
    fu[0]->SetParameter(0,lastOne->GetParameter(0));
    fu[0]->SetParameter(1,lastOne->GetParameter(1));
    fu[0]->SetParameter(2,lastOne->GetParameter(2));
    
    fu[1]->SetParameter(0,lastOne->GetParameter(3));
    fu[1]->SetParameter(1,lastOne->GetParameter(4));
    fu[1]->SetParameter(2,lastOne->GetParameter(5));
    
    fu[2]->SetParameter(0,lastOne->GetParameter(6));
    fu[2]->SetParameter(1,lastOne->GetParameter(7));
    fu[2]->SetParameter(2,(TMath::Power((2),0.5)*lastOne->GetParameter(5)));
    
   for(Int_t i = 1; i<n_peaks; i++){
       fu[i+2]->SetParameter(0,abs(lastOne->GetParameter(i+7)));
       fu[i+2]->SetParameter(1,(lastOne->GetParameter(4) + (i+1)*(lastOne->GetParameter(7)-lastOne->GetParameter(4))));
       fu[i+2]->SetParameter(2,(TMath::Power((i+2),0.5)*lastOne->GetParameter(5)));
    }
    
    lastOne->SetNpx(1000);
    
    hcharge->Draw("hist");

    
    hcharge->GetXaxis()->SetRangeUser(-5000,200000);
    hcharge->StatOverflows(kTRUE);
    lastOne->SetRange(-1000,xmax);
    
    
    for(Int_t i = 0; i<(2+n_peaks); i++){
        fu[i]->SetLineColor(kGray+1);
        fu[i]->SetNpx(1000);
        fu[i]->Draw("SAME");
    }
    lastOne->Draw("LP SAME");
    string name = histogram + ".root";
//     c1->Print(name.c_str());
    cout << "1th peak = " << lastOne->GetParameter(4) << endl;
    cout << "2th peak = " << lastOne->GetParameter(7) << endl;
    cout << "sphe charge = " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
    cout << " SNR = " << lastOne->GetParameter(4)/sqrt(pow(lastOne->GetParameter(2),2)+pow(lastOne->GetParameter(5),2)) << endl;
    out <<  lastOne->GetParameter(4) << " " << lastOne->GetParameter(7) << " " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
    
    // ____________________________ Finish of sphe fit ____________________________ //
    

    //_________________ Drawing lines for the cross-talk probability _________________ //
    
    sphe_charge = lastOne->GetParameter(4);
    sphe_charge2 = lastOne->GetParameter(7);
    
    Double_t delta1 = (sphe_charge2 - sphe_charge)/deltaminus;
    Double_t delta2 = deltaplus*(sphe_charge2 - sphe_charge);
//     Double_t delta2 = sphe_charge+(sphe_charge2 - sphe_charge)/2;
    
    Double_t ymax = hcharge->GetMaximum();

    TLine *l1 = new TLine(delta1,0,delta1,ymax);
    TLine *l2 = new TLine(delta2,0,delta2,ymax);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    
    l1->Draw("");
    l2->Draw("");
    
    
    if(darknoise){
        // ____________________________ Start dark noise CT analysis ____________________________ //
        
        cout << "\n\n\nMaking ratio between 1 sphe and 2 sphe: " << endl;
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram 
        // could be taken, but the effort 
        // is bigger.
        // Also, remember to multiply by 
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);
        
        cout << "Total 1 = " << oneIntegral << " \t total 2 =  " << twoIntegral << endl;
        cout << "Ratio = " << twoIntegral/oneIntegral << endl;
        
        cout << "Total events normalized = " << lastOne->Integral(xmin,xmax) << endl;
        
        cout << "CT probability = " << twoIntegral/(oneIntegral+twoIntegral) << endl;
        
        
        
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        
        cout << "\n\n\n Another possibility would be CT = " << (twoIntegral+threeIntegral)/(oneIntegral+twoIntegral+threeIntegral) << endl;
        
        // ____________________________ Finish Poisson analysis ____________________________ //
        
    }
    
    
    
    if(is_poisson_test){
        
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram 
        // could be taken, but the effort 
        // is bigger.
        // Also, remember to multiply by 
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);       
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        

        Double_t lambda = -TMath::Log(zeroIntegral/totalIntegral);
        
        cout << "Lambda = " << lambda << endl;
        
        Double_t ct = hcharge->GetMean()/(lambda*oneMean);
        
        cout << "Cross-talk = " << ct << endl;
        
        
        
//         cout << "teste -> " << hcharge->GetMean() << endl;
        // ____________________________ Finish Poisson analysis ____________________________ //
    }
}

    
};









// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //





    




















    





























































class SPHE{
       
    
    
    
public:
    
    
TFile *fout;

TTree *tout;

Bool_t creation = true; //to verify creation of ttree branches


TFile *fwvf;
TTree *twvf;


Double_t value = 0;
Double_t desv = 0;
Int_t channel = 1;


// ____________________ Variables to calculate ____________________ //

TH1D *hbase = new TH1D("hbase","histogram for baseline",5*800,-400,400);
TH1D *hbase_smooth = new TH1D("hbase_smooth","histogram for baseline smoothed",5*800,-400,400);
// TH1D *hcharge = new TH1D("hcharge","",100000,-50000,50000);
  TH1D *hcharge = new TH1D("hcharge","",50000,0,0);
TH1D *hzero = new TH1D("hzero","",120000,-200000,2*1300000);
TH1D *hnobase = new TH1D("hnobase","",120000,-200000,2*1300000);

TH1D *hstat = new TH1D("hstat","",120000,-200000,2*1300000);


Double_t dtime = 2.;

Double_t mean = 0;
Double_t stddev = 0;
Double_t tolerance; // n sigmas
Double_t baseLimit;

Double_t timeLimit; // time after LED signal
Double_t timeLow ; // integration time before peak
Double_t timeHigh ; // integration time after peak

Double_t start = 0;
Double_t finish = memorydepth*dtime;

Double_t noiseLow;
Double_t noiseHigh;

Double_t forBaseline = 2000;
Double_t baseCharge = 0;


Int_t interactions; // for moving avarage
Int_t midpoint = 0; // midpoint that takes the avarage
Int_t width = 0;
Double_t sum=0;

vector<Double_t> peakPosition;
vector<Double_t> peakMax;

Double_t charge = 0;
Double_t treeCharge = 0;
Double_t treePeak = 0;
Double_t strikes = 0;
Double_t ptsLow = 0;
Double_t ptsHigh = 0;
Bool_t discard = false;

TH1D *hfilter = new TH1D("hfilter","",5*800,-400,400);
Double_t devFilter;


Int_t cut = 0;

Int_t gap = 2000; // at the and, we wont look at the final 4000 ns
Double_t lowerThreshold = -9999;
Double_t maxHits = 15;
Double_t nSigmas_integration = 0.5;
Int_t atleast = 120; //at least this time in the first integration and second integration
Int_t lowPtsLimit = 20;
Bool_t just_a_test = false; //if true, will only look 5000 events;
Bool_t led_calibration = false; // for led integration only

Int_t just_this = 5000;


// ____________________ matching filter ____________________

vector<Double_t> mySample;
vector<Double_t> matched;
vector<Double_t> matchingAreas;

Double_t shift = 0;
Double_t matched_value = 0;
Double_t higherValue = 0;
Double_t tempPosition;


Double_t area_off = 0;
Double_t matchTrigger = 1.;

Double_t from = 0;
Double_t to = 0;


Bool_t matching = false;
Bool_t badBaseline = false; // this need to be used in the first Gamma run...

// ____________________ Variables to show events ____________________ //

const static Int_t nshow = 100;
Double_t my_events[nshow];
Double_t eventNow = 0;
Int_t aux_events = 0;


TGraph *g_smooth;
TGraph *g_normal;
TGraph *g_points;
    
vector<Double_t> temp_peak;
vector<Double_t> peak_smooth;
vector<Double_t> timeg;

vector<Double_t> selected_peaks;
vector<Double_t> selected_time;
Int_t teste = 0;


Bool_t fillg = true;

Int_t mark = 0;

Double_t xi = 0; // time interval for sampling (ns)
Double_t xf = 1000;
Double_t center = 300;
Int_t nsample = (xf-xi)/2+1;
Double_t counter = 0;

TH2D *hsample = new TH2D("sample","peak vs time", nsample, xi, xf, 3000,-1500,1500);

vector<Double_t> ymean;
vector<Double_t> ynormal;
vector<Double_t> xmean;

// variables to find mean waveform 
Bool_t get_wave_form = false;
vector<Double_t> mean_waveforms;
Double_t wvf[memorydepth];
Double_t wvfcharge;
Bool_t valid;
Int_t naverages;
Double_t mean_before = 100;
Double_t mean_after = 400;
Int_t low_cut = 800; // to select sphe waveforms, inside this region, peak should not be bigger then val_cut
//it was 920 before
Int_t high_cut = 2000;
Double_t val_cut = 6;
TMultiGraph *gm = new TMultiGraph();
TGraph *gwaveforms;
Double_t shifter = 12;

Double_t sphe_charge;
Double_t sphe_charge2;
Double_t delta;
Double_t deltaplus = 1.4;
Double_t deltaminus = 1.2;
Double_t sphe_std;

Double_t sphe_charge_ch0;
Double_t sphe_charge_ch1;

Double_t sphe_charge2_ch0;
Double_t sphe_charge2_ch1;

Double_t sphe_std_ch0;
Double_t sphe_std_ch1;

string charge_status = "";

Bool_t darkNoise = false;

Double_t baselineTime=2000;
Double_t too_big = 1000;
Double_t waiting = 0;
// ____________________________________________________________________________________________________ //

void giveMeSphe_darkCount(string name){
  gROOT->SetBatch(kTRUE);
  fout = new TFile(Form("sphe_histograms_darkCount_Ch%i.root",channel),"RECREATE");
  tout = new TTree("t1","baseline info");
  darkNoise = true;
  fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
  twvf = new TTree("t1","mean waveforms");
  
  if(led_calibration==false){
    makeHistogram(name);  
    fout->WriteObject(tout,"t1","TObject::kOverwrite");
    fout->Close();
  }
  else{
    makeSimpleHistogram(name);
    // get_wave_form = false;
  }
  if(get_wave_form==false){
    fwvf->Close();
    system(Form("rm sphe_waveforms_Ch%i.root",channel));
  }
  else{
    fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
  }
  channel=channel + 1;
}


// ____________________________________________________________________________________________________ //
void giveMeSphe(string name){
  gROOT->SetBatch(kTRUE);
  
  fout = new TFile("sphe_histograms.root","RECREATE");
  tout = new TTree("t1","baseline info");
  darkNoise = false;
  fwvf = new TFile("sphe_waveforms.root","RECREATE");
  twvf = new TTree("t1","mean waveforms");
  
  makeHistogram(name);

  fout->WriteObject(tout,"t1","TObject::kOverwrite");
  fout->Close();
  if(get_wave_form==false){
    fwvf->Close();
    system(Form("rm sphe_waveforms_Ch%i.root",channel));
  }
  else{
    fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
  }
  channel=channel + 1;
}

// ____________________________________________________________________________________________________ //
void integrateSignal(){
    charge = 0;
    Int_t npeaks = peakPosition.size();
    value = 0;
    desv = 0;

    Int_t discardedPeaks = 0;

    Int_t aux_sample = 0;
    
    
    if(channel==1){
        sphe_charge = sphe_charge_ch0; // wave0
        sphe_charge2 = sphe_charge2_ch0; // wave0
        delta = sphe_charge2 - sphe_charge;
        sphe_std = sphe_std_ch0;
    }
    if(channel==2){
        sphe_charge = sphe_charge_ch1; // wave0
        sphe_charge2 = sphe_charge2_ch1; // wave0
        delta = sphe_charge2 - sphe_charge;
        sphe_std = sphe_std_ch1;
    }

    Double_t statcharge[npeaks];
    Double_t statpeak[npeaks];
    Double_t trackLow[npeaks];
    Double_t trackHigh[npeaks];
    Double_t trackStrikes[npeaks];
    Double_t statpos[npeaks];
    Bool_t discard_this[npeaks];
    
    vector<vector<Double_t>> temp_waveforms(npeaks);
    vector<Double_t> waveforms(memorydepth);
//     vector<TGraph*> gwaveforms(npeaks);
    vector<Bool_t> notAGoodWaveform(npeaks);


    Int_t auxstat = 0;
    Bool_t notGood = false;
    Double_t lastTime = 0 ;
    Double_t peakStd = 0;
    
    Double_t thismanypointsLow = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t thismanypointsHigh = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t thismanypointsBase = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t totalmany = 0; //this helps speacially for matching...
    Bool_t next_is_bad = false;
    Double_t waitingInterval = -1;
    for(Int_t i = 0; i<npeaks; i++){
      statcharge[i] = 0;
      statpeak[i] = 0;
      trackStrikes[i] = 0;
      discard_this[i] = false;
      notGood = false;
      thismanypointsLow = 0;
      thismanypointsHigh = 0;
      thismanypointsBase = 0;
      totalmany = 0;
      
      
      if((peakPosition[i]+2*timeHigh>=peakPosition[i+1] && i+1!=npeaks) || next_is_bad || peakPosition[i]<waitingInterval){

        selected_peaks.push_back(temp_peak.at(peakPosition[i]/dtime));
        selected_time.push_back(peakPosition[i]);
        discard_this[i] = true;
        for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
            if(temp_peak.at(j)>=statpeak[i]){
                statpeak[i] = temp_peak.at(j);
            }
        }
        if(statpeak[i]>too_big){
         
          waitingInterval = peakPosition[i]+waiting;
        }
        if(peakPosition[i]>waitingInterval) next_is_bad = true;
        else next_is_bad = false;
        continue;
      }
      else{
        next_is_bad = false;
      }
      
      
      
      notAGoodWaveform[i] = false;

            
                
      Int_t strikes = 0;
      // integration for the back part
      for(Int_t j = (peakPosition.at(i))/dtime; j>= (peakPosition[i]-timeLow)/dtime; j--){
        charge += temp_peak.at(j);
        statcharge[i]+= temp_peak.at(j);
        if(temp_peak.at(j)>=statpeak[i]){
          statpeak[i] = temp_peak.at(j);
        }
        selected_peaks.push_back(temp_peak.at(j));
        selected_time.push_back(timeg.at(j));
        thismanypointsLow++;
        totalmany++;
        if(temp_peak[j]<lowerThreshold){
          strikes++;
        }
      }
      
      
      // integration for the front part
      for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
        charge += temp_peak.at(j);
        
        statcharge[i]+= temp_peak.at(j);
        if(temp_peak.at(j)>=statpeak[i]){
          statpeak[i] = temp_peak.at(j);
        }
        selected_peaks.push_back(temp_peak.at(j));
        selected_time.push_back(timeg.at(j));
        thismanypointsHigh++;
        totalmany++;
        if(temp_peak[j]<lowerThreshold){
          strikes++;
        }
      }
      
      if(statpeak[i]>too_big){
        waitingInterval = peakPosition[i]+waiting;
        discard_this[i] = true;
      }
      if(strikes>=maxHits){
//         cout << "strikes " << strikes << ".... " << eventNow << endl;
        discard_this[i] = true;
      }
      
      // for get_wave_form 
      if(peakPosition.at(i) - mean_before>=0 && peakPosition.at(i)+mean_after<memorydepth*dtime){
        for(Int_t j = peakPosition.at(i)/dtime - mean_before/dtime; j <= peakPosition.at(i)/dtime+mean_after/dtime; j++){
          temp_waveforms[i].push_back(temp_peak.at(j));
//           if(temp_waveforms[i].size()>200/dtime && temp_waveforms[i].size()<400/dtime){
//             if(temp_peak.at(j)<-20) notAGoodWaveform[i]=true;
//           }
        }
      }
      else{
        notAGoodWaveform[i] = true;
      }
      
      trackHigh[i] = thismanypointsHigh;
      trackLow[i] = thismanypointsLow;
      
      
      
      
      if(snap()){
        
        cout << "charge = " << statcharge[i]*dtime << " at " << peakPosition.at(i) << " with " << thismanypointsLow << " + " << thismanypointsHigh << " or this " << thismanypointsBase << " discard ? " << discard_this[i] << endl;
      }
      
      
      
        
        
        
    }
    for(Int_t i = 0; i<npeaks; i++){
        if(discard_this[i]==false){
            
            statcharge[i] = dtime*statcharge[i];
            charge = dtime*charge;
            hcharge->Fill(statcharge[i]);
            hnobase->Fill(statcharge[i]);
            hstat->Fill(statcharge[i]);
            treeCharge = statcharge[i];
            treePeak = statpeak[i];
            ptsHigh = trackHigh[i];
            ptsLow = trackLow[i];
            charge_status += " charge = " + to_string(statcharge[i]);
            tout->Fill();
//             if (charge<-800){
//               cout << "\n\n------>event = " << eventNow << " " << charge << " " << peakPosition[i] << endl;
//             }
//             if(statcharge[i]>10000 && statcharge[i]<11000 && snap()){
            //                 cout << ".............................. charge = " << statcharge[i] << " at " << peakPosition.at(i) << " event = " << eventNow << endl;
            //             }
            
            if(get_wave_form){
              Double_t newbase = 0;
              if(notAGoodWaveform[i]==false){
                if(statcharge[i]>=delta/deltaminus && statcharge[i]<=delta*deltaplus){
                  naverages++;
                  valid = true;
                  for(Int_t j = low_cut/4; j<high_cut/4; j++){
                    if(temp_waveforms[i][j]>val_cut){
                      naverages--;
                      valid = false;
                      break;
                    }
                  }
                }
                else{
                  valid = false;
                }
                for(Int_t j = 0; j<waveforms.size(); j++){
                  if(j<temp_waveforms[i].size()){
                    
                    waveforms[j] = temp_waveforms[i][j];
                    wvf[j]=waveforms[j];
                    if(valid) mean_waveforms[j]+=waveforms[j];
                  }
                  else{
                    waveforms[j] = 0;
                    wvf[j]=0;
                    mean_waveforms[j]+=waveforms[j];
                  }
                  
                }
                wvfcharge = statcharge[i];
                twvf->Fill();
                
                //gwaveforms[i] = new TGraph(waveforms.size(),&timeg[0],&waveforms[0]);
                //gm->Add(gwaveforms[i],"LP");
              }
            }
            
            
        }
        else{
          //             charge_status += " Discarded - > " + to_string(statcharge[i]*dtime) + " ";
//             charge_status += " strikes - > " + to_string(trackStrikes[i]) + " ";
        }
    }
    
    

    

 
}

// ____________________________________________________________________________________________________ //
void searchForPeaks(){
    // For each peak found, we i am looking for the maximum above tolerance
    Int_t n = peak_smooth.size();
    Int_t npeaks = 0;
   
    higherValue = 0;
    tempPosition = 0;
    Bool_t rebase = false;
    Bool_t normalbase = true; // make it true to calculate
    Double_t gapstart = baselineTime;
    Double_t gapend = 6000;
    Double_t actual_finish = finish;
//     for(Int_t i = gapstart/dtime; i<=gapend/dtime; i++){
//       if(temp_peak[i]>200){
//         actual_finish = gapstart;
//       }
//     }

    for(Int_t i = 0; i<n; i++){
      
      if(i>timeLow/dtime && i<(n-timeHigh/dtime) && (i>=start/dtime && i<=actual_finish/dtime)){ //searching out of the beginning
        // if(i>baselineTime/dtime && temp_peak[i]>300){ // not sure why this was here
        if(i>baselineTime/dtime){
            break;
        }
        if(peak_smooth[i]>(mean+tolerance*stddev) && peak_smooth[i-1]<(mean+tolerance*stddev) && (i<=gapstart || i>=gapend)){
          npeaks++;
          
          peakMax.push_back(peak_smooth.at(i));
          peakPosition.push_back(timeg.at(i));

          if(snap()){
            cout << "npeaks = " << npeaks << " " << peakPosition.at(npeaks-1) << " " << peakMax.at(npeaks-1) << " " << eventNow << endl;
            
          }
//           if(i<=baselineTime/dtime){
//               rebase = true;
//           }
//           if(rebase && normalbase){
//             rebase =false;
//             normalbase=false;
//             Int_t binmax = hbase->GetMaximumBin();
//             Double_t newB1 = hbase->GetXaxis()->GetBinCenter(binmax);
//             Double_t newB2 = hbase->GetMean();
//             Double_t newB = (newB1<newB2)? newB1 : newB2;
//             
//             for(Int_t k = 0; k<temp_peak.size();k++){
//               temp_peak[k] = temp_peak[k]-newB;
// //               mean = newB;
//               if(eventNow==my_events[nshow-1]){
// //                 cout << eventNow << " " << k << " " << temp_peak[k] << " " << newB << " " << hbase->GetMean() << endl;
//               }
//             }
//           }
          
        }
        
        
  
      }
      else if(i>finish/dtime){
        break;
      }
    }
    
}

vector<double> delay_line(vector<double> v, double delay_time){
    if(delay_time==0) return v;
    vector<double> res(v.size());
    for(int i=0; i<v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
    }
    return res;
}

// ____________________________________________________________________________________________________ //
void lookWaveform(){

    vector<Double_t> shifted=delay_line(temp_peak, shifter);//cusp(MIN[i], h);
    smoothWithMovingAvarage(shifted); // first calculate avarage

 
    searchForPeaks(); //search for peaks with the moving avarage
    if(snap()){ // if the events match, create the graphs
      
      g_smooth = new TGraph(timeg.size(),&timeg[0],&peak_smooth[0]);
      g_normal = new TGraph(timeg.size(),&timeg[0],&temp_peak[0]);
      //         cout << mean << " " << stddev << endl;
    }
    integrateSignal();
    
    if(snap()){drawMySamples();} // draw sample graphs
}

// ____________________________________________________________________________________________________ //
void makeHistogram(string filename){
    
    if(matching == true){
        getMySample();
        
        if(shift==0){
            shift = (to-from)/2; //make sure to not use "2." we want an integer here
        }
        if(area_off==0){area_off = to - from;}
        
        cout << area_off << " " << shift << endl;
    }
    
    if(creation){
        tout->Branch("value",&value,"value/D");
        tout->Branch("desv",&desv,"desv/D");
        tout->Branch("channel",&channel,"channel/D");
        tout->Branch("treeCharge",&treeCharge,"treeCharge/D");
        tout->Branch("treePeak",&treePeak,"treePeak/D");
        tout->Branch("ptsLow",&ptsLow,"ptsLow/D");
        tout->Branch("ptsHigh",&ptsHigh,"ptsHigh/D");
        tout->Branch("strikes",&strikes,"strikes/D");
//         creation = false;
        
        twvf->Branch("charge",&wvfcharge,"charge/D");
        twvf->Branch("wvf",&wvf,Form("wvf[%i]/D",memorydepth));
        twvf->Branch("valid",&valid,"valid/O");
    }
    
    cout << "reading: " << filename << endl;
    string rootfile = filename + ".root";
    
    TFile *f1 = new TFile(rootfile.c_str(),"READ");
    TTree *t1 = (TTree*)f1->Get("t1");
    ADC_DATA ch;
    TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
    bch->SetAddress(&ch);
    Int_t nentries = t1->GetEntries();
    
    TH1D *hdark = new TH1D("hdark","",1000,-10000,30000);

    //_____________________________ Start creating random data to check _____________________________ //
    TRandom *rmd = new TRandom();

    my_events[0] = 1;
    for(Int_t i = 1; i<nshow; i++){
//         my_events[i] = static_cast<Int_t>(rmd->Uniform(my_events[i-1]+1,my_events[i-1]+50));
        my_events[i] = i+1;
    }
    sort(my_events,my_events+nshow);
    


    eventNow = 0;
    aux_events = 0;
    //_____________________________ End creating random data to check _____________________________ //
    
    
    
    // ____________________ Start resetting globals ____________________ //    
    mean = 0;
    stddev = 0;
    midpoint = 0; // midpoint that takes the avarage
    width = 0;
    sum=0;
    charge = 0;
    discard = false;
    eventNow = 0;
    aux_events = 0;
    teste = 0;
    fillg = true;
    hsample->Reset();
    counter = 0;
    // ____________________ Finish resetting globals ____________________ //
    
    

    
    Double_t aux = 0;
    Int_t j = 0;
    
    temp_peak.clear();
    timeg.clear();
    peak_smooth.clear();
    peakPosition.clear();
    peakMax.clear();
    selected_peaks.clear();
    selected_time.clear();
    
    
    matched.clear();
    
    ymean.resize(nsample);
    ynormal.resize(nsample);
    xmean.resize(nsample);
    for(Int_t i = 0; i<ymean.size(); i++){
        ymean.at(i) = 0;
        xmean.at(i) = 0;
    }

    
    hbase->Reset();
    hstat->Reset();
    hbase_smooth->Reset();
    hcharge->Reset();
    hzero->Reset();
    hnobase->Reset();
    
    mean_waveforms.clear();
    mean_waveforms.resize(memorydepth,0);
    naverages = 0;
    
    for(Int_t i = 0; i<memorydepth; i++){
        timeg.push_back(dtime*i);
    }
    
    if(just_a_test){nentries = just_this;}
    //     aux = 0;
    for(Int_t i = 0; i<nentries; i++){
      bch->GetEvent(i);
      
      for(Int_t i = 0; i<memorydepth; i++){
        temp_peak.push_back(ch.wvf[i]);
        
        if((i<=5000/dtime) && abs(ch.wvf[i])<6){
          hbase->Fill(ch.wvf[i]);
        }
      }
      
      
      aux = ch.event;
      if(static_cast<Int_t>(eventNow)%200==0){
        cout << eventNow << "\r" << flush;
      }
      eventNow =  ch.event;
      // here is were all the calculation is done !!!
      lookWaveform();
      
      if(snap() && aux_events+1<nshow){
        aux_events = aux_events+1;
      }
      
      hbase->Reset();
      hbase_smooth->Reset();
      
      charge_status = "";
      temp_peak.clear();
      peak_smooth.clear();
      peakPosition.clear();
      peakMax.clear();
      selected_peaks.clear();
      selected_time.clear();
      baseCharge = 0;
      matched.clear();
      
    }
    
    
    cout << "____________________________ " << counter/nsample << " \n \n \n" << endl;
    
    if(get_wave_form){
        for(Int_t i = 0; i<memorydepth; i++){
//             cout << i << " " << mean_waveforms[i] << " ";
             
            mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
//             cout << mean_waveforms[i] << endl;
            
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;
        
    }
    
    TCanvas *cwvf = new TCanvas();
    cwvf->cd();
//     gm->Draw("A");
    
    TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);
    
    if(get_wave_form){
//         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
    }

    fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

    f1->Close();

}

// ____________________________________________________________________________________________________ //
void smoothWithMovingAvarage(vector<Double_t> shifted){
    
    Int_t n = shifted.size();
    if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
    }
    else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
    }
    
    
    for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
            peak_smooth.push_back(shifted.at(i));
        }
        else if(i>cut/dtime){
            for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
                sum = sum+shifted.at(j);
//                 cout << sum << endl;
            }
            peak_smooth.push_back(sum/width);
            
            if(timeg.at(i)<=baselineTime && abs(peak_smooth.at(i))<baseLimit){
                hbase_smooth->Fill(peak_smooth.at(i));
            }

        }
        else{
            peak_smooth.push_back(0);
        }

        
        
        sum=0;
    }
    mean = hbase_smooth->GetMean();
    stddev = hbase_smooth->GetStdDev();
    
    
    
}
// ____________________________________________________________________________________________________ //
void drawMySamples(){
    
    TCanvas *c1 = new TCanvas();
    c1->cd(1);
    g_smooth->SetLineColor(kRed);
    g_smooth->SetLineWidth(3);
    
    g_normal->SetLineColor(kBlue);
    g_normal->SetTitle(charge_status.c_str());
    
    g_normal->Draw("AL");
    g_smooth->Draw("L SAME");
    
    TLine *lmean = new TLine(timeLimit,mean,memorydepth*dtime,mean);
    TLine *ldev = new TLine(timeLimit,mean+tolerance*stddev,memorydepth*dtime,mean+tolerance*stddev);
    
    lmean->SetLineColor(kGreen);
    ldev->SetLineColor(kGreen);
    
    lmean->SetLineWidth(2);
    ldev->SetLineWidth(2);
    
    lmean->Draw();
    ldev->Draw();
    Int_t n = selected_peaks.size();
    if(n!=0){
        g_points = new TGraph(n,&selected_time[0],&selected_peaks[0]);
        g_points->SetMarkerColor(kBlack);
        g_points->Draw("P* SAME");
    }
    
    
    
    
    // Some fancy drawing now; 
    if(matchingAreas.size()!=0){
        Int_t nAreas = matchingAreas.size(); //always two points per area of course
        TLine *larea;
        Double_t max = g_normal->GetYaxis()->GetXmax();
        for(Int_t i = 0; i<nAreas; i++){
            larea = new TLine(matchingAreas[i],0,matchingAreas[i],max);
            cout << "\t\t\t" << matchingAreas[i] << " ";
            if((i+1)%2==0)cout << "\n" << endl;
            larea->SetLineColor(kRed);
            larea->SetLineWidth(2);
            larea->Draw();
        }
    }
    
    
    
    string sampleName = to_string(static_cast<Int_t>(my_events[aux_events]))+"_"+to_string(channel);
    
    fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");

    
}



// ____________________________________________________________________________________________________ //
Bool_t snap(){
    if(eventNow == my_events[aux_events]){
        return true;
    }
    else{
        return false;
    }
    return false; // just to be sure oO
}

void getMySample(){
    TFile *fsample = new TFile("mysample.root","READ");
    if (fsample->IsZombie()) {
        matching = false;
        std::cout << "\n\n\n\n\n Error opening file \n\n\n\n\n" << std::endl;
        return;
    }
    TTree *tsample = (TTree*)fsample->Get("t1");
    
    Double_t values;
    tsample->SetBranchAddress("values",&values);
    Int_t nsample = tsample->GetEntries();
    for(Int_t i = 0; i < nsample; i++){
        tsample->GetEntry(i);
        mySample.push_back(values);
    }
    fsample->Close();
    return;
    
}

Bool_t checkAreas(Double_t totalmany){ //return true if the points did not touched the area...
    Int_t nAreas = matchingAreas.size(); //always two points per area of course
    Int_t n = selected_time.size();
    for(Int_t i = n; i>n-totalmany; i--){
        for(Int_t j = 0; j<nAreas; j++){
            if(selected_time[i-1]>=matchingAreas[j] && selected_time[i-1]<=matchingAreas[j+1]){
                return true;
            }
            j++;
        }
    }
    return false;
}

void makeSimpleHistogram(string filename){


    sphe_charge = sphe_charge_ch0; // wave0
    sphe_charge2 = sphe_charge2_ch0; // wave0
    delta = sphe_charge2 - sphe_charge;
    sphe_std = sphe_std_ch0;
  
  if(creation){
    twvf->Branch("charge",&wvfcharge,"charge/D");
    twvf->Branch("wvf",&wvf,Form("wvf[%i]/D",memorydepth));
    twvf->Branch("valid",&valid,"valid/O");
  }
  cout << "reading: " << filename << endl;
  string rootfile = filename + ".root";
  
  TFile *f1 = new TFile(rootfile.c_str(),"READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  ADC_DATA ch;
  TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
  bch->SetAddress(&ch);
  Int_t nentries = t1->GetEntries();
  Double_t charge = 0;
  Bool_t noise = false;
  Int_t noise_hits = 0;
  Double_t max = -1e12;
  mean_waveforms.clear();
  mean_waveforms.resize(memorydepth,0);
  naverages = 0;
    
  ofstream ftmp;
  ftmp.open("valid_events.log",ios::out);
  for(Int_t i = 0; i<nentries; i++){
    bch->GetEvent(i);
    noise = false;
    max = -1e12;
    for(Int_t j = start/dtime; j<finish/dtime; j++){
      charge += ch.wvf[j];
      if(ch.wvf[j]>=max){
        max = ch.wvf[j];
      }
      if(ch.wvf[j]>too_big) noise=true;
      if(ch.wvf[j]<lowerThreshold){
        noise_hits++;
        if(noise_hits>=maxHits){
          noise = true;
          break;
        }
      }
    }
    for(Int_t j = 0; j<start/dtime-800/4; j++){
      if(ch.wvf[j]<lowerThreshold){
        noise_hits++;
        if(noise_hits>=maxHits){
          noise = true;
          break;
        }
      }
    }
    for(Int_t j = finish/dtime+800/4; j<memorydepth; j++){
      if(ch.wvf[j]<lowerThreshold){
        noise_hits++;
        if(noise_hits>=maxHits){
          noise = true;
          break;
        }
      }
    }

    
    
    if(noise==false && ch.selection==0){
      valid = false;
      for(Int_t j = 0; j<memorydepth; j++){
        wvf[j] = ch.wvf[j];
      }
      if(charge*dtime>=delta/deltaminus  && charge*dtime<=delta*deltaplus){
        valid = true;
        for(Int_t j = 0; j<memorydepth; j++){
          mean_waveforms[j]+=ch.wvf[j];
        }
        naverages++;
      }

      wvfcharge = charge*dtime;
      twvf->Fill();
      hcharge->Fill(charge*dtime);
      ftmp << i << "\n";
    }
    charge=0;
  }
  ftmp.close();

  if(get_wave_form){
    for(Int_t i = 0; i<memorydepth; i++){
      timeg.push_back(dtime*i);
      mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
    }
    cout << "A total of " << naverages << " waveforms where found "<< endl;
        
  }

  
  TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);
    if(get_wave_form){
//         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
    }
    fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

    f1->Close();
    
}


};




class MeanSignal{
  
public:
  
  Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
  Int_t nbits = 14; // DIGITIZER bits
  
  vector<Int_t> channels = {1,2};
  Double_t minval = 0;
  Double_t maxval = 1e12;
  Double_t avoid_saturation = 1e6;
  Double_t fprompt = 0.8; // well... fprompt
  string mustbe = "bigger"; // must be 'bigger' or 'lower' then fprompt set.
  
Bool_t checkFprompt(Double_t fprompt_data, Double_t fprompt_ref, string mustbe){
  if(mustbe == "bigger"){
    if(fprompt_data>=fprompt_ref) return true;
    else return false;
  }
  else{
//     cout << "Test ... " << fprompt_data << " " << fprompt_ref << endl;
    if(fprompt_data<=fprompt_ref){
      return true;
    }
      else return false;
  }
}
void mean_signal(string filename){
  
  string rootfile = filename + ".root";
  
  TFile *f1 = new TFile(rootfile.c_str(),"READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  vector<ADC_DATA> ch(channels.size());
  vector<TBranch*> bch(channels.size());
  for(Int_t k = 0; k<channels.size();k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
  }
  Int_t nentries = t1->GetEntries();
  vector<Double_t> norm(channels.size(),0);
  
  vector<vector<Double_t>> avg(channels.size(),std::vector<Double_t>(memorydepth,0));
  vector<vector<Double_t>> avgn(channels.size(),std::vector<Double_t>(memorydepth,0));
  vector<Double_t> time(memorydepth);
  for(Int_t j = 0; j<memorydepth; j++){
    time[j]+=j*dtime;
  }
  
  TFile *fout = new TFile("averaged_waveforms.root","RECREATE");
  vector<TH2D*> hpersistence(channels.size());
  for(Int_t k = 0; k<channels.size();k++) hpersistence[k] = new TH2D(Form("hpersistence_%i",channels[k]),Form("hpersistence_%d",channels[k]),5000,0,20000,500,-500,avoid_saturation);

  for(Int_t i = 0; i<nentries; i++){
    for(Int_t k = 0; k<channels.size();k++){
      bch[k]->GetEvent(i);
      
      if(ch[k].charge>minval && ch[k].charge<maxval && ch[k].peak<avoid_saturation){
        if(checkFprompt(ch[k].fprompt,fprompt,mustbe)){
          norm[k]+=1;
          for(Int_t j = 0; j<memorydepth; j++){
            hpersistence[k]->Fill(j*dtime,ch[k].wvf[j]);
            avg[k][j]+=ch[k].wvf[j];
          }
        }
      }
    }
  }
  
  vector<TGraph*> gavg(channels.size());
  vector<TGraph*> gavgn(channels.size());
  vector<Double_t> maxvalue(channels.size());
  for(Int_t k = 0; k<channels.size(); k++){
    maxvalue[k] = *std::max_element(begin(avg[k]),end(avg[k]))/norm[k];
    for(Int_t j = 0; j<memorydepth; j++){
      avg[k][j]=avg[k][j]/norm[k];
      avgn[k][j]=avg[k][j]/maxvalue[k];
    }
    cout << "Ch" << channels[k] << " total waveforms = " << norm[k] << endl;
    gavg[k] = new TGraph(memorydepth,&time[0],&avg[k][0]);
    gavgn[k] = new TGraph(memorydepth,&time[0],&avgn[k][0]);
    fout->WriteObject(gavg[k],Form("average_ch%i",channels[k]),"TObject::kOverwrite");
    fout->WriteObject(gavgn[k],Form("average_normalized_ch%i",channels[k]),"TObject::kOverwrite");
    fout->WriteObject(hpersistence[k],Form("hpersistence_%i",channels[k]),"TObject::kOverwrite");
  }
  
  
  
}
  
};
















class Resolution{
  
public:

vector<Int_t> channels = {1,2};
TFile *fout;
TTree *tout;
Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
Int_t nbits = 14; // DIGITIZER bits

Double_t intmin = 5600;
Double_t intmax = 20000;
Double_t offset = 400; // just so it is easier
Double_t intstep = 1000;
Int_t nints = (intmax-(intmin+offset))/intstep+4;
vector<Double_t> integrations;
vector<vector<Double_t>> resCharge;
vector<Double_t> lowlim = {450,530,630,700,730,820,850,850,850,850,850,850,850,850,850,850,850,850};
vector<Double_t> average = {500,550,600,800,850,850,900,950,950,950,950,950,950,950,950,950,950,950};
Double_t ampl = 7e4;
void trackResolution(string filename){
  integrations.resize(nints);
  resCharge.resize(channels.size(),std::vector<Double_t>(nints,0));
  
  string rootfile = filename + ".root";
  
  TFile *f1 = new TFile(rootfile.c_str(),"READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  vector<ADC_DATA> ch(channels.size());
  vector<TBranch *> bch(channels.size());
  for(Int_t k = 0; k<channels.size(); k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
  }
  Int_t nentries = t1->GetEntries();
  
  fout = new TFile("resolution.root","RECREATE");
  tout = new TTree("resolution","tout");
  

  integrations[0] = offset+100;
  integrations[1] = integrations[0]+100;
  integrations[2] = integrations[1]+200;
  integrations[3] = integrations[2]+200;
  integrations[4] = offset+intstep;
  for(Int_t i = 5; i<nints; i++) integrations[i] = integrations[i-1]+intstep;
  vector<TH1D*> hspecs(nints);
  vector<TF1*> fa(nints);


  for(Int_t i = 0; i<nints; i++){
    hspecs[i] = new TH1D(Form("spec_%.0f",integrations[i]),Form("t = %.0f ns",integrations[i]),400,-100,1900);
    fa[i] = new TF1(Form("fa_run%.0f",integrations[i]),"([0]/(2*[2]))*exp((x-[1])/[2]+[3]*[3]/(2*[2]*[2]))*TMath::Erfc(((x-[1])/[3]+[3]/[2])/TMath::Power(2.,0.5))",lowlim[i],1300);
  }
  
  
  
  ifstream fsphe;
  fsphe.open("sphe.txt",ios::in);
  Double_t temp;
  vector<Double_t> sphes(channels.size());
  for(Int_t k = 0; k<channels.size(); k++){
    fsphe >> temp >> temp >> temp;
      sphes[k] = temp;
  }
  fsphe.close();
  
  
  Double_t charge1,charge2;
  Double_t photoelec;
  for(Int_t i = 0; i<nentries; i++){
//   for(Int_t i = 0; i<1; i++){
    cout << "reading event " << i << "\r" << flush;
    photoelec = 0;
    for(Int_t k = 0; k<channels.size(); k++){
      bch[k]->GetEvent(i);
      integrate(ch[k].wvf,resCharge[k]);
    }
    if(ch[0].fprompt>0.5 && ch[1].fprompt>0.5 && ch[0].peak<2800 && ch[1].peak<2800){
      for(Int_t j = 0; j<nints; j++){
        for(Int_t k = 0; k<channels.size(); k++){
          photoelec += resCharge[k][j]/sphes[k];
        }
        hspecs[j]->Fill(photoelec);
        photoelec=0;
      }
    }


    
  }
  f1->Close();
  
  vector<Double_t> sigmu(nints);
  vector<Double_t> ersigmu(nints);
  vector<Double_t> mu(nints);
  vector<Double_t> ermu(nints);
  THStack *hs = new THStack();
  for(Int_t j = 0; j<nints; j++){
    fa[j]->SetParameters(ampl,average[j],150,30);
    fa[j]->SetParName(0,"A");
    fa[j]->SetParName(1,"#mu");
    fa[j]->SetParName(2,"#tau");
    fa[j]->SetParName(3,"#sigma");
    hspecs[j]->Fit(Form("fa_run%.0f",integrations[j]),"R0Q");
    hspecs[j]->Fit(Form("fa_run%.0f",integrations[j]),"R0Q");
    fa[j]->SetRange(0,1300);
    sigmu[j] = fa[j]->GetParameter(3)/fa[j]->GetParameter(1);
    ersigmu[j] = sigmu[j]*sqrt(pow(fa[j]->GetParError(3)/fa[j]->GetParameter(3),2)+pow(fa[j]->GetParError(1)/fa[j]->GetParameter(1),2));
    mu[j] = fa[j]->GetParameter(1);
    ermu[j] = fa[j]->GetParError(1);
    fout->WriteObject(hspecs[j],Form("spec_%.0f",integrations[j]),"TObject::kOverwrite");
    fout->WriteObject(fa[j],Form("spectrum_run%.0f",integrations[j]),"TObject::kOverwrite");
    hspecs[j]->SetLineWidth(2);
    hs->Add(hspecs[j],"hist");
    
  }
  
//   gStyle->SetPalette(kVisibleSpectrum);
  TGraphErrors *g = new TGraphErrors(nints,&integrations[0],&sigmu[0],0,&ersigmu[0]);
  g->GetXaxis()->SetTitle("Integration time (ns)");
  g->GetYaxis()->SetTitle("#sigma/#mu");
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.7);
  TCanvas *c1 = new TCanvas();
  g->Draw("AP");
  fout->WriteObject(g,"resolution","TObject::kOverwrite");
  TGraphErrors *g2 = new TGraphErrors(nints,&integrations[0],&mu[0],0,&ermu[0]);
  g2->GetXaxis()->SetTitle("Integration time (ns)");
  g2->GetYaxis()->SetTitle("#mu");
  g2->SetMarkerStyle(21);
  g2->SetMarkerSize(0.7);
  TCanvas *c2 = new TCanvas();
  g2->Draw("AP");
  fout->WriteObject(g2,"peaks","TObject::kOverwrite");
  TCanvas *c3 = new TCanvas();
  hs->Draw("nostack plc");
  c3->BuildLegend();
    
}

void integrate(Double_t wvf[], vector<Double_t> &charges){
  Double_t refCharge = 0;
  Int_t index_ints = 0;
  for(Int_t i = 0; i<memorydepth; i++){
    if(i>=intmin/dtime && i<(intmin+integrations[index_ints])/dtime){
      refCharge+=wvf[i]*dtime;
      if(i==memorydepth-1){
        charges[index_ints] = refCharge;
        index_ints++;
      }
      
    }
    else if(i==(intmin+integrations[index_ints])/dtime){
//       cout << "filling i = " << i*dtime << " "  << integrations[index_ints] << " with " << refCharge << endl;
      charges[index_ints] = refCharge;
      refCharge+=wvf[i]*dtime;
      index_ints++;
    }
    
  }
  
}
  
  
  
};









































class TimeDistribuction{
public:
  
  vector<Int_t> channels={1,2};
  Double_t dtime = 4;
  Int_t interactions = 15;
  Bool_t just_a_test = false;
  Int_t just_this = 200;
  Int_t nshow = 20;
  Double_t lowerThreshold = 10;
  Double_t gap = 48; // in ns
  vector<Double_t> gtime;
  
  
vector<double> delay_line(Double_t v[], Int_t delay_time){
  vector<double> res(memorydepth);
  for(int i=0; i<memorydepth; i++){
    res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
  }
  return res;
}

void drawGraphs(TFile *f, TGraph * g1, TGraph * g2, Int_t count,Int_t k){ 
  TCanvas *c1 = new TCanvas();
  g1->SetLineColor(kBlue);
  g2->SetLineColor(kRed);
  g1->Draw("ALP");
  g2->Draw("LP SAME");
  f->WriteObject(c1,Form("g%i_ch%i",count,channels[k]),"TObject::kOverwrite");
}
  
    
// ____________________________________________________________________________________________________ //
void time_distribuction(string filename){
  gROOT->SetBatch(kTRUE);
  Int_t nchannels = channels.size();
  vector<vector<Double_t>> peaks(channels.size(),std::vector<Double_t>(memorydepth,0));
  vector<vector<Double_t>> shifted(channels.size(),std::vector<Double_t>(memorydepth,0));
  vector<TGraph *> gsample(nshow*nchannels);
  vector<TGraph *> gshifted(nshow*nchannels);
  vector<TH1D*> htd(nchannels);
  vector<Int_t> posaux(nchannels);
  vector<vector<Double_t>> position(nchannels);
  
  string rootfile = filename + ".root";
  TFile *f1 = new TFile(rootfile.c_str(),"READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  
  TFile *fout = new TFile("time_distribuction.root","RECREATE");
  
  vector<ADC_DATA>  ch(nchannels);
  vector<TBranch *> bch(nchannels);
  for(Int_t k = 0; k<nchannels;k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
    htd[k] = new TH1D(Form("h_ch%i",channels[k]),Form("h_ch%i",channels[k]),memorydepth,0,dtime*memorydepth);
  }
  for(Int_t j = 0; j<memorydepth; j++){
    gtime.push_back(dtime*j);
  }
  
  
  Int_t nentries = t1->GetEntries();
  if(just_a_test) nentries = just_this;
  for(Int_t i = 0; i<nentries; i++){
//     cout << "scanning event: " << i  << endl;
    cout << "scanning event: " << i << "\r" << flush;
    for(Int_t k = 0; k<nchannels;k++){
      bch[k]->GetEvent(i);
    }
    
    
    for(Int_t k = 0; k<nchannels;k++){
      shifted[k]=delay_line(ch[k].wvf, 40);
      smoothWithMovingAvarage(shifted[k]);
    }
    for(Int_t k = 0; k<nchannels;k++){
      for(Int_t j = 1; j<memorydepth-1; j++){ 
        if(shifted[k][j]>=lowerThreshold ){
          if(shifted[k][j]>shifted[k][j-1]){
            position[k].push_back(searchMax(ch[k].wvf,shifted[k],j)*dtime);
            posaux[k]++;
          }
        }
      }
    }
    if(i<nshow){
      for(Int_t k = 0; k<nchannels;k++){
        gsample[nchannels*i+k] = new TGraph(memorydepth,&gtime[0],&ch[k].wvf[0]);
        gshifted[nchannels*i+k] = new TGraph(memorydepth,&gtime[0],&shifted[k][0]);
      }
    }
  }
  
  cout << "\n";
  Int_t aux = 0;
  
  for(Int_t k = 0; k<channels.size(); k++){
    for(Int_t i = 0; i<position[k].size();i++){
      htd[k]->Fill(position[k][i]);
    }
    fout->WriteObject(htd[k],Form("h_ch%i",channels[k]),"TObject::kOverwrite");
  }

  for(Int_t k = 0; k<channels.size(); k++){
    aux=0;
    for(Int_t i = 0; i<nshow; i++){
//       cout << "scanning event: " << channels.size()*i+k << "\r" << flush;
      drawGraphs(fout,gsample[channels.size()*i+k],gshifted[channels.size()*i+k],aux,k);
      aux++;
    }
  }
  gROOT->SetBatch(kFALSE);
  
  
}

Int_t searchMax(Double_t wvf[],vector<Double_t> ref,Int_t &j){
  Double_t res = -1e12;
  Int_t maxpos = 0;
  for(Int_t i = j; ;i++){
    if(i<memorydepth){
      if(wvf[i]>res){
        res = wvf[i];
        maxpos = i;
      }
    }
    if(ref[i]<=lowerThreshold){
      j = i;
      break;
    }
  }
  return maxpos;
}
  
  
  // ____________________________________________________________________________________________________ //
void smoothWithMovingAvarage(vector<Double_t> &shifted){
  vector<Double_t> peak_smooth;
  Int_t n = shifted.size();
  Double_t sum = 0;
  Int_t midpoint;
  Double_t width;
  if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
    midpoint = interactions/2+1;    //midpoint will be 5 here
    width = interactions+1;
  }
  else{
    midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
    width = interactions;
  }
  
  for(Int_t i = 0; i < n; i++){
    
    if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
      peak_smooth.push_back(shifted.at(i));
    }
    else if(i>0){
      for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
        sum = sum+shifted.at(j);
        //                 cout << sum << endl;
      }
      peak_smooth.push_back(sum/width);
    }
    else{
      peak_smooth.push_back(0);
    }
    sum=0;
  }
  
  for(Int_t i = 0; i<memorydepth; i++){
    shifted[i] = peak_smooth[i];
  }
  
  
}
  
  
  
  
  
  
  
  
};




















































































// backup of previous if

// if(statcharge[i]>=delta/2. && statcharge[i]<=delta*1.5 && notAGoodWaveform[i]==false){
//   naverages++;
//   for(Int_t j = 0; j<waveforms.size(); j++){
//     if(j<temp_waveforms[i].size()){
//       if(j<=50){
//         newbase += temp_waveforms[i][j];
//       }
//       else if(j==51){
//         
//         newbase/=51;
//         for(Int_t k=0; k<=51; k++){
//           waveforms[k] = temp_waveforms[i][k]-newbase;
//           wvf[k]=waveforms[k];
//           mean_waveforms[k]+=waveforms[k];
//         }
//       }
//       else{
//         waveforms[j] = temp_waveforms[i][j]-newbase;
//         wvf[j]=waveforms[j];
//         mean_waveforms[j]+=waveforms[j];
//       }
//     }
//     else{
//       waveforms[j] = 0;
//       wvf[j]=0;
//       mean_waveforms[j]+=waveforms[j];
//     }
//     
//   }
//   twvf->Fill();
//   
//   gwaveforms[i] = new TGraph(waveforms.size(),&timeg[0],&waveforms[0]);
//   gm->Add(gwaveforms[i],"LP");
// }
              
