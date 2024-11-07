
class CT_MyFunctionObject{
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




vector<Double_t> probs;
vector<Double_t> er_probs;
vector<Double_t> phs;
Int_t n_CT = 3;
Int_t manytimes = 0;
const Int_t nparameters_CT = 3;

Double_t factorial(Int_t k){
  Double_t fac = 1;
  for(int i = 1; i <=k; ++i)
  {
    fac *= i;
  }
  return fac;
}
void evaluateProbs(vector<Double_t> &y, Double_t par[]){
  Double_t sum = 0;
  Double_t B = 0;
  for(Int_t k = 0; k<n_CT; k++){
    for(Int_t i = 0; i<=k;i++){
      if(i==0 && k ==0){
        B = 1;
      }
      else if(i==0 && k>0){
        B=0;
      }
      else{
        B = factorial(k-1)/(factorial(i)*factorial(i-1)*factorial(k-i));
      }
//       cout << k << " " << i << " "  << B << endl;
      sum += B*pow(par[1]*(1-par[2]),i)*pow(par[2],k-i);
    }
    y[k] = par[0]*TMath::Exp(-par[1])*sum;
//     cout << y[k] << " ";
    sum = 0;
  }
//   cout << endl;
  
}

void fcn_CT(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   //calculate chisquare
  Double_t chisq = 0;
  Double_t delta =0;
// 
  vector<Double_t> y(n_CT);
  evaluateProbs(y,par);
  for(Int_t i=0; i<n_CT; i++) {
      delta  = (probs.at(i)-y.at(i))/er_probs.at(i);
      //             cout << i << " " << z[i] << " " << y[i] << endl;
      chisq += delta*delta;
  }
  f = chisq;
  manytimes++;
  
  cout << manytimes << ", " << f;
  
  for(Int_t i = 0; i < nparameters_CT; i++){
    cout << ", " << par[i];
  }
  cout << "\n";
}



























class CT_Calibration

{
    
    
    
public:
    
    
    Double_t dtime = 2; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits
    Int_t linhasEvento = 9000;
    
    
    
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
    
    Int_t rebin = 4;
    
    Bool_t fixZero = false;
    
    
    Double_t sphe_charge = 0; // wave0
    Double_t sphe_charge2 = 0; // wave0
    
    
    // This was add here to try fitting dark noise data with already found values
    Bool_t darknoise = false;
    
    Bool_t is_poisson_test = false; // if running tests of poisson statistics

    
    
    
    
    // _______________ Parameters for smoothing  _______________/
    Int_t smooth_factor = 35;
    




    Double_t startFirst = 5000;
    Double_t finishFirst = 15000;












// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //




















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
    
    TFile *f1 = new TFile(rootFile.c_str(),"READ");
    TH1D *hcharge = (TH1D*)f1->Get(histogram.c_str());
    
    hcharge->Rebin(rebin);


    // ____________________________ Start of sphe fit ____________________________ //
    hcharge->GetYaxis()->SetTitle("Normalized count");
    hcharge->GetYaxis()->SetTitleOffset(1.0);
    hcharge->GetXaxis()->SetTitle("Charge (ADC*nsec)");
    
    Int_t nbins = hcharge->GetNbinsX();
    Double_t hmax = hcharge->GetBinCenter(nbins);
    Double_t hmin = hcharge->GetBinCenter(1);
    Double_t steps = (hmax-hmin)/(nbins-1); // this is the same as to do "GetBinCenter(2) - GetBinCenter(1)"
    
    
    TCanvas *c1 = new TCanvas("c1","Carga");
    c1->SetLogy();
    gPad->SetGrid(1,1);
    gPad->SetTicks(1,1);
    gStyle->SetOptFit();
    
    Double_t scale = 1/(hcharge->Integral());
    hcharge->Scale(scale);
    
    ofstream out;
    out.open("sphe.txt", ios::app);

    // ____________________________ Start of sphe fit ____________________________ //
    
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
    
    CT_MyFunctionObject MyFunc;
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
        aux=aux+2;
    }

    
    lastOne->SetParName(4,"#mu");
    lastOne->SetParName(5,"#sigma");
    lastOne->SetParName(7,"#mu2");
    
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

    
//     hcharge->GetXaxis()->SetRangeUser(-5000,15000);
    hcharge->StatOverflows(kTRUE);
    lastOne->SetRange(xmin,xmax);
    
    
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
    out <<  lastOne->GetParameter(4) << " " << lastOne->GetParameter(7) << " " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
    
    // ____________________________ Finish of sphe fit ____________________________ //
    
    //_________________ Drawing lines for the cross-talk probability _________________ //
    
    sphe_charge = lastOne->GetParameter(4);
    sphe_charge2 = lastOne->GetParameter(7);
    
    Double_t valley = 1e5;
    Double_t valley_pos = 0;
    for(Int_t i = 0; i<sphe_charge/2; i++){
      if(valley>lastOne->Eval(i)){
        valley = lastOne->Eval(i);
        valley_pos = i;
      }
    }
    Double_t delta1 = (sphe_charge2 - sphe_charge)/2;
    Double_t delta2 = sphe_charge+(sphe_charge2 - sphe_charge)/2;
    
    Double_t ymax = hcharge->GetMaximum();

    TLine *l1 = new TLine(valley_pos,0,valley_pos,ymax);
    TLine *l2 = new TLine(delta2,0,delta2,ymax);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    
    l1->Draw("");
    //     l2->Draw("");
    

    
    
    
    
    cout << "\n\n\nMaking analysis of CT: " << endl;
    
    
    
    Int_t valley_bin = (valley_pos-hmin)/steps; 
    //         cout << "valley real = " << valley_pos << " with bin = " << hcharge->GetBinCenter(valley_bin+1) << endl;
    
    Double_t zeroAmp = lastOne->GetParameter(0);
    Double_t zeroMean = lastOne->GetParameter(1);
    Double_t zeroSigma = lastOne->GetParameter(2);
    
    TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
    zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
    
    Double_t falseZero = zeroGaus->Integral(xmin,xmax);
    Double_t zeroIntegral = hcharge->Integral(1,valley_bin+1,"width"); // +1 is correct.. +2 could work.. small change here
    //         cout << falseZero << " " << zeroIntegral << endl;
    
    
    
    Double_t oneAmp = lastOne->GetParameter(3);
    Double_t oneMean = lastOne->GetParameter(4);
    Double_t oneSigma = lastOne->GetParameter(5);
    
    
    TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
    oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
    
    Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
    //         Double_t oneIntegral = oneGaus->Integral(valley_pos,xmax);
    
    Double_t totalevents = hcharge->Integral(1,nbins+1,"width");
    Double_t probZero = zeroIntegral/totalevents;
    Double_t probOne = oneIntegral/totalevents;
    
    cout << "total events = " << totalevents << endl;
    cout << "zero events = " << zeroIntegral << endl;
    cout << "sphe events = " << oneIntegral << endl;
    cout << "test = " << fu[1]->Integral(xmin,xmax) << endl;
    cout << "\n" << endl;
    cout << "prob zero = " << probZero << endl;
    cout << "prob one = " << probOne << endl;
    
    Double_t L = -TMath::Log(probZero);
    Double_t p = 1+probOne/(probZero*TMath::Log(probZero));
    cout << "L = " << L << endl;
    cout << "p = " << p << endl;
    
    
    if(n_CT>n_peaks+2) n_CT = n_peaks+2;
    probs.resize(n_CT);
    er_probs.resize(n_CT);
    phs.resize(n_CT);
    
    probs[0] = probZero;
    er_probs[0] = probZero*pow(1./zeroIntegral+1./totalevents,0.5);
    probs[1] = probOne;
    er_probs[1] = probOne*pow(1./oneIntegral+1./totalevents,0.5);
    phs[0] = 0;
    phs[1] = 1;
    
    for(Int_t i = 2; i<n_CT; i++){
      probs[i] = fu[i]->Integral(xmin,xmax)/totalevents;
      er_probs[i] = probs[i]*pow(1./fu[i]->Integral(xmin,xmax)+1./totalevents,0.5);
      phs[i] = i;
    }
    
    
    TCanvas *c2 = new TCanvas();
    c2->cd();
    TMultiGraph *gm = new TMultiGraph();
    TGraphErrors *gprobs = new TGraphErrors(n_CT,&phs[0],&probs[0],0,&er_probs[0]);
    gprobs->GetXaxis()->SetRangeUser(-1,8);
    gprobs->SetNameTitle("Experimental probabilites","Experimental probabilites");
    gprobs->SetMarkerStyle(21);
    gprobs->SetMarkerSize(0.7);
    
    
    
    
    
    
  TMinuit *gMinuit = new TMinuit(nparameters_CT);  //initialize TMinuit with a maximum of 3 params
  gMinuit->SetFCN(fcn_CT); // setting fcn 
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  
    
  // Set starting values and step sizes for parameters
  static Double_t vstart[nparameters_CT] = {1,L,p};
  static Double_t step[nparameters_CT] =   {0.1,L/10.,p/10.};
  gMinuit->mnparm(0, "scale",vstart[0], step[0], 0.003,10,ierflg);
  gMinuit->mnparm(1, "L",vstart[1], step[1],0,0,ierflg);
  gMinuit->mnparm(2, "p",vstart[2], step[2], 0,0,ierflg);
  
  
  
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET PRI", arglist ,1,ierflg);
  
  // Now ready for minimization step
  arglist[0] = 100000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //     gMinuit->mnprin(3,amin);
  
  Double_t param[nparameters_CT];
  Double_t Erparam[nparameters_CT];
  
  for(Int_t i = 0; i < nparameters_CT; i++){
    gMinuit->GetParameter(i,param[i],Erparam[i]);
    cout << param[i] << " +/- " << Erparam[i] << endl;        
  }
  
  vector<Double_t> fitted(n_CT);
  evaluateProbs(fitted,param);
  
  TGraph *gfit = new TGraph(n_CT,&phs[0],&fitted[0]);
  gfit->SetNameTitle("Fitted probabilites","Fitted probabilites");
  gfit->SetMarkerStyle(23);
  gfit->SetMarkerColor(kRed);
  gm->Add(gprobs,"P");
  gm->Add(gfit,"P");
  gm->Draw("A P");
  gm->GetXaxis()->SetTitle("n-phe");
  gm->GetYaxis()->SetTitle("Probability");
  c2->BuildLegend();
  
  
  // DANTE method... 
  // 
  cout << "\n\n\nDANTE method: " << endl;
  Double_t lambda = -log(fu[0]->Integral(xmin,xmax)/totalevents);
  Double_t meanexp = hcharge->GetMean();
  Double_t gain1 = oneMean-zeroMean;
  Double_t gain2 = lastOne->GetParameter(7)-oneMean;
  cout << "lambda = " << lambda << endl;
  cout << "meanexp = " << meanexp << endl;
  cout << "Method 1 (amplitude) = " << meanexp/(lambda*gain1) << endl;
  cout << "Method 2 (difference) = " << meanexp/(lambda*gain2) << endl;
//   cout << lambda << endl;
  
}

    
};









// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //















// ____________________________________________________________________________________________________ //


















































class CT_SPHE{
       
    
    
    
public:
    
    
TFile *fout;

TTree *tout;

Bool_t creation = true; //to verify creation of ttree branches




Double_t value = 0;
Double_t desv = 0;
Double_t channel = 0;


// ____________________ Variables to calculate ____________________ //

TH1D *hbase = new TH1D("hbase","histogram for baseline",5*800,-400,400);
TH1D *hbase_smooth = new TH1D("hbase_smooth","histogram for baseline smoothed",5*800,-400,400);
TH1D *hcharge = new TH1D("hcharge","",120000,-200000,2*1300000);
TH1D *hzero = new TH1D("hzero","",120000,-200000,2*1300000);
TH1D *hnobase = new TH1D("hnobase","",120000,-200000,2*1300000);

TH1D *hstat = new TH1D("hstat","",120000,-200000,2*1300000);


Double_t dtime = 2.;
Int_t npoints = 9000;

Double_t mean = 0;
Double_t stddev = 0;
Double_t tolerance; // n sigmas
Double_t baseLimit;

Double_t timeLimit; // time after LED signal
Double_t timeLow ; // integration time before peak
Double_t timeHigh ; // integration time after peak

Double_t start = 0;
Double_t finish = npoints*dtime;

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
Bool_t filled = false;
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
Double_t wvf[9000];
Int_t naverages;
Double_t mean_before = 100;
Double_t mean_after = 400;
TMultiGraph *gm = new TMultiGraph();
TGraph *gwaveforms;

Double_t sphe_charge;
Double_t sphe_charge2;
Double_t delta;
Double_t sphe_std;

Double_t sphe_charge_ch0;
Double_t sphe_charge_ch1;

Double_t sphe_charge2_ch0;
Double_t sphe_charge2_ch1;

Double_t sphe_std_ch0;
Double_t sphe_std_ch1;

string charge_status = "";

Bool_t darkNoise = false;

Double_t new_base = 0;
Double_t itegration_time = 0;
Double_t maximum_pos;


Double_t chargeOut = 50000;
// ____________________________________________________________________________________________________ //

void giveMeSphe_darkCount(string name){
    gROOT->SetBatch(kTRUE);
    fout = new TFile("sphe_histograms_darkCount.root","RECREATE");
    tout = new TTree("t1","baseline info");
    darkNoise = true;
    
    makeHistogram(name);
    channel=channel + 1;
    
    fout->WriteObject(tout,"t1","TObject::kOverwrite");
    
}


// ____________________________________________________________________________________________________ //
void giveMeSphe(string name){
    gROOT->SetBatch(kTRUE);
    
    fout = new TFile("sphe_histograms_CT.root","RECREATE");
    tout = new TTree("t1","baseline info");
    darkNoise = false;
    
    makeHistogram(name);
    channel=channel + 1;
    
    fout->WriteObject(tout,"t1","TObject::kOverwrite");
    
}
// ____________________________________________________________________________________________________ //
Bool_t checkForValidBaseline(){
    Double_t mydev = 0;
    value = charge;
    mydev = hfilter->GetStdDev();
    desv = mydev;
//     t1->Fill();
    if(desv>devFilter){
            return false;
    }
    else{
        return true;
    }
        
}
// ____________________________________________________________________________________________________ //
void integrateSignal(){
    charge = 0;
    Int_t npeaks = 1;
    hfilter->Reset();
    value = 0;
    desv = 0;

    strikes = 0;
    Int_t discardedPeaks = 0;

    Int_t aux_sample = 0;
    
    
    if(channel==0){
        sphe_charge = sphe_charge_ch0; // wave0
        sphe_charge2 = sphe_charge2_ch0; // wave0
        delta = sphe_charge2 - sphe_charge;
        sphe_std = sphe_std_ch0;
    }
    if(channel==1){
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
    vector<Double_t> waveforms(npoints);
    vector<TGraph*> gwaveforms(npeaks);
    vector<Bool_t> notAGoodWaveform(npeaks);
    //for matching
    Bool_t inBadArea[npeaks];

    Int_t auxstat = 0;
    Bool_t notGood = false;
    Double_t lastTime = 0 ;
    Double_t peakStd = 0;
    
    Double_t thismanypointsLow = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t thismanypointsHigh = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t thismanypointsBase = 0; //this is to check the amount of points, noise may cause low points that is terrible
    Double_t totalmany = 0; //this helps speacially for matching...
    
    for(Int_t i = 0; i<npeaks; i++){
      
      // if is the first peak
      // or the pois is beyond the previous peak position + low+high
      statcharge[i] = 0;
      statpeak[i] = 0;
      trackStrikes[i] = 0;
      discard_this[i] = false;
      inBadArea[i] = false;
      notGood = false;
      thismanypointsLow = 0;
      thismanypointsHigh = 0;
      thismanypointsBase = 0;
      totalmany = 0;
      
      notAGoodWaveform[i] = false;
      
      
      for(Int_t j = start/dtime; j<=finish/dtime; j++){
        charge += temp_peak.at(j);
        statcharge[i]+= temp_peak.at(j);
        if(temp_peak.at(j)>=statpeak[i]){
          statpeak[i] = temp_peak.at(j);
        }
        selected_peaks.push_back(temp_peak.at(j));
        selected_time.push_back(timeg.at(j));
        thismanypointsLow++;
        totalmany++;
        if(peak_smooth.at(j)<lowerThreshold){ // if there is any point below -X ADCs I discard
          trackStrikes[i]+=1;
        }
      }
      
      
      
      
      //             if(trackStrikes[i]>=maxHits || statcharge[i] == 0){
      //                 discard = true;
      if(charge*dtime>chargeOut && discard){
                      discard_this[i] = discard;
      }
      else{
                      discard_this[i] = false;
      }
      //                 discardedPeaks++;
      //             }
      //             else if(matching){ // if strikes are good, then check the area...
      //                 inBadArea[i] = checkAreas(totalmany);
      //                 
      //             }
      
      
      
      if(snap()){
        // //                 if(discard_this[i]==true){//discarting charge because pulse was a big fat noise
        // // //                     selected_peaks.erase(selected_peaks.end()-(thismanypointsBase),selected_peaks.end());
        // // //                     selected_time.erase(selected_time.end()-(thismanypointsBase),selected_time.end());
        // //                     cout << " Discarded -> ";
        // // //                     charge_status += " Discarded - > " + to_string(statcharge[i]*dtime);
        // //                 }
        // //                 cout << "charge = " << statcharge[i]*dtime << " at " << peakPosition.at(i) << " with " << thismanypointsLow << " + " << thismanypointsHigh << " or this " << thismanypointsBase << endl;
      }
      
    }

        
        
        
    
    for(Int_t i = 0; i<npeaks; i++){
      
        statcharge[i] = dtime*statcharge[i];
        strikes = trackStrikes[i];
        charge = dtime*charge;
        treeCharge = statcharge[i];
        treePeak = maximum_pos;
        ptsHigh = trackHigh[i];
        ptsLow = trackLow[i];
        filled = discard_this[i];
        tout->Fill();
                  
        if(discard_this[i]==false){
            
            if(matching){
                if(inBadArea[i]){
                    charge_status += " AREA - > " + to_string(statcharge[i]*dtime) + " ";
                    continue; //go to next peak... this one was no good...
                }
            }
            
            
            hcharge->Fill(statcharge[i]);
            hnobase->Fill(statcharge[i]);
            hstat->Fill(statcharge[i]);
            
            charge_status += " charge = " + to_string(statcharge[i]*dtime);
//             if(statcharge[i]>10000 && statcharge[i]<11000 && snap()){
//                 cout << ".............................. charge = " << statcharge[i] << " at " << peakPosition.at(i) << " event = " << eventNow << endl;
//             }
            
            if(snap()){
              gwaveforms[i] = new TGraph(temp_peak.size(),&timeg[0],&temp_peak[0]);
              gm->Add(gwaveforms[i],"LP");
                
            }
            
            
        }
        else{
            charge_status += " Discarded - > " + to_string(statcharge[i]*dtime) + " ";
            charge_status += " strikes - > " + to_string(trackStrikes[i]) + " ";
        }
    }
    
    

    

 
}

// ____________________________________________________________________________________________________ //
void searchForPeaks(){
    // For each peak found, we i am looking for the maximum above tolerance
    Int_t n = peak_smooth.size();
    Int_t npeaks = 0;
    Bool_t found_a_peak = false;
    
    Double_t time_zero = 0;
    
    Bool_t wait_a_low = false;
    
    Int_t m = mySample.size();
    Double_t sum = 0;
    Double_t sum2 = 0;
    Double_t sum3 = 0;
    
    // matching parameters
    Int_t tempTime1 = 0;
    Int_t tempTime2 = 0;
    matchingAreas.clear();
    higherValue = 0;
    tempPosition = 0;
    Double_t matchtime = 0;
    
    Int_t temp_i_peaks = 0;
    Bool_t notIgnoringPeaks = true;
    Int_t factorPeaks = 1;
    
    Bool_t notIgnoringMatching = true;
    Int_t temp_i_matching = 0;

    
    itegration_time = (finish-start)/dtime;
    Double_t maximum = -1e7;
    for(Int_t i = 0; i<n; i++){
        temp_peak[i] = temp_peak[i]-new_base;
        if(i>=start/dtime && i<=finish/dtime){
          if(maximum<=temp_peak[i]){
            maximum = temp_peak[i];
            maximum_pos = i*dtime;
          }
          
        }
    }
    discard = false;
    if(maximum_pos == 0){
      cout << maximum << " " << maximum_pos << endl;
    }
    if(maximum_pos<=3680){
      discard = true;
//       cout << " discarding peak at " << endl;
    }
    
    
    
}

// ____________________________________________________________________________________________________ //
void lookWaveform(){

    
    smoothWithMovingAvarage(); // first calculate avarage

    searchForPeaks(); //search for peaks with the moving avarage

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
        tout->Branch("filled",&filled,"filled/O");
        creation = false;
        
    }
    
    cout << "reading: " << filename << endl;
    string rootfile = filename + "_copy.root";
    
    TFile *f1 = new TFile(rootfile.c_str(),"READ");
    TTree *t1 = (TTree*)f1->Get("t1");
    Int_t nentries = t1->GetEntries();
    Int_t nevents = nentries/npoints;
    
    Double_t adc,peak,time, event;
    
    t1->SetBranchAddress("adc",&adc);
    t1->SetBranchAddress("peak",&peak);
    t1->SetBranchAddress("time",&time);
    t1->SetBranchAddress("event",&event);
    
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
    mean_waveforms.resize(npoints,0);
    naverages = 0;
    
    for(Int_t i = 0; i<npoints; i++){
        timeg.push_back(dtime*i);
    }
    
    if(just_a_test){nentries = just_this*npoints;}
//     aux = 0;
    for(Int_t i = aux*npoints; i<nentries; i++){
        t1->GetEntry(i);
        
        
        if(aux==event && i!=(nentries-1)){
        }
        else{ // new event
            if( i == (nentries-1) ){
                if(time<=timeLimit || time>=8000){
                    hbase->Fill(peak);
                }
                else if(time<forBaseline && abs(peak)<(5*baseLimit/15) && time>500  && badBaseline){
                    baseCharge += peak;
                }
                temp_peak.push_back(peak);
            }
            aux = event;
            if(static_cast<Int_t>(eventNow)%200==0){
                cout << eventNow << "\r" << flush;
            }
            eventNow =  eventNow + 1;
            if(badBaseline)hcharge->Fill(baseCharge*dtime);
            // here is were all the calculation is done !!!
            lookWaveform();
            
            if(snap()){
                aux_events = aux_events+1;
            }
            
            if(event!=nevents){
                hbase->Reset();
                hbase_smooth->Reset();
            }
            charge_status = "";
            temp_peak.clear();
            peak_smooth.clear();
            peakPosition.clear();
            peakMax.clear();
            selected_peaks.clear();
            selected_time.clear();
            baseCharge = 0;
            matched.clear();
    
            j=0;
        }
        
        
        if(i != (nentries-1)){
            
            if(time<=timeLimit || time>=8000){
                hbase->Fill(peak);
            }
            else if(time<forBaseline && abs(peak)<(10*baseLimit/15) && time>500 && badBaseline){
                baseCharge += peak;
            }
            temp_peak.push_back(peak);
        }
        
        j++;

    }
    
    
    cout << "____________________________ " << counter/nsample << " \n \n \n" << endl;
    
    if(get_wave_form){
        for(Int_t i = 0; i<npoints; i++){
//             cout << i << " " << mean_waveforms[i] << " ";
             
            mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
//             cout << mean_waveforms[i] << endl;
            
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;
        
    }
    
    TCanvas *cwvf = new TCanvas();
    cwvf->cd();
    gm->Draw("A");
    
    TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);


    fout->WriteObject(cwvf,Form("all valid waveforms of ch%.0f",channel),"TObject::kOverwrite");
    fout->WriteObject(hcharge,filename.c_str(),"TObject::kOverwrite");

    f1->Close();

    mark++;
}

// ____________________________________________________________________________________________________ //
void smoothWithMovingAvarage(){
    
    Int_t n = temp_peak.size();
    if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
    }
    else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
    }
    
    new_base = hbase->GetMean();

    for(Int_t i = 0; i < n; i++){
        
        if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
            peak_smooth.push_back(temp_peak.at(i)-new_base);
        }
        else if(i>cut/dtime){
            for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
                sum = sum+temp_peak.at(j)-new_base;
//                 cout << sum << endl;
            }
            peak_smooth.push_back(sum/width);
            
            if(timeg.at(i)<=timeLimit && abs(peak_smooth.at(i))<baseLimit){
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
    
    if(snap()){ // if the events match, create the graphs
 
        g_smooth = new TGraph(n,&timeg[0],&peak_smooth[0]);
        g_normal = new TGraph(n,&timeg[0],&temp_peak[0]);
//         cout << mean << " " << stddev << endl;
    }
    
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
    
    TLine *lmean = new TLine(timeLimit,mean,npoints*dtime,mean);
    TLine *ldev = new TLine(timeLimit,mean+tolerance*stddev,npoints*dtime,mean+tolerance*stddev);
    
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
    
    
    
    string sampleName = to_string(static_cast<Int_t>(my_events[aux_events]))+"_"+to_string(mark);
    
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


};
