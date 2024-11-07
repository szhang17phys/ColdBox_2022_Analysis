#include "MYCODES.h"


// class MeanSignal{
//   
//     
// public:
// 
// 
// Int_t linhasEvento = 9000;
// Int_t dtime = 2;
// 
// Int_t Max_LowLim = 3000;
// Int_t Max_HighLim = 4400;
// Double_t cutoff = 200;
// Double_t fcut = 0;
// 
// Bool_t all = false;
// Int_t all_channel = 0;
// 
// TFile *fout = new TFile("mean_signal.root","RECREATE");
// 
// 
// 
// void mean_signal(string filename){
//     string rootfile = filename + ".root";
//     string rootfile_copy = filename + "_copy.root";
//     
//     TFile *fcopy = new TFile(rootfile_copy.c_str(),"READ");
//     TFile *f = new TFile(rootfile.c_str(),"READ");
//     TTree *t1 = (TTree*)fcopy->Get("t1");
//     TTree *t2 = (TTree*)f->Get("t1");
//     
//     Int_t nentries = t1->GetEntries();
//     
//     Double_t peak;
//     Double_t fprompt;
//     
//     if(all){
//         if(all_channel==0){
//             t1->SetBranchAddress("peak0",&peak);
//             t2->SetBranchAddress("Fprompt0",&fprompt);
//         }
//         else if(all_channel==1){
//             t1->SetBranchAddress("peak1",&peak);
//             t2->SetBranchAddress("Fprompt1",&fprompt);
// 
//         }
//     }
//     else{        
//         t1->SetBranchAddress("peak",&peak);
//         t2->SetBranchAddress("Fprompt",&fprompt);
//     }
//     
//     Double_t y_mean[linhasEvento];
//     Double_t y[linhasEvento];
//     Double_t x_mean[linhasEvento];
//     
//     std::fill(y_mean,y_mean+linhasEvento,0); // set array with 0 in every entry
//     std::fill(x_mean,x_mean+linhasEvento,0);
//     
//     Int_t aux_mean = 0;
//     Bool_t checkAmplitude = true;
//     Bool_t getmean = true;
//     Double_t maxamplitude = 0;
//     Int_t contador = 0;
//     
//     Int_t aux_t2 = 0;
//     for(Int_t i = 0; i<nentries; i++){
//         
//         if(checkAmplitude){
//             checkAmplitude = false;
//             maxamplitude = 0;
//             for(Int_t j = i + Max_LowLim/dtime; j<=i+Max_HighLim/dtime; j++){
//                 t1->GetEntry(j);
//                 if(peak>=maxamplitude){
//                     maxamplitude = peak;
//                 }
//             }
//             if(maxamplitude>=cutoff && maxamplitude<14800){
//                 getmean = true;
//             }
//             else{
//                 getmean = false;
//             }
//             
//             
//             // here for the Fprompt analysis
//             t2->GetEntry(aux_t2);
//             if(fprompt<fcut){
//                 getmean = false;
//             }
//             aux_t2++;
//         }
//         
//         t1->GetEntry(i);
//         
//         if(getmean){
//             contador++;
//             y_mean[aux_mean] += peak;
//             x_mean[aux_mean] = dtime*aux_mean;
//         }
//         aux_mean++;
//         if(aux_mean==(linhasEvento)){
//             checkAmplitude = true;
//             aux_mean=0;
//         }
//     }
//     cout << contador/linhasEvento << endl;
//     contador = contador/linhasEvento;
//     Double_t maxvalue = *std::max_element(y_mean,y_mean+linhasEvento);
//     for(Int_t i = 0; i<linhasEvento; i++){
//         y[i] = y_mean[i]/(1.*contador);
//         y_mean[i] = y_mean[i]/maxvalue;
//     }
//     
//     f->Close();
//     
//     
//     TCanvas *c1 = new TCanvas();
//     c1->cd();
//     TGraph *g = new TGraph(linhasEvento,x_mean,y_mean);
//     g->Draw("ALP");
//     
//     TCanvas *c2 = new TCanvas();
//     c2->cd();
//     TGraph *g_normal = new TGraph(linhasEvento,x_mean,y);
//     g_normal->Draw("ALP");
//     
//     if(all){
//         filename = filename + "_" + to_string(all_channel);
//     }
//     
//     fout->WriteObject(g,filename.c_str(),"TObject::kOverwrite");
//     
//     filename = filename + "_normal";
//     
//     fout->WriteObject(g_normal,filename.c_str(),"TObject::kOverwrite");
//     
// }
// 
// };



















class TakeCharge{
  
    
public:


Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
Int_t nbits = 14; // DIGITIZER bits

vector<Int_t> channels = {1,2};
Double_t integration_start = 500;
Double_t integration_end = 1060;
Double_t hmin = 0;
Double_t hmax = 1000;
Double_t nbins = 200;



void takeCharge(string filename){
  
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
  vector<Double_t> charge(channels.size(),0);
  
  TFile *fout = new TFile("charges.root","RECREATE");
  vector<TH1D*> hcharge(channels.size());
  for(Int_t k = 0; k<channels.size();k++) hcharge[k] = new TH1D(Form("hcharge_%i",channels[k]),Form("hcharge_%d",channels[k]),nbins,hmin,hmax);
  
  for(Int_t i = 0; i<nentries; i++){
    for(Int_t k = 0; k<channels.size();k++){
      bch[k]->GetEvent(i);
      for(Int_t j = integration_start/dtime; j<integration_end/dtime; j++){
        charge[k] += dtime*ch[k].wvf[j];
      }
      hcharge[k]->Fill(charge[k]);
      charge[k]=0;
    }
  }
  
  
  
  
  f1->Close();
  
  for(Int_t k = 0; k<channels.size();k++) fout->WriteObject(hcharge[k],Form("hcharge_%i",channels[k]),"TObject::kOverwrite");
  fout->Close();
  
}

};












