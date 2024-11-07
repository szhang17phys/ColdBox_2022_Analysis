#include <fstream>
#include <iostream>

//Used to get the difference of different signal regions---
//Shuaixiang Zhang; Aug 2nd, 2022---

void subtract_sig(){
  int headbin; // to store headers
  int nbytes_headers = 4; // 4 bytes (32 bits) for each head 
  int nbytes_data = 2; // 2 bytes (16 bits) per sample
  int memorydepth = 0; // size of waveforms 
  uint32_t valbin = 0; // to read data
  vector<Double_t> raw; // waveform as vector

  int waveNum = 0;//To record num of waveforms; Added by szhang; 20220720---

  Double_t nbits = 14; // ADC is a 14 bits, 2 Vpp
  Double_t samplingRate = 500.e6; // 250 MSamples/s for DT5725
  Double_t nADCs = pow(2,nbits); // number of digital channels
  Double_t inVolts = 2./nADCs; // Multiply by this number to have the amplitude in volts;
  Double_t dtime = (1/samplingRate)*1e9; // steps in nanoseconds
  
  Bool_t first_line = true; //so we can set the length of the vector
  TH2D *h; 

  int binNum = 1050;

  TH1F *result1 = new TH1F("result1", "result1", binNum, -100000, 600000);
  TH1F *result2 = new TH1F("result2", "result2", binNum, -100000, 600000);
  TH1F *difference = new TH1F("difference", "Signal Integrals",binNum, -100000,600000);
//  TH1F *dis = new TH1F("dis", "Distribution of max ADC",1005, -100000,2000000);
 
  int sumBgd1 = 0;
  int sumBgd2 = 0;  
  int region1 = 0;
  int region2 = 0;
  int binContent1 = 0;
  int binContent2 = 0;
  int diff = 0;
  int counts = 0;

  ifstream fin;
  fin.open("/afs/cern.ch/user/s/shuaixia/private/PDS_Data/cathodexarapuca/wave5.dat", ios::in | ios::binary);

  if(fin.good() && fin.is_open()){ // Ok
    cout << "Reading file" << endl;
  }
  else{ // emergency shutdown
    cout << "File did not open!!" << endl;
    return;      
  }
  while(!fin.fail()){ 
    for(Int_t ln=0;ln<6;ln++){ // 4 bytes (32 bits) for each head 
      fin.read((char *) &headbin, nbytes_headers);
      // header0 will be EventSize, so: you can do
      if(ln==0){
        memorydepth = headbin;
        // the result is in bytes, so 2*NSamples+4*headers
        memorydepth = (memorydepth-4*6)/2;
      }
    }
    if(first_line){
      printf("Waveform size: %d \n",memorydepth);
      raw.resize(memorydepth);
      first_line=false;
      h = new TH2D("h","h",memorydepth,0,memorydepth*dtime,nADCs,0,nADCs);
    }
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
//        h->Fill(j*dtime,raw[j]);
        //printf("%d %.0f \n",j,raw[j]);
    }
    
//To integrate the signal and background---
    for(int i=1750; i<2000; ++i){
	sumBgd1 += raw[i];
    }
    for(int i=3250; i<3500; ++i ){
	sumBgd2 += raw[i];
    }

    for(int i=2500; i<2750; ++i){
	region1 +=raw[i];
    }
    region1 = region1 - (sumBgd1+sumBgd2)/2.0;  


    for(int i=3250; i<4000; ++i){
	region2 +=raw[i];
    }
    region2 = region2 - (sumBgd1+sumBgd2)/2.0;  

    result1->Fill(region1);
    result2->Fill(region2);

    sumBgd1 = 0;
    sumBgd2 = 0;
    region1 = 0;
    region2 = 0;
//-----------------------------------------

    waveNum += 1;
    printf("Num of waveforms: %d\n", waveNum);
    
//    if(waveNum>0)
//	break;

  }

//To get the difference of above two histos------
  for(int i=0; i< binNum; ++i){
      binContent1 = result1->GetBinContent(i);
      binContent2 = result2->GetBinContent(i);
      diff = binContent1 - binContent2;
//      if(diff > 0){
          counts += diff;
          difference->SetBinContent(i, diff);
//      }
  }


  TCanvas *c1 = new TCanvas();
  gPad->SetLogy();
  result1->Draw("colz");
  result1->GetXaxis()->SetTitle("[ADC*ns]");
  result1->GetYaxis()->SetTitle("Counts");
  result1->SetLineColor(kRed);
  result2->SetLineColor(kBlue);
//  result2->Draw();
  result1->Draw();

  c1->SaveAs("sigIntegral_1.png");


  TCanvas *c2 = new TCanvas();
  gPad->SetLogy();
  gStyle->SetOptStat("nemr");
  difference->Draw("colz");
  difference->GetXaxis()->SetTitle("(Subtraction)[ADC*ns]");
  difference->GetYaxis()->SetTitle("Counts");
  c2->SaveAs("sigIntegral_1.png");

  cout<<"------------------------"<<endl;
  cout<<"Total Counts: "<<counts<<endl;
  cout<<"------------------------"<<endl;


}
