#include <fstream>
#include <iostream>

//Used to select and show single waveforms---
//Shuaixiang Zhang; Jul 30, 2022---

void select_single_waveform(){
  int headbin; // to store headers
  int nbytes_headers = 4; // 4 bytes (32 bits) for each head 
  int nbytes_data = 2; // 2 bytes (16 bits) per sample
  int memorydepth = 0; // size of waveforms 
  uint32_t valbin = 0; // to read data
  vector<Double_t> raw; // waveform as vector

  int waveNum = 0;//To record num of waveforms; Added by szhang; 20220720---
  int select1 = 30;
  int select2 = 70;
  int select3 = 140;
  int select4 = 210;
  int select5 = 280;
  int select6 = 350;
  int select7 = 420;
  int select8 = 490;
  int select9 = 560;
  int select10 = 630;
  int select11 = 4619;
  int select12 = 4620;
  int select13 = 4621;
  int select14 = 4874;
  int select15 = 4875;


  Double_t nbits = 14; // ADC is a 14 bits, 2 Vpp
  Double_t samplingRate = 500.e6; // 250 MSamples/s for DT5725
  Double_t nADCs = pow(2,nbits); // number of digital channels
  Double_t inVolts = 2./nADCs; // Multiply by this number to have the amplitude in volts;
  Double_t dtime = (1/samplingRate)*1e9; // steps in nanoseconds

  
  Bool_t first_line = true; // so we can set the length of the vector
  TH2D *h1;
  TH2D *h2;
  TH2D *h3;
  TH2D *h4;
  TH2D *h5;
  TH2D *h6;
  TH2D *h7;
  TH2D *h8;
  TH2D *h9;
  TH2D *h10;
  TH2D *h11;
  TH2D *h12;
  TH2D *h13;
  TH2D *h14;
  TH2D *h15;
 
  
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
      h1 = new TH2D("h1","h1",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h2 = new TH2D("h2","h2",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h3 = new TH2D("h3","h3",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h4 = new TH2D("h4","h4",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h5 = new TH2D("h5","h5",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h6 = new TH2D("h6","h6",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h7 = new TH2D("h7","h7",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h8 = new TH2D("h8","h8",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h9 = new TH2D("h9","h9",memorydepth,0,memorydepth*dtime,nADCs,6500,6700);
      h10 = new TH2D("h10","h10",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
      h11 = new TH2D("h11","h11",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
      h12 = new TH2D("h12","h12",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
      h13 = new TH2D("h13","h13",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
      h14 = new TH2D("h14","h14",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
      h15 = new TH2D("h15","h15",memorydepth,0,memorydepth*dtime,nADCs,6400,6700);
     // h = new TH2D("h", "h", memorydepth, 0, memorydepth*dtime, nADCs, 6000, 10000);
    }

    waveNum += 1;
    printf("Num of waveforms: %d\n", waveNum);
    if(waveNum>select15)
	break;

    for(int j = 0; j < memorydepth; j++){
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
	if(waveNum == select1)
	    h1->Fill(j*dtime,raw[j]);
	if(waveNum == select2)
	    h2->Fill(j*dtime,raw[j]);
	if(waveNum == select3)
	    h3->Fill(j*dtime,raw[j]);
	if(waveNum == select4)
	    h4->Fill(j*dtime,raw[j]);
	if(waveNum == select5)
	    h5->Fill(j*dtime,raw[j]);
	if(waveNum == select6)
	    h6->Fill(j*dtime,raw[j]);
	if(waveNum == select7)
	    h7->Fill(j*dtime,raw[j]);
	if(waveNum == select8)
	    h8->Fill(j*dtime,raw[j]);
	if(waveNum == select9)
	    h9->Fill(j*dtime,raw[j]);
	if(waveNum == select10)
	    h10->Fill(j*dtime,raw[j]);
	if(waveNum == select11)
	    h11->Fill(j*dtime,raw[j]);
	if(waveNum == select12)
	    h12->Fill(j*dtime,raw[j]);
	if(waveNum == select13)
	    h13->Fill(j*dtime,raw[j]);
	if(waveNum == select14)
	    h14->Fill(j*dtime,raw[j]);
	if(waveNum == select15)
	    h15->Fill(j*dtime,raw[j]);
        //printf("%d %.0f \n",j,raw[j]);
    }



  }

  TCanvas *c1 = new TCanvas();
  h1->Draw("colz");
  h1->GetXaxis()->SetTitle("Time (ns)");
  h1->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c2 = new TCanvas();
  h2->Draw("colz");
  h2->GetXaxis()->SetTitle("Time (ns)");
  h2->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c3 = new TCanvas();
  h3->Draw("colz");
  h3->GetXaxis()->SetTitle("Time (ns)");
  h3->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c4 = new TCanvas();
  h4->Draw("colz");
  h4->GetXaxis()->SetTitle("Time (ns)");
  h4->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c5 = new TCanvas();
  h5->Draw("colz");
  h5->GetXaxis()->SetTitle("Time (ns)");
  h5->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c6 = new TCanvas();
  h6->Draw("colz");
  h6->GetXaxis()->SetTitle("Time (ns)");
  h6->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c7 = new TCanvas();
  h7->Draw("colz");
  h7->GetXaxis()->SetTitle("Time (ns)");
  h7->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c8 = new TCanvas();
  h8->Draw("colz");
  h8->GetXaxis()->SetTitle("Time (ns)");
  h8->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c9 = new TCanvas();
  h9->Draw("colz");
  h9->GetXaxis()->SetTitle("Time (ns)");
  h9->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c10 = new TCanvas();
  h10->Draw("colz");
  h10->GetXaxis()->SetTitle("Time (ns)");
  h10->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c11 = new TCanvas();
  h11->Draw("colz");
  h11->GetXaxis()->SetTitle("Time (ns)");
  h11->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c12 = new TCanvas();
  h12->Draw("colz");
  h12->GetXaxis()->SetTitle("Time (ns)");
  h12->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c13 = new TCanvas();
  h13->Draw("colz");
  h13->GetXaxis()->SetTitle("Time (ns)");
  h13->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c14 = new TCanvas();
  h14->Draw("colz");
  h14->GetXaxis()->SetTitle("Time (ns)");
  h14->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c15 = new TCanvas();
  h15->Draw("colz");
  h15->GetXaxis()->SetTitle("Time (ns)");
  h15->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

//  c1->SaveAs("wave1.png");

}
