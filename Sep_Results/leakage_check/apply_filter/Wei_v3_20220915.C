#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Total variation noise filter: https://hal.archives-ouvertes.fr/hal-00675043v2/document
void TV1D_denoise(vector<Double_t>& input, vector<Double_t>& output, const int width, const Double_t lambda){
  if (width>0) {                /*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
    Double_t umin=lambda, umax=-lambda;    /*u is the dual variable*/
    Double_t vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
    int kplus=0, kminus=0;     /*last positions where umax=-lambda, umin=lambda, respectively*/
    const Double_t twolambda=2.0*lambda;    /*auxiliary variable*/
    const Double_t minlambda=-lambda;        /*auxiliary variable*/
    for (;;) {                /*simple loop, the exit test is inside*/
      while (k==width-1) {    /*we use the right boundary condition*/
        if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
          do output[k0++]=vmin; while (k0<=kminus);
          umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
        } else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
          do output[k0++]=vmax; while (k0<=kplus);
          umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
        } else {
          vmin+=umin/(k-k0+1);
          do output[k0++]=vmin; while(k0<=k);
          return;
        }
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
        do output[k0++]=vmin; while (k0<=kminus);
        vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
        umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
        do output[k0++]=vmax; while (k0<=kplus);
        vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
        umin=lambda; umax=minlambda;
      } else {     /*no jump necessary, we continue*/
        k++;
        if (umin>=lambda) {        /*update of vmin*/
          vmin+=(umin-lambda)/((kminus=k)-k0+1);
          umin=lambda;
        }
        if (umax<=minlambda) {    /*update of vmax*/
          vmax+=(umax+lambda)/((kplus=k)-k0+1);
          umax=minlambda;
        }
      }
    }
  }
}

void Wei_v3_20220915(){
  //=============
  // User config
  //=============
  Bool_t first_line = true; // so we can set the length of the vector
  Bool_t filternoise = true; // apply noise filter if true (most likely you need)
  int example_waveform = 19470; // the example waveform number you want to plot: can't exceed maximum waveforms

  //==========================
  // Digitizer config
  //     usually no change
  //==========================
  int headbin; // to store headers
  int nbytes_headers = 4; // 4 bytes (32 bits) for each head
  int nbytes_data = 2; // 2 bytes (16 bits) per sample
  int memorydepth = 0; // size of waveforms
  uint32_t valbin = 0; // to read data

  int waveNum = 0;//To record num of waveforms; Added by szhang; 20220720---

  Int_t nbits = 14; // ADC is a 14 bits, 2 Vpp
  Double_t samplingRate = 250.e6; // 250 MSamples/s for DT5725
  Int_t nADCs = pow(2, nbits); // number of digital channels
  Double_t inVolts = 2./nADCs; // Multiply by this number to have the amplitude in volts;
  Double_t dtime = (1/samplingRate)*1e9; // steps in nanoseconds

  //===================
  // General Variables
  //===================
  vector<Double_t> raw; // waveform as vector
  vector<Double_t> denoised; // denoised waveform as vector
  Double_t baseline_ADC;
  Int_t baseline_ADC_bin;

  double level[8] = {0};
  double delta_ADC[8] = {0};
  double Count0[8] = {0};
  double Count1[8] = {0};
  double Count2[8] = {0};
  double Count3[8] = {0};
  double Count123[8] = {0};

  double Gain1[8] = {0};
  double Gain2[8] = {0};
  double Gain3[8] = {0};
  double Gain123[8] = {0};

  // delta ADC above baseline as threshold
  delta_ADC[0] = 50;
  level[0] = 3385;//for wave5, baseline=3335---
//  level[0] = 2590;//for wave6, baseline=2540---
  for(int i=0; i<8; ++i){
      level[i] = level[0]+25*i;
      delta_ADC[i] = delta_ADC[0] + 25*i;
  }

  // Single waveform
  TH1D *h1;
  TH1D *h2; // denoised

  //=================
  // Read No Laser
  //=================
  ifstream fin;
  fin.open("/afs/cern.ch/work/s/shiw/public/ColdBoxVD/Coldbox_Sep2022_leakage_check/20220915_v3_laser1_705_laser2_737_laser3_740/xarapuca_V3_cosmic_trigger50_no_laser/wave5.dat", ios::in | ios::binary);

  if(fin.good() && fin.is_open()) cout << "Reading file" << endl;
  else{
    cout << "File did not open!!" << endl;
    exit(1);
  }

  // Loop over waveforms
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
      denoised.resize(memorydepth);
      first_line=false;
      h1 = new TH1D("h1", "h1", memorydepth, 0, memorydepth);
      h2 = new TH1D("h2", "h2", memorydepth, 0, memorydepth);
    }

    // Store each waveform
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        denoised[j] = raw[j];
    }

    // Apply noise filter
    if (filternoise) TV1D_denoise(raw, denoised, memorydepth, 50);

    // Evaluate single waveform baseline
    baseline_ADC = 0.;
    baseline_ADC_bin = -1;
    TH1F *ADC_hist = new TH1F("ADC_hist", "ADC_hist", nADCs/4, 0, nADCs); // this bin can't be too fine, 4 adc binning
    for(int j = 0; j < memorydepth; j++) ADC_hist->Fill(denoised[j]);
    baseline_ADC_bin = ADC_hist->GetMaximumBin();
    baseline_ADC = ADC_hist->GetXaxis()->GetBinCenter(baseline_ADC_bin);
    if(waveNum == example_waveform) {
      cout << "baseline_ADC: "<< baseline_ADC <<endl;
      TCanvas *tmp = new TCanvas();
      ADC_hist->Draw();
      tmp->SaveAs("example_ADC_hist.png");
      delete tmp;
    }
    delete ADC_hist;

    // Calculate up-crossing in each waveform
    for(int i=0; i<8; ++i){
      for(int j=0; j<1750; j++){
        if( ( denoised[j] - (baseline_ADC+delta_ADC[i]) )<0 && ( denoised[j+1] - (baseline_ADC+delta_ADC[i]) )>0 ){
          Count0[i] += 1;
        }
      }
    }

    waveNum += 1;
    if (waveNum % 10000 == 0) printf("Wavenum: %d\n", waveNum);

    // Plot single waveform
    if(waveNum == example_waveform) {
      for(int j=0; j<2500; j++){
        h1->Fill(j, raw[j]);
        h2->Fill(j, denoised[j]);
      }
    } // end plot single waveform

  } // End loop over waveforms

  cout<<"==========Result of No Laser:=========== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count0[i]<<endl;
  }

  fin.close();//Very Important------

//======Read v3 laser 1 (20220915; high power)===========
  fin.open("/afs/cern.ch/work/s/shiw/public/ColdBoxVD/Coldbox_Sep2022_leakage_check/20220915_v3_laser1_705_laser2_737_laser3_740/xarapuca_V3_cosmic_trigger50_laser_1/wave5.dat", ios::in | ios::binary);

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
    }
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        denoised[j] = raw[j];
    }

    // Apply noise filter
    if (filternoise) TV1D_denoise(raw, denoised, memorydepth, 50);

    // Evaluate single waveform baseline
    baseline_ADC = 0.;
    baseline_ADC_bin = -1;
    TH1F *ADC_hist = new TH1F("ADC_hist", "ADC_hist", nADCs/4, 0, nADCs); // this bin can't be too fine, 4 adc binning
    for(int j = 0; j < memorydepth; j++) ADC_hist->Fill(denoised[j]);
    baseline_ADC_bin = ADC_hist->GetMaximumBin();
    baseline_ADC = ADC_hist->GetXaxis()->GetBinCenter(baseline_ADC_bin);
    delete ADC_hist;

    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if( ( denoised[j] - (baseline_ADC+delta_ADC[i]) )<0 && ( denoised[j+1] - (baseline_ADC+delta_ADC[i]) )>0 ){
	        Count1[i] += 1;
	    }
        }
    }

  }

  cout<<"=====Result of v3 laser 1 (high power):======== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count1[i]<<endl;
  }

  fin.close();

//======Read v2 laser 2 (low power)========================
  fin.open("/afs/cern.ch/work/s/shiw/public/ColdBoxVD/Coldbox_Sep2022_leakage_check/20220915_v3_laser1_705_laser2_737_laser3_740/xarapuca_V3_cosmic_trigger50_laser_2/wave5.dat", ios::in | ios::binary);

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
    }
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        denoised[j] = raw[j];
    }

    // Apply noise filter
    if (filternoise) TV1D_denoise(raw, denoised, memorydepth, 50);

    // Evaluate single waveform baseline
    baseline_ADC = 0.;
    baseline_ADC_bin = -1;
    TH1F *ADC_hist = new TH1F("ADC_hist", "ADC_hist", nADCs/4, 0, nADCs); // this bin can't be too fine, 4 adc binning
    for(int j = 0; j < memorydepth; j++) ADC_hist->Fill(denoised[j]);
    baseline_ADC_bin = ADC_hist->GetMaximumBin();
    baseline_ADC = ADC_hist->GetXaxis()->GetBinCenter(baseline_ADC_bin);
    delete ADC_hist;

    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if( ( denoised[j] - (baseline_ADC+delta_ADC[i]) )<0 && ( denoised[j+1] - (baseline_ADC+delta_ADC[i]) )>0 ){
	        Count2[i] += 1;
	    }
        }
    }

  }

  cout<<"=====Result of v3 laser 2 (high power):======== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count2[i]<<endl;
  }

  fin.close();

//======Read v3 laser 3 (high power)========================
  fin.open("/afs/cern.ch/work/s/shiw/public/ColdBoxVD/Coldbox_Sep2022_leakage_check/20220915_v3_laser1_705_laser2_737_laser3_740/xarapuca_V3_cosmic_trigger50_laser_3/wave5.dat", ios::in | ios::binary);

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
    }
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        denoised[j] = raw[j];
    }

    // Apply noise filter
    if (filternoise) TV1D_denoise(raw, denoised, memorydepth, 50);

    // Evaluate single waveform baseline
    baseline_ADC = 0.;
    baseline_ADC_bin = -1;
    TH1F *ADC_hist = new TH1F("ADC_hist", "ADC_hist", nADCs/4, 0, nADCs); // this bin can't be too fine, 4 adc binning
    for(int j = 0; j < memorydepth; j++) ADC_hist->Fill(denoised[j]);
    baseline_ADC_bin = ADC_hist->GetMaximumBin();
    baseline_ADC = ADC_hist->GetXaxis()->GetBinCenter(baseline_ADC_bin);
    delete ADC_hist;

    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if( ( denoised[j] - (baseline_ADC+delta_ADC[i]) )<0 && ( denoised[j+1] - (baseline_ADC+delta_ADC[i]) )>0 ){
	        Count3[i] += 1;
	    }
        }
    }

  }

  cout<<"=====Result of v3 laser 3 (high power):======== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count3[i]<<endl;
  }

  fin.close();


//======Read all lasers of v3========================
  fin.open("/afs/cern.ch/work/s/shiw/public/ColdBoxVD/Coldbox_Sep2022_leakage_check/20220915_v3_laser1_705_laser2_737_laser3_740/xarapuca_V3_cosmic_trigger50_laser_1_and_2_and_3/wave5.dat", ios::in | ios::binary);

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
    }
    for(int j = 0; j < memorydepth; j++) {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        denoised[j] = raw[j];
    }

    // Apply noise filter
    if (filternoise) TV1D_denoise(raw, denoised, memorydepth, 50);

    // Evaluate single waveform baseline
    baseline_ADC = 0.;
    baseline_ADC_bin = -1;
    TH1F *ADC_hist = new TH1F("ADC_hist", "ADC_hist", nADCs/4, 0, nADCs); // this bin can't be too fine, 4 adc binning
    for(int j = 0; j < memorydepth; j++) ADC_hist->Fill(denoised[j]);
    baseline_ADC_bin = ADC_hist->GetMaximumBin();
    baseline_ADC = ADC_hist->GetXaxis()->GetBinCenter(baseline_ADC_bin);
    delete ADC_hist;

    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if( ( denoised[j] - (baseline_ADC+delta_ADC[i]) )<0 && ( denoised[j+1] - (baseline_ADC+delta_ADC[i]) )>0 ){
	        Count123[i] += 1;
	    }
        }
    }

  }

  cout<<"====Result of All lasers of v3 (high power):=== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count123[i]<<endl;
  }

  fin.close();


//======Gain calculations==========
  for(int i=0; i<8; ++i){
      Gain1[i] = (Count1[i] - Count0[i])/20000*5000/3500;
      Gain2[i] = (Count2[i] - Count0[i])/20000*5000/3500;
      Gain3[i] = (Count3[i] - Count0[i])/20000*5000/3500;
      Gain123[i] = (Count123[i] - Count0[i])/20000*5000/3500;
  }

//======Drawing====================
  TGraph *dis0 = new TGraph(8, delta_ADC, Count0);
  TGraph *dis1 = new TGraph(8, delta_ADC, Count1);
  TGraph *dis2 = new TGraph(8, delta_ADC, Count2);
  TGraph *dis3 = new TGraph(8, delta_ADC, Count3);
  TGraph *dis123 = new TGraph(8, delta_ADC, Count123);

  TCanvas *c = new TCanvas();
  c->SetLogy(); c->SetGridx(); c->SetGridy();
  TMultiGraph *mg = new TMultiGraph();

  dis0->SetMarkerSize(0.8);
  dis0->SetMarkerStyle(21);
  dis0->SetMarkerColor(1);
  dis0->SetLineColor(1);

  dis1->SetMarkerSize(0.8);
  dis1->SetMarkerStyle(21);
  dis1->SetMarkerColor(2);
  dis1->SetLineColor(2);

  dis2->SetMarkerSize(0.8);
  dis2->SetMarkerStyle(21);
  dis2->SetMarkerColor(3);
  dis2->SetLineColor(3);

  dis3->SetMarkerSize(0.8);
  dis3->SetMarkerStyle(21);
  dis3->SetMarkerColor(41);
  dis3->SetLineColor(41);

  dis123->SetMarkerSize(0.8);
  dis123->SetMarkerStyle(3);
  dis123->SetMarkerColor(6);
  dis123->SetLineColor(6);

  mg->Add(dis0);
  mg->Add(dis1);
  mg->Add(dis2);
  mg->Add(dis3);
  mg->Add(dis123);
  mg->GetXaxis()->SetTitle("#Delta (ADC)");
  mg->GetYaxis()->SetTitle("Counts");
  mg->SetMaximum(Count123[0]*2.);
  mg->SetMinimum(Count0[7]/2.);
  mg->Draw("AP");

  TLegend *leg = new TLegend(0.55, 0.65, 0.9, 0.9);
  leg->AddEntry(dis0, "No Laser: 3350ADC");
  leg->AddEntry(dis1, "v3 Laser 1 [388mW]");
  leg->AddEntry(dis2, "v3 Laser 2 [410mW]");
  leg->AddEntry(dis3, "v3 Laser 3 [250mW]");
  leg->AddEntry(dis123, "v3 (Laser1+2+3) [1048mW]");
  leg->Draw();

  c->SaveAs("up-crossing_count.png");


//======Drawing Gain==================
  TCanvas *c2 = new TCanvas();
  c2->SetGridx(); c2->SetGridy();
  TMultiGraph *mg2 = new TMultiGraph();

  TGraph *diff1 = new TGraph(8, delta_ADC, Gain1);
  TGraph *diff2 = new TGraph(8, delta_ADC, Gain2);
  TGraph *diff3 = new TGraph(8, delta_ADC, Gain3);
  TGraph *diff123 = new TGraph(8, delta_ADC, Gain123);

  diff1->SetMarkerSize(0.8);
  diff1->SetMarkerStyle(21);
  diff1->SetMarkerColor(2);
  diff1->SetLineColor(2);

  diff2->SetMarkerSize(0.8);
  diff2->SetMarkerStyle(21);
  diff2->SetMarkerColor(3);
  diff2->SetLineColor(3);

  diff3->SetMarkerSize(0.8);
  diff3->SetMarkerStyle(21);
  diff3->SetMarkerColor(41);
  diff3->SetLineColor(41);

  diff123->SetMarkerSize(0.8);
  diff123->SetMarkerStyle(3);
  diff123->SetMarkerColor(6);
  diff123->SetLineColor(6);

  mg2->Add(diff1);
  mg2->Add(diff2);
  mg2->Add(diff3);
  mg2->Add(diff123);

  TLegend *leg2 = new TLegend(0.55, 0.65, 0.9, 0.9);
  leg2->AddEntry(diff1, "v2 SwitchA [1458mW]");
  leg2->AddEntry(diff2, "v2 SwitchB [581mW]");
  leg2->AddEntry(diff3, "v2 SwitchA & SwitchB [2039mW]");
  leg2->AddEntry(diff123, "All Lasers of v2 & v3 on [4221mW]");

  mg2->GetXaxis()->SetTitle("#Delta (ADC)");
  mg2->GetYaxis()->SetTitle("Counts per 10#mus");

  mg2->Draw("AP");
  leg2->Draw();
  c2->SaveAs("extra_up-crossing_per10us_count.png");

  TCanvas *c3 = new TCanvas();
  h1->SetLineColor(2);
  h1->Draw("HIST");
  h1->GetXaxis()->SetTitle("Time tick (4ns)");
  h1->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
  h2->SetLineColor(4);
  h2->SetLineStyle(7);
  h2->Draw("HIST SAME");
  c3->SaveAs("example_single_waveform.png");

}
