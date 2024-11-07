#include <fstream>
#include <iostream>

#include "filter.h"

//This code will show overlaid figure of initial waveform 
//as well as waveform after filter!---
//Shuaixiang Zhang; Sep 27, 2022---

void filter_effect(){
  int headbin; // to store headers
  int nbytes_headers = 4; // 4 bytes (32 bits) for each head 
  int nbytes_data = 2; // 2 bytes (16 bits) per sample
  int memorydepth = 0; // size of waveforms 
  uint32_t valbin = 0; // to read data

  float time[2500] = {0};//to record the x axis---
  float raw[2500] = {0}; // waveform as array---
  float filter_data[2500] = {0};//record data after filter---

  int waveNum = 0;//To record num of waveforms; Added by szhang; 20220720---

  int select14 = 19070;

  int y1 = 3100;//range of TH2D---
  int y2 = 3700;



  Double_t nbits = 14; // ADC is a 14 bits, 2 Vpp
  Double_t samplingRate = 500.e6; // 250 MSamples/s for DT5725
  Double_t nADCs = pow(2,nbits); // number of digital channels
  Double_t inVolts = 2./nADCs; // Multiply by this number to have the amplitude in volts;
  Double_t dtime = (2/samplingRate)*1e9; // steps in nanoseconds

  
  Bool_t first_line = true;// so we can set the length of the vector---
  
  ifstream fin;
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v2_and_v3_all_lasers/wave5.dat", ios::in | ios::binary);

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
//      raw.resize(memorydepth);
      first_line=false;
      
    }

    waveNum += 1;
    printf("Num of waveforms: %d\n", waveNum);


    for(int j = 0; j < memorydepth; j++){
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
            break;
        }
	time[j] = 4.0*(j+1);
        raw[j] = valbin;
    }

    
    if(waveNum == select14){
	TV1D_denoise(raw, filter_data, memorydepth, 50.0);
        break;
    }
    
 
  }


//===Plotting====================
    TGraph *graph_raw = new TGraph(2500, time, raw);
    TGraph *graph_filter = new TGraph(2500, time, filter_data);

    TCanvas *c14 = new TCanvas();

    graph_raw->GetXaxis()->SetTitle("Time [ns]");
    graph_raw->GetYaxis()->SetTitle("Amplitude [ADC Counts]");   
    graph_raw->GetXaxis()->SetRangeUser(0, 10000);
    graph_raw->GetYaxis()->SetRangeUser(3100, 4000);
    graph_raw->SetLineColor(1);
    graph_filter->SetLineColor(2);
    graph_filter->SetLineWidth(2);

    graph_raw->Draw();
    graph_filter->Draw("same");

    TLegend *leg = new TLegend(0.4, 0.75, 0.6, 0.9);
    leg->AddEntry(graph_raw, "Raw Waveform");
    leg->AddEntry(graph_filter, "Filtered Waveform");
    leg->Draw();
  
//  c1->SaveAs("wave1.png");

}
