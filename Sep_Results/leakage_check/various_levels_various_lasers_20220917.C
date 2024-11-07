#include <fstream>
#include <iostream>

void various_levels_various_lasers_20220917(){
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
  
  Bool_t first_line = true; // so we can set the length of the vector


//======General Variables====================  
  double level[8] = {0};
  double Count0[8] = {0};
  double Count1[8] = {0};
  double Count2[8] = {0};
  double Count3[8] = {0};
  double Count123[8] = {0};

  double Gain1[8] = {0};
  double Gain2[8] = {0};
  double Gain3[8] = {0};
  double Gain123[8] = {0};


//  level[0] = 3385;//for wave5.dat---
  level[0] = 2595;//for wave6.dat, baseline is 2545---
  for(int i=1; i<8; ++i){
      level[i] = level[0]+25*i;
  }


//======Read No Laser========================
  ifstream fin;
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/argon_only/wave6.dat", ios::in | ios::binary);

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
//        h->Fill(j*dtime,raw[j]);
        //printf("%d %.0f \n",j,raw[j]);
    }
    
    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if((raw[j] - level[i])<0 && (raw[j+1] - level[i])>0){
	        Count0[i] += 1;
	    }
        }   
    }  


    waveNum += 1;
    printf("Wavenum: %d\n", waveNum);

  }

  cout<<"==========Result of No Laser:=========== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count0[i]<<endl;
  }

  fin.close();//Very Important------

//======Read Laser 1========================
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_1_and_argon_read/wave6.dat", ios::in | ios::binary);

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
    }
    
    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if((raw[j] - level[i])<0 && (raw[j+1] - level[i])>0){
	        Count1[i] += 1;
	    }
        }   
    }  

  }

  cout<<"==========Result of Laser 1:=========== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count1[i]<<endl;
  }

  fin.close();

//======Read Laser 2========================
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_2_and_argon_read/wave6.dat", ios::in | ios::binary);

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
    }
    
    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if((raw[j] - level[i])<0 && (raw[j+1] - level[i])>0){
	        Count2[i] += 1;
	    }
        }   
    }  

  }

  cout<<"==========Result of Laser 2:=========== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count2[i]<<endl;
  }

  fin.close();

//======Read Laser 3========================
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_3_and_argon_read/wave6.dat", ios::in | ios::binary);

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
    }
    
    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if((raw[j] - level[i])<0 && (raw[j+1] - level[i])>0){
	        Count3[i] += 1;
	    }
        }   
    }  

  }

  cout<<"==========Result of Laser 3:=========== "<<endl;
  for(int i=0; i<8; ++i){
      cout<<"Counts of Level "<<level[i]<<": "<<Count3[i]<<endl;
  }

  fin.close();


//======Read Laser 1 and 2 and 3========================
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_123_and_argon_read/wave6.dat", ios::in | ios::binary);

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
    }
    
    for(int i=0; i<8; ++i){
        for(int j=0; j<1750; j++){
            if((raw[j] - level[i])<0 && (raw[j+1] - level[i])>0){
	        Count123[i] += 1;
	    }
        }   
    }  

  }

  cout<<"==========Result of Laser 1 & 2 & 3:=========== "<<endl;
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

  for(int i=0; i<8; ++i){
      level[i] = 25*(i+2);
  }
//======Drawing====================
  TGraph *dis0 = new TGraph(8, level, Count0);
  TGraph *dis1 = new TGraph(8, level, Count1);
  TGraph *dis2 = new TGraph(8, level, Count2);
  TGraph *dis3 = new TGraph(8, level, Count3);
  TGraph *dis123 = new TGraph(8, level, Count123);

  TCanvas *c = new TCanvas();
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

  TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.9);
  leg->AddEntry(dis0, "No Laser: 2545ADC");
  leg->AddEntry(dis1, "v3 Laser 1 [388mW]");
  leg->AddEntry(dis2, "v3 Laser 2 [410mW]");
  leg->AddEntry(dis3, "v3 Laser 3 [250mW]");
  leg->AddEntry(dis123, "v3 (Laser1+2+3) [1048mW]");

  mg->GetXaxis()->SetTitle("#Delta (ADC)");
  mg->GetYaxis()->SetTitle("Counts");

  mg->Draw("AP");
  leg->Draw();


//======Drawing Gain==================
  TCanvas *c2 = new TCanvas();
  TMultiGraph *mg2 = new TMultiGraph();

  TGraph *diff1 = new TGraph(8, level, Gain1);
  TGraph *diff2 = new TGraph(8, level, Gain2);
  TGraph *diff3 = new TGraph(8, level, Gain3);
  TGraph *diff123 = new TGraph(8, level, Gain123);

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

  TLegend *leg2 = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg2->AddEntry(diff1, "v3 Laser 1 [388mW]");
  leg2->AddEntry(diff2, "v3 Laser 2 [410mW]");
  leg2->AddEntry(diff3, "v3 Laser 3 [250mW]");
  leg2->AddEntry(diff123, "v3 (Laser1+2+3) [1048mW]");

  mg2->GetXaxis()->SetTitle("#Delta (ADC)");
  mg2->GetYaxis()->SetTitle("Counts per 10#mus");

  mg2->Draw("AP");
  leg2->Draw();


//  c1->SaveAs("disADC6.png");

}
