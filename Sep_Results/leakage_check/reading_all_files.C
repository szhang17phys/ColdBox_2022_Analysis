#include <fstream>
#include <iostream>

void reading_all_files(){
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

  //here each points last for 4ns!
  Double_t dtime = (2/samplingRate)*1e9; // steps in nanoseconds
  Bool_t first_line = true; // so we can set the length of the vector
  TH2D *h; // I don't really know at this point
  TH2D *h1;



//======v3 20220917 argon only==================  
  ifstream fin;
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/argon_only/wave5.dat", ios::in | ios::binary);

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
      h = new TH2D("h","h",memorydepth,0,memorydepth*dtime,nADCs,1000,16000);
      h1 = new TH2D("h1","h1",memorydepth,0,memorydepth*dtime,40000,3100,3600);
    }
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }

//    waveNum += 1;
//    printf("Num of waveforms: %d\n", waveNum);

  }

  fin.close();
  cout<<"======v3 20220917 argon only Read!======"<<endl;




//======v3 20220917 laser 1 on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_1_and_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }

//    waveNum += 1;
//    printf("Num of waveforms: %d\n", waveNum);

  }

  fin.close();
  cout<<"======v3 20220917 laser 1 Read!======"<<endl;



//======v3 20220917 laser 2 on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_2_and_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v3 20220917 laser 2 Read!======"<<endl;




//======v3 20220917 laser 3 on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_3_and_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v3 20220917 laser 3 Read!======"<<endl;






//======v3 20220917 laser 1&2&3 on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v3_laser1_388_laser2_410_laser3_250/cathode_v3_laser_123_and_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v3 20220917 laser 123 Read!======"<<endl;





//======v2 20220917 switchA on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v2_swtichA_1458_switchB_581/cathode_v2_switchA_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v2 20220917 switchA Read!======"<<endl;




//======v2 20220917 switchB on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v2_swtichA_1458_switchB_581/cathode_v2_switchB_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v2 20220917 switchB Read!======"<<endl;




//======v2 20220917 switchA & switchB on==================  
  fin.open("/afs/cern.ch/work/s/shuaixia/public/Coldbox_Data_2022/Sep_Test/light_leakage_check/20220917_v2_swtichA_1458_switchB_581/cathode_v2_switchA_and_switchB_argon_read/wave5.dat", ios::in | ios::binary);

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
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======v2 20220917 switchA&switchB Read!======"<<endl;





//======v2&v3 20220917==================  
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
      raw.resize(memorydepth);
      first_line=false;
    }
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        h1->Fill(j*dtime,raw[j]);
      }
  }

  fin.close();
  cout<<"======All laser of v2&v3 20220917 Read!======"<<endl;




//===Drawing==============================

  TCanvas *c1 = new TCanvas();
  h->Draw("colz");
  h->GetXaxis()->SetTitle("Time (ns)");
  h->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  TCanvas *c2 = new TCanvas();
  h1->Draw("colz");
  h1->GetXaxis()->SetTitle("Time (ns)");
  h1->GetYaxis()->SetTitle("Amplitude (ADC Channels)");


}
