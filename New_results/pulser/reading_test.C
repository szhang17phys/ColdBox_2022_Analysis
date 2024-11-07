#include <fstream>
#include <iostream>

void reading_test(){
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
  TH2D *h; // I don't really know at this point
  
  ifstream fin;
  fin.open("/afs/cern.ch/user/s/shuaixia/private/PDS_Data/pulser/wave1.dat", ios::in | ios::binary);

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
      h = new TH2D("h","h",memorydepth,0,memorydepth*dtime,nADCs,0,nADCs); //if you are brave, don't divide by 5 so you have high definition plot :)
    }
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, nbytes_data); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        h->Fill(j*dtime,raw[j]);
        //printf("%d %.0f \n",j,raw[j]);
      }

    waveNum += 1;
    printf("Num of waveforms: %d\n", waveNum);

  }

  TCanvas *c1 = new TCanvas();
  h->Draw("colz");
  h->GetXaxis()->SetTitle("Time (ns)");
  h->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  c1->SaveAs("wave6.png");

}
