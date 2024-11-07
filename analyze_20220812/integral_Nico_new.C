#include <fstream>
#include <iostream>

void integral_Nico_new(){
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

  TH1F *intPeak = new TH1F("intPeak", "Integral of signal Peak",600,0,15000);
  TH2D *amp_int = new TH2D("amp_int", "Amplitude and Integral", 240, 0, 120, 400, 0, 10000);

  int maxADC = 0;
  int num_1_2 = 0;
  int sumBkg = 0;//to integral over bkg window---
  int sumSig = 0;
  int baseline = 0;
  int peak_height = 0;
  int peakNum = 0;


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
    }
    
    for(int i=2550; i<4000; i++){
        if(raw[i] > maxADC){
	    maxADC = raw[i];
	}
    }//used to remove peak3---   

    waveNum += 1;
//    cout<<"Waveform: "<<waveNum<<endl;

    if(maxADC<6580){// to omit waveform of type3---
	num_1_2 += 1;
	
	for(int i=1100; i<2000; i++){
	    sumBkg += raw[i];
	}
	for(int i=2400; i<2550; i++){
	    if(raw[i] > peak_height){
		peak_height = raw[i];
	    }//to locate the maximal ADC---
	    sumSig += raw[i];
	} 
        sumSig = sumSig - sumBkg/6.0;
	intPeak->Fill(sumSig);//Fill the peak integral

        //to output waveforms with integral larger than 4900---
        if(sumSig<4900 && sumSig>4000){
	    cout<<"Waveform larger than 4900: "<<waveNum<<endl;
	    peakNum += 1;
	}	

	baseline = sumBkg/900.0;
	peak_height = peak_height - baseline;
	amp_int->Fill(peak_height, sumSig);


	sumSig = 0;
	sumBkg = 0;
	peak_height = 0;
	baseline = 0;

    }


    maxADC = 0;
  }

  cout<<"========================================"<<endl;
  cout<<"Num of total waveform: "<<waveNum<<endl;
  cout<<"Num of peak 2: "<<peakNum<<endl;
  cout<<"========================================"<<endl;


  TCanvas *c1 = new TCanvas();
//  gPad->SetLogy();
  intPeak->Draw("colz");
  intPeak->GetXaxis()->SetTitle("[ADC*ns]");
  intPeak->GetYaxis()->SetTitle("Counts");
  intPeak->SetLineColor(kRed);
//  c1->SaveAs("integral_Nico.png");

  TCanvas *c2 = new TCanvas();
  amp_int->Draw("colz");
  amp_int->GetXaxis()->SetTitle("Amplitude[ADC]");
  amp_int->GetYaxis()->SetTitle("Integral of Peak[ADC*ns]");



}



