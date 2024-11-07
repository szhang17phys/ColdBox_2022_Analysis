
class ADC_DATA{
public:
  Double_t peak;
  Double_t peakpos;
  Double_t charge;
  Double_t fprompt;
  Double_t event;
  Double_t time;
  Double_t wvf[memorydepth];
  Int_t selection;
    
  ADC_DATA();
  const char *tobranch = Form("peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[%i]/D:selection/O",memorydepth);
};


ADC_DATA::ADC_DATA(){
  
  for(Int_t i = 0; i<memorydepth; i++){
    wvf[i]=0;
  }
  peak = 0;
  peakpos = 0;
  charge = 0;
  fprompt = 0;
  event = 0; 
  time = 0;
  selection = 0; 
}
  
