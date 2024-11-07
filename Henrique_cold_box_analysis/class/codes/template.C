#define memorydepth 1252
#include "/Documents/apc_root/cold_box_analysis/class/MYCODES.h"

void a(){

   vector<Int_t> channels = {1};
  TFile *f1 = new TFile("analyzed.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  vector<ADC_DATA> ch(channels.size());
  vector<TBranch*> bch(channels.size());
  for(Int_t k = 0; k<channels.size();k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
  }
  Int_t nentries = t1->GetEntries();
  for(Int_t i = 0; i<nentries; i++){
    for(Int_t k = 0; k<channels.size();k++){
      bch[k]->GetEvent(i);
      // code here
      // for(Int_t j = 0; j<memorydepth; j++){
      //   cout << ch[k].wvf[j] << endl;
      // }
    }
  }
}
