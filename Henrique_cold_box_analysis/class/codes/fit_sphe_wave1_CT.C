#include "/Documents/apc_root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1_CT(){
  CT_Calibration Cal;
  
  Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
  
  Cal.rebin = 1;
  
  Cal.n_peaks = 7;
  n_CT = 6;
  Cal.peak0 = 2800;
  Cal.mean0 = 0;
  Cal.sigma0 = 511;
  
  Cal.peak1 = 1500;
  Cal.mean1 = 13000;
  Cal.sigma1 = 700;
  
  Cal.startpoint = 500;
  
  Cal.xmin = -10000;
  Cal.xmax = 300000;
  
  
  Cal.startFirst = 5000;
  Cal.finishFirst = 26000;
  
  Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
  Cal.fit_sphe_wave("analyzed_1");  

    

}
