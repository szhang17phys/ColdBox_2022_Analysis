#include "/home/henrique/Dropbox/Unicamp/Doutorado/Root/Programs/italy/ADC_LAr_SuperCell/class/MYCODES.h"

void fit_sphe_wave2(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 4;
    
    Cal.n_peaks = 3; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 = 0.05; // amplitude baseline
    Cal.mean0 = 100; // average baseline
    Cal.sigma0 = 100; // sigma baseline
     
    Cal.peak1 = 0.1; // same for 1 p.e.
    Cal.mean1 = 1800;
    Cal.sigma1 = 200;
    
    Cal.startpoint = 0.015; // amplitude for 2 p.e. 
        
    Cal.xmin = -10000; // range for graph display (not fit)
    Cal.xmax = 80000;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch2.root";
    Cal.fit_sphe_wave("analyzed_2");    
}
