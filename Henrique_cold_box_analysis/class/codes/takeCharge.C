#include "/Documents/apc_root/cold_box_analysis/class/MYCODES.h"


void takeCharge(){
    
    TakeCharge C;
    
    
    C.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    C.channels = {1};
    C.integration_start = 6000;
    C.integration_end = 6200;
    C.hmin = -40000;
    C.hmax = 600000;
    C.nbins = 2000;
    
    
    C.takeCharge("analyzed");
    
}
