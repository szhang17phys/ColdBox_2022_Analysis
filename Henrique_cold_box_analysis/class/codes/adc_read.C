#include "/Documents/apc_root/cold_box_analysis/class/MYCODES.h"

void adc_read(){
    
    Read r;
    
    r.dtime = 4;
    r.nbits = 14;
    
    r.baselineTime = 5500; // time limit to start looking for baseline
    r.chargeTime = 4*5000; // last time to integrate 
    
    Int_t linhas = 5000;

    r.OnlyOneEvent = false;
    r.stopEvent = 100;
    
    r.adc_read("0_wave0_36V70_120ADC_Ch0_24",linhas);
//     r.adc_read("0_wave1_36V70_130ADC_Ch0_new_amp",linhas);
    
}

