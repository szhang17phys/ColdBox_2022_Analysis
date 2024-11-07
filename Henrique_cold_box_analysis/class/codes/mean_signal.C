#include "/Documents/apc_root/cold_box_analysis/class/MYCODES.h"

void mean_signal(){
    
    MeanSignal MS;
    
    MS.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    MS.minval = 400000; //charge minimal and maximal values
    MS.maxval = 550000;
    MS.avoid_saturation = 2800; // to avoid peak saturation
    MS.fprompt = 0.5; // well... fprompt
    MS.mustbe = "bigger"; // must be 'bigger' or 'lower' then fprompt set.
    
    MS.mean_signal("analyzed"); 
}
