#define memorydepth 5000
#include "/home/henrique/Dropbox/Unicamp/Doutorado/Root/Programs/italy/ADC_LAr_SuperCell/class/MYCODES.h"



void takeResolution(){
    
    Resolution res;
    
    res.dtime = 4.;
    res.channels = {1,2};
    
    res.trackResolution("wave_35V70_200ADC");
    
}
