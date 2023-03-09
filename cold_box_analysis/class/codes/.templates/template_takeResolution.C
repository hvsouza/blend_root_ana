#define memorydepth 5000
#include "__USER_PATH__/MYCODES.h"



void takeResolution(){
    
    Resolution res;
    
    res.dtime = 4.;
    res.channels = {1,2};
    
    res.trackResolution("wave_35V70_200ADC");
    
}
