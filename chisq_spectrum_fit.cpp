#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <math.h>

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27

std::vector<double> refHist = {24.0,276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

typedef struct evt {
    double time;
    double energy;
} evt;

std::vector<double> createHist(double threshold, double cutoff, double height, std::vector<evt>& events) {
    std::vector<double> hist;
    hist.resize(185, 0.0);
    for(auto it = events.begin(); it < events.end(); it++) {
        double weight = 0.0;
        if(it->time < 40 || it->time > 225) {
            continue;
        }
        if(it->energy*JTONEV < threshold) {
            continue;
        }
        if(it->energy*JTONEV >= threshold && it->energy*JTONEV < cutoff) {
            weight = height/(cutoff - threshold)*(it->energy*JTONEV - threshold);
//            weight = it->energy/(GRAV*MASS_N)/0.1;
        }
        if(it->energy*JTONEV >= cutoff) {
             weight = height + (1.0 - height)/(0.345 - cutoff)*(it->energy*JTONEV - cutoff);
        }
        if(weight > 1.0) {
            printf("Boo!\n");
        }
        hist[int(it->time)-40] += weight;
    }
    return hist;
}

//std::vector<double> createHistSqrt(std::vector<evt>& events) {
//    std::vector<double> hist;
//    hist.resize(185, 0.0);
//    for(auto it = events.begin(); it < events.end(); it++) {
//        double weight = sqrt(it->energy*JTONEV/0.345);
//        hist[int(it->time)-40] += weight;
//    }
//    return hist;
//}

double calcChisq(std::vector<double>& hn, std::vector<double>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double cn = std::accumulate(hm.begin(), hm.end(), 0);
    double cm = std::accumulate(hn.begin(), hn.end(), 0);
//    printf("%f %f\n", cn, cm);
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        chisqSum += (1/(cn*cm))*pow(cn*hn[i] - cm*hm[i],2)/(hn[i]+hm[i]);
    }
    return chisqSum;
}
//def chisq(hn, hm):
//    if len(hn) != len(hm):
//        raise ValueError('Mismatch between histogram lengths')
//    cn = np.sum(hm)
//    cm = np.sum(hn)
//    chisqArray = (1/(cn*cm))*((cn*hn-cm*hm)**2)/(hn+hm)
//    return np.sum(chisqArray))

int main(int argc, char** argv) {
    //buff_len = 4 + 3*8 + 4
    const size_t buff_len = 4 + 3*8 + 4;
    char* buf = new char[buff_len];
    if(argc != 2) {
        printf("Error! Usage: ./chisq_spectrum_fit fname\n");
    }

    std::ifstream binfile(argv[1], std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        printf("Error! Could not open file %s\n", argv[1]);
    }
    
    std::vector<evt> events;
    evt event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        event.time = *((double *)(&buf[0] + sizeof(unsigned int)));
        event.energy =  *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
        events.push_back(event);
    }
    
//    std::vector<double> hist1 = createHist(0.021*GRAV*MASS_N*JTONEV, 0.1*GRAV*MASS_N*JTONEV, 1.0, events);
//    std::vector<double> hist1 = createHist(15.377918, 17.428307, 1.0, events); //"Optimal"
    std::vector<double> hist1 = createHist(0.0, 0.0, 1.0, events);
    printf("%f\n", calcChisq(hist1, refHist)/hist1.size());
    for(auto it = hist1.begin(); it < hist1.end(); it++) {
        printf("%f,", *it);
    }
    printf("\n");
    
//    for(int i = 0; i < 21; i++) {
//        for(int j = 0; j < 21; j++) {
////            for(int k = 0; k < 11; k++) {
//                std::vector<double> hist1 = createHist(GRAV*MASS_N*(0.01+0.2*i/20.0)*JTONEV, GRAV*MASS_N*(0.1+0.2*j/20.0)*JTONEV, 1.0, events);
//                double chisq = calcChisq(hist1, refHist);
////                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
//                printf("%f %f %f\n", GRAV*MASS_N*(0.01+0.2*i/20.0)*JTONEV, GRAV*MASS_N*(0.1+0.2*j/20.0)*JTONEV, chisq/hist1.size());
////            }
//        }
//    }
}