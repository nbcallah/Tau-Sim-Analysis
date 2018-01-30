#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <complex>

extern "C" {
    #include "xorshift.h"
}

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34

#define NBORON 1.37e29
#define ABORON -0.1e-15
#define SIGMABORON 2200*3.835e-25

#define NCARBON 1.133e29
#define ACARBON 6.6460e-15
#define SIGMACARBON 5.551e-28

#define NZINC 2.527e28
#define AZINC 5.68e-15
#define SIGMAZINC 5.241e-28

#define NSULFUR 2.527e28
#define ASULFUR 2.847e-15
#define SIGMASULFUR 1.556e-28

//std::vector<double> refHist = {24.0,276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};
std::vector<double> refHist = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

typedef struct evt {
    double energy;
    float times[50];
    float ePerp[50];
} evt;

typedef struct weightedBin {
    double wgt;
    double wgtSqr;
    double num;
} weightedBin;

bool absorb(double ePerp, double u, double fracCarb) {
//    printf("%f\n", contam);
    double v = 2*M_PI*HBAR*HBAR*(NCARBON*fracCarb)*ACARBON/MASS_N + 2*M_PI*HBAR*HBAR*(NBORON*(1-fracCarb))*ABORON/MASS_N;
    double w = HBAR*(NBORON*(1-fracCarb))*SIGMABORON/2;
    double alpha = sqrt((v-ePerp)*(v-ePerp) + w*w);
    double prob = 1 - 
        (ePerp - sqrt(ePerp)*sqrt(2*alpha - 2*(v-ePerp)) + alpha)
        /
        (ePerp + sqrt(ePerp)*sqrt(2*alpha - 2*(v-ePerp)) + alpha);
//    printf("%e %e %e %e %e\n", ePerp, prob, v, w, alpha);
    return u < prob ? true : false;
}

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b) {
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0)};
    if(a.size() != 4 || b.size() != 4) {
        return res;
    }
    res[0] = a[0]*b[0] + a[1]*b[2];
    res[1] = a[0]*b[1] + a[1]*b[3];
    res[2] = a[2]*b[0] + a[3]*b[2];
    res[3] = a[2]*b[1] + a[3]*b[3];
    return res;
}

std::complex<double> k(double ePerp, std::complex<double> u) {
    return std::sqrt((2*MASS_N/(HBAR*HBAR))*(ePerp - u));
}

std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1) {
    return knm1/kn;
}

std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z) {
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0)};
    res[0] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1-kn)*z);
    res[1] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1+kn)*z);
    res[2] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1+kn)*z);
    res[3] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1-kn)*z);
    return res;
}

bool absorbMultilayer(double ePerp, double u, double thick) {
    const double vcarbon = (2*M_PI*(HBAR*HBAR)/MASS_N)*ACARBON*NCARBON;
    const double wcarbon = (HBAR/2)*NCARBON*SIGMACARBON;
    const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
    const double wboron = (HBAR/2)*NBORON*SIGMABORON;
    const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
    const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
    
    std::vector<std::complex<double>> pots = {std::complex<double>(0, 0),
                                              std::complex<double>(vcarbon, -wcarbon),
                                              std::complex<double>(vboron, -wboron),
                                              std::complex<double>(vzns, -wzns)};
    std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
    std::vector<double> zs = {0.0, thick*1e-9, thick*1e-9 + 20e-9, 10000e-9};
    
    for(int i = pots.size(); i > 0; i--) {
        mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
    }
    
    double refl = (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();    
    if(u < (1-refl)) {
        return true;
    }
    return false; 
}

std::vector<weightedBin> createHistQuant(double fracCarb, double threshold, double saturation, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        if(events[i].energy*JTONEV >= threshold && events[i].energy*JTONEV < saturation) {
            weight = 1.0/(saturation - threshold)*(events[i].energy*JTONEV - threshold);
        }
        else {
             weight = 1.0;
        }

        for(int j = 0; j < 50; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
            if(absorb(events[i].ePerp[j], randU01s[i*100 + j], fracCarb)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                break;
            }
        }
    }
    return hist;
}

std::vector<weightedBin> createHistQuantEdE(double fracCarb, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        weight = events[i].energy*JTONEV/34.5;

        for(int j = 0; j < 50; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
            if(absorb(events[i].ePerp[j], randU01s[i*100 + j], fracCarb)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                break;
            }
        }
    }
    return hist;
}

std::vector<weightedBin> createHistQuantMultilayerEdE(double thick, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        weight = events[i].energy*JTONEV/34.5;

        for(int j = 0; j < 50; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thick)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                break;
            }
        }
    }
    return hist;
}

std::vector<weightedBin> createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        if(events[i].energy*JTONEV >= threshold && events[i].energy*JTONEV < saturation) {
            weight = 1.0/(saturation - threshold)*(events[i].energy*JTONEV - threshold);
        }
        else {
             weight = 1.0;
        }

        for(int j = 0; j < 50; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
            if(events[i].ePerp[j]*JTONEV > absCut && randU01s[i*100 + j] < absProb) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                break;
            }
        }
    }
    return hist;
}

double calcChisq(std::vector<double>& hn, std::vector<weightedBin>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double cn = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double cm = std::accumulate(hn.begin(), hn.end(), 0);
//    printf("%f %f\n", cn, cm);
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        chisqSum += (1/(cn*cm))*pow(cn*hn[i] - cm*hm[i].wgt,2)/(hn[i]+hm[i].wgt);
    }
    return chisqSum;
}

double calcChisqWgt(std::vector<double>& hn, std::vector<weightedBin>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double cn = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double cm = std::accumulate(hn.begin(), hn.end(), 0);
    
//    printf("cn: %f cm:%f\n", cn, cm)
    
//    double cntilde = cn;
//    double cmtilde = cm
//                    *std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt*bin.wgt;})
//                    /std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
        
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        double cntilde = cn;
        double cmtilde = cm*hm[i].wgtSqr/hm[i].wgt;
        chisqSum += (1/(cntilde*cmtilde))*pow(cntilde*hn[i] - cmtilde*hm[i].num,2)/(hn[i]+hm[i].num);
    }
    return chisqSum;
}

double calcChisqNate(std::vector<double>& hn, std::vector<weightedBin>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double sumN = std::accumulate(hn.begin(), hn.end(), 0);
    double sumM = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        chisqSum += pow(hn[i]/sumN - hm[i].wgt/sumM,2)/(hn[i]/(sumN*sumN) + (hm[i].wgt*hm[i].wgt)/(sumM*sumM*hm[i].num));
    }
    return chisqSum;
}

void writeYOD(std::vector<weightedBin>& hm, char* oFile) { //function to write yoda histograms for professor minimizer
    std::ofstream out;
    out.open(oFile, std::ios::out);
    out << std::scientific;
    
    double sumW = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double sumW2 = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgtSqr;});
    double num = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.num;});
    double sumWx = 0.0;
    double sumWx2 = 0.0;
    for(int i = 0; i < hm.size(); i++) {
        sumWx += i*hm[i].wgt/sumW;
        sumWx2 += i*i*hm[i].wgt/sumW;
//        sumWx += i*hm[i].wgt;
//        sumWx2 += i*i*hm[i].wgt;
    }
    
    out << "BEGIN YODA_HISTO1D_V2 /\n"
            "Path: /\n";
    out << "ScaledBy: " << 1.0/sumW << "\n";
    out << "Title: \n"
            "Type: Histo1D\n"
            "---\n"
//            "# Mean: 1.000000e+00\n"
//            "# Area: 1.000000e+00\n"
            "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    out << "Total\tTotal\t1.000000e+00\t" << sumW2/(sumW*sumW) << "\t" << sumWx << "\t" << sumWx2 << "\t" << num << "\n";
//    out << "Total\tTotal\t" << sumW << "\t" << sumW2 << "\t" << sumWx << "\t" << sumWx2 << "\t" << num << "\n";
    out << "Underflow\tUnderflow\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\n"
            "Overflow\tOverflow\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\n"
            "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    for(int i = 0; i < hm.size(); i++) {
        out << (double)i << "\t" << (double)(i+1) << "\t" << hm[i].wgt/sumW << "\t" << hm[i].wgtSqr/(sumW*sumW) << "\t" << i*hm[i].wgt/sumW << "\t" << i*i*hm[i].wgt/sumW << "\t" << hm[i].num << "\n";
//        out << (double)i << "\t" << (double)(i+1) << "\t" << hm[i].wgt << "\t" << hm[i].wgtSqr << "\t" << i*hm[i].wgt << "\t" << i*i*hm[i].wgt << "\t" << hm[i].num << "\n";
    }
    out << "END YODA_HISTO1D\n";
    out.close();
}

void writeYODDouble(std::vector<double>& hm, char* oFile) { //function to write yoda histograms for professor minimizer
    std::ofstream out;
    out.open(oFile, std::ios::out);
    out << std::scientific;

    double sumW = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, double bin)->double{return sum + bin;});
    double sumW2 = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, double bin)->double{return sum + bin;});
    double num = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, double bin)->double{return sum + bin;});
    double sumWx = 0.0;
    double sumWx2 = 0.0;
    for(int i = 0; i < hm.size(); i++) {
        sumWx += i*hm[i]/sumW;
        sumWx2 += i*i*hm[i]/sumW;
//        sumWx += i*hm[i];
//        sumWx2 += i*i*hm[i];
    }

    out << "BEGIN YODA_HISTO1D_V2 /\n"
            "Path: /\n";
    out << "ScaledBy: " << 1.0/sumW << "\n";
    out << "Title: \n"
            "Type: Histo1D\n"
            "---\n"
//            "# Mean: 1.000000e+00\n"
//            "# Area: 1.000000e+00\n"
            "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
//    out << "Total\tTotal\t" << sumW << "\t" << sumW2 << "\t" << sumWx << "\t" << sumWx2 << "\t" << num << "\n";
    out << "Total\tTotal\t1.000000e+00\t" << sumW2/(sumW*sumW) << "\t" << sumWx << "\t" << sumWx2 << "\t" << num << "\n";
    out << "Underflow\tUnderflow\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\n"
            "Overflow\tOverflow\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00\n"
            "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    for(int i = 0; i < hm.size(); i++) {
        out << (double)i << "\t" << (double)(i+1) << "\t" << hm[i]/sumW << "\t" << hm[i]/(sumW*sumW) << "\t" << i*hm[i]/sumW << "\t" << i*i*hm[i]/sumW << "\t" << hm[i] << "\n";
//        out << (double)i << "\t" << (double)(i+1) << "\t" << hm[i] << "\t" << hm[i] << "\t" << i*hm[i] << "\t" << i*i*hm[i] << "\t" << hm[i] << "\n";
    }
    out << "END YODA_HISTO1D\n";
    out.close();
}

int main(int argc, char** argv) {
    initxorshift();
    //buff_len = 4 + 3*8 + 4
    const size_t buff_len = 4 + 1*8 + 100*4 + 4;
    char* buf = new char[buff_len];
    if(argc != 2) {
        printf("Error! Usage: ./chisq_spectrum_fit fname\n");
    }

    std::ifstream binfile(argv[1], std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        printf("Error! Could not open file %s\n", argv[1]);
        return 1;
    }
    
    std::vector<evt> events;
    evt event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double)), 50*sizeof(float));
        std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double) + 50*sizeof(float)), 50*sizeof(float));
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events!\n", events.size());
    
    double* randU01s = new double[events.size()*100];
    double* randDeathTimes = new double[events.size()];
//    randU01s = (double *)malloc(sizeof(double)*events.size()*100);
//    randU01s = (double *)malloc(sizeof(double)*events.size()*100);
    
    for(unsigned long i = 0; i < events.size()*100; i++) {
        randU01s[i] = nextU01();
    }
    for(unsigned long i = 0; i < events.size(); i++) {
        randDeathTimes[i] = -877.7*log(nextU01());
    }
    
    
/*    int nBins = 20;
    #pragma omp parallel for collapse(4)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            for(int k = 0; k < nBins+1; k++) {
                for(int l = 0; l < nBins+1; l++) {
//                    double absProb = 0.15 + 0.15*k/(double)nBins;
//                    double absCut = 0.5+3*l/(double)nBins;
//                    double thresh = GRAV*MASS_N*(0.01+0.10*i/(double)nBins)*JTONEV;
//                    double sat = 20+14.5*j/(double)nBins;
                    double absProb = 0.22 + 0.12*k/(double)nBins;
                    double absCut = 1.5+2*l/(double)nBins;
                    double thresh = 2.0+4.0*i/(double)nBins;
                    double sat = 24+10.5*j/(double)nBins;
                    std::vector<weightedBin> hist1 = createHist(absProb, absCut, thresh, sat, events, randU01s, randDeathTimes);
                    double chisqWgt = calcChisqWgt(refHist, hist1);
                    double chisqUnWgt = calcChisq(refHist, hist1);
                    double chisqNate = calcChisqNate(refHist, hist1);
    //                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
                    printf("%f %f %f %f %f %f %f\n", absProb, absCut, thresh, sat, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
                    fflush(stdout);
                }
            }
        }
    }*/
    
    /*int nBins = 20;
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            for(int k = 0; k < nBins+1; k++) {
                double thresh = 2.0 + 4.0*i/(double)nBins;
                double sat = 24 + 10.5*j/(double)nBins;
                double contam = 0.31 + 0.04*k/(double)nBins;
                std::vector<weightedBin> hist1 = createHistQuant(contam, thresh, sat, events, randU01s, randDeathTimes);
                double chisqWgt = calcChisqWgt(refHist, hist1);
                double chisqUnWgt = calcChisq(refHist, hist1);
                double chisqNate = calcChisqNate(refHist, hist1);
//                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
                printf("%f %f %f %f %f %f\n", contam, thresh, sat, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
                fflush(stdout);
            }
        }
    }*/
    
    /*int nBins = 20;
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            double thresh = 2.0 + 4.0*i/(double)nBins;
            double contam = 0.31 + 0.04*j/(double)nBins;
            std::vector<weightedBin> hist1 = createHistQuantEdE(contam, thresh, events, randU01s, randDeathTimes);
            double chisqWgt = calcChisqWgt(refHist, hist1);
            double chisqUnWgt = calcChisq(refHist, hist1);
            double chisqNate = calcChisqNate(refHist, hist1);
//                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
            printf("%f %f %f %f %f\n", contam, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
            fflush(stdout);
        }
    }*/
    
    int nBins = 20;
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            double thresh = 2.0 + 4.0*i/(double)nBins;
            double thick = 0 + 20*j/(double)nBins;
            std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(thick, thresh, events, randU01s, randDeathTimes);
            double chisqWgt = calcChisqWgt(refHist, hist1);
            double chisqUnWgt = calcChisq(refHist, hist1);
            double chisqNate = calcChisqNate(refHist, hist1);
//                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
            printf("%f %f %f %f %f\n", thick, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
            fflush(stdout);
        }
    }
    
    /*int nBins = 100;
    for(int i = 0; i < nBins+1; i++) {
        double contam = 0.2 + 0.2*i/(double)nBins;
        std::vector<weightedBin> hist1 = createHistQuantEdE(contam, 0.0, events, randU01s, randDeathTimes);
        double chisqWgt = calcChisqWgt(refHist, hist1);
        double chisqUnWgt = calcChisq(refHist, hist1);
        double chisqNate = calcChisqNate(refHist, hist1);
//                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
        printf("%f %f %f %f %f\n", contam, 0.0, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
        fflush(stdout);
    }*/

//    std::vector<weightedBin> hist1 = createHist(0.155, 0.0, 15.377918, 17.428307, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHist(0.26, 2.1, 5.125973, 30.755837, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHist(0.225, 1.1, 7.176362, 22.9, events, randU01s, randDeathTimes); //Nate's chisq
//    std::vector<weightedBin> hist1 = createHist(0.27, 2.3, 5.125973, 34.5, events, randU01s, randDeathTimes); //Nate's chisq EdE
//      std::vector<weightedBin> hist1 = createHist(0.2551, 1.892, 5.701, 29.19, events, randU01s, randDeathTimes); //Prof II
//      std::vector<weightedBin> hist1 = createHist(0.285, 2.6, 4.100778, 34.5, events, randU01s, randDeathTimes); //Prof II
//      std::vector<weightedBin> hist1 = createHistQuant(0.1, 4.100778, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.4, 4.100778, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.444444, 3.8, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.363636, 5.8, 34.5, events, randU01s, randDeathTimes); //Test of quantum frac
      //std::vector<weightedBin> hist1 = createHistQuant(0.35, 4.4, 34.5, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuant(0.31, 6.0, 24.0, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuantEdE(0.35, 0.0, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuantEdE(0.32, 0.0, events, randU01s, randDeathTimes); //quant EdE opt, no cutoff
/*      double chisq = calcChisqWgt(refHist, hist1);
      double chisq2 = calcChisq(refHist, hist1);
      printf("%f %f\n\n", chisq/hist1.size(), chisq2/hist1.size());
      for(auto it = hist1.begin(); it < hist1.end(); it++) {
          printf("%f,", it->wgt);
      }
      printf("\n");*/

/*    char paramFile[256];
    char histFile[256];
    double eff;
    double thr;
    double cut;
    double sat;
    for(int i = 0; i < 1296; i++) {
        sprintf(paramFile, "professor/scan/%04d/params.dat", i);
        sprintf(histFile, "professor/scan/%04d/hist.yoda", i);
        std::ifstream params(paramFile, std::ios::in);
        std::string input;
        std::string paramVal;

        std::getline(params, input);
        paramVal = input.substr(input.find_first_of("-0123456789"));
//        printf("%s\n", paramVal.c_str());
        eff = atof(paramVal.c_str());

        std::getline(params, input);
        paramVal = input.substr(input.find_first_of("-0123456789"));
//        printf("%s\n", paramVal.c_str());
        thr = atof(paramVal.c_str());

        std::getline(params, input);
        paramVal = input.substr(input.find_first_of("-0123456789"));
//        printf("%s\n", paramVal.c_str());
        cut = atof(paramVal.c_str());

        std::getline(params, input);
        paramVal = input.substr(input.find_first_of("-0123456789"));
//        printf("%s\n", paramVal.c_str());
        sat = atof(paramVal.c_str());

        params.close();
        printf("%f %f %f %f\n", eff, thr, cut, sat);
        std::vector<weightedBin> hist1 = createHist(eff, thr, cut, sat, events, randU01s, randDeathTimes);
        writeYOD(hist1, histFile);
    }
    writeYODDouble(refHist, "professor/hist.yoda");*/

    delete[] randU01s;
    delete[] buf;

    return 0;
}
