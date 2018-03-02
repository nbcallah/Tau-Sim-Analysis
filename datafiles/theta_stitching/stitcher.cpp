#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>

#define NRECORDS 50

typedef struct evt_in {
    double energy;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt_in;

typedef struct evt_out {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt_out;

typedef struct evt_th {
    double energy;
    double theta;
} evt_th;


std::vector<evt_in> readDataFile(const char* fName) {
    std::vector<evt_in> events;
    
    const size_t buff_len = sizeof(unsigned int) + 1*sizeof(double) + 2*NRECORDS*sizeof(float) + sizeof(unsigned int);
    char* buf = new char[buff_len];
    std::ifstream binfile(fName, std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! Could not open file %s\n", fName);
        return events;
    }
    
    evt_in event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        if(*((unsigned int *)&buf[0]) != 1*8 + 2*NRECORDS*4) {
            fprintf(stderr, "Error! Aliased read on %s\n", fName);
            exit(2);
        }
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double)), NRECORDS*sizeof(float));
        std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events from %s!\n", events.size(), fName);

    delete[] buf;
    return events;
}

std::vector<evt_th> readThetaFile(const char* fName) {
    std::vector<evt_th> events;
    
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + sizeof(unsigned int);
    char* buf = new char[buff_len];
    std::ifstream binfile(fName, std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! Could not open file %s\n", fName);
        return events;
    }
    
    evt_th event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        if(*((unsigned int *)&buf[0]) != 2*8) {
            fprintf(stderr, "Error! Aliased read on %s\n", fName);
            exit(2);
        }
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events from %s!\n", events.size(), fName);

    delete[] buf;
    return events;
}

void writeData(const char* fName, std::vector<evt_out> events) {
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float) + sizeof(unsigned int);
    char* buf = new char[buff_len];
    
    std::ofstream binfile(fName, std::ios::out | std::ios::binary);
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! Could not open file %s\n", fName);
        return;
    }
    for(auto it = events.begin(); it < events.end(); it++) {
        *((unsigned int *)(&buf[0])) = 2*8 + 2*NRECORDS*4;
        *((double *)(&buf[0] + sizeof(unsigned int))) = it->energy;
        *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double))) = it->theta;
        std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(it->times[0]), NRECORDS*sizeof(float));
        std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(it->ePerp[0]), NRECORDS*sizeof(float));
        *((unsigned int *)(&buf[0] + sizeof(unsigned int) + sizeof(double) + 2*NRECORDS*sizeof(float))) = 2*8 + 2*NRECORDS*4;
        binfile.write(buf, buff_len);
    }
    
    printf("Wrote %lu Events to %s!\n", events.size(), fName);
    
    binfile.close();
    
    delete[] buf;
}

int main(int argc, char** argv) {
    if(argc != 4) {
        printf("Error! Usage: ./minimal data_in theta data_out\n");
        return 1;
    }

    std::vector<evt_in> data_in = readDataFile(argv[1]);
    std::vector<evt_th> theta = readThetaFile(argv[2]);
    std::vector<evt_out> data_out;
    
    assert(data_in.size() == theta.size());
    
    data_out.resize(data_in.size());
    
    for(unsigned long i = 0; i < data_in.size(); i++) {
        assert(fabs((data_in[i].energy - theta[i].energy)/theta[i].energy) < 5e-5);
        data_out[i].energy = data_in[i].energy;
        data_out[i].theta = theta[i].theta;
        std::memcpy((void *)&(data_out[i].times[0]), (void *)&(data_in[i].times[0]), NRECORDS*sizeof(float));
        std::memcpy((void *)&(data_out[i].ePerp[0]), (void *)&(data_in[i].ePerp[0]), NRECORDS*sizeof(float));
    }
    
    writeData(argv[3], data_out);
    
    return 0;
}
