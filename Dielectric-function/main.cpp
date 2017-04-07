#include <iostream>
#include <seels.h>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
    string basisFile = "BASISFILE";
    string waveFile = "WAVEFILE";

    SEELS Simulate(0,10,0.1);
    Simulate.readBASISFILE(basisFile);
    Simulate.readWAVEFILE(waveFile);

    Simulate.calculateSpectrum(vec3(0,0,0));

    return 0;
}
