#include <iostream>
#include <seels.h>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
    string basisFile = "BASISFILE";
    string waveFile = "WAVEFILE";

    SEELS Simulate(0,10,0.01);
    Simulate.readBASISFILE(basisFile);
    Simulate.readWAVEFILE(waveFile);

    Simulate.setFermilevel(15);
    Simulate.setTemperature(0.1);

    vec3 k,q;
    k = vec3(0.00125664,0,0);
    Simulate.writeEnergyRangeToFile("ENERGYFILE");

    for(int i = 0; i<10; i++){
        q = i*k;
        Simulate.calculateSpectrum(q);
        Simulate.addKPointToFile("KPOINTFILE", q);
        Simulate.addSpectrumToFile("SPECTRUMFILE");
    }

    return 0;
}
