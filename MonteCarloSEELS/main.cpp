#include <iostream>
#include <mcseels.h>
using namespace std;

int main(int argc, char *argv[]) {
    float energyMin, energyMax, energyBin;
    float momentumMin, momentumMax, momentumBin;

    //Define energy spectrum
    energyMin = 0;
    energyMax = 10;
    energyBin = 0.2;

    //Define momentum dispersion
    momentumMin = 0;
    momentumMax = 0.5;
    momentumBin = 0.1;

    MCSEELS MC;

    bool basis = MC.readBasisfile("BASISFILE");
    bool waves = MC.readWavefile("WAVEFILE");

    if(!basis || !waves){
        cout << "An error occured during loading files." << endl;
        return -1;
    }

    MC.setFermiEnergy(10.5);
    MC.setTemperature(1);

    MC.defineEnergyRange(energyMin, energyMax, energyBin);
    MC.defineMomentumRange(momentumMin, momentumMax, momentumBin);

    bool spectrum = MC.setSpectrum();

    if(!spectrum) {
        cout << "An error occured during creating empty spectrum matrix" << endl;
        return -1;
    }

    MC.calculateSpectrum(1);
    MC.writeSpectrum("SPECTRA");



    return 0;
}
