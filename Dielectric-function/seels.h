#ifndef SEELS_H
#define SEELS_H

#include <string>
#include <../vec3-class/vec3.h>
#include <../wavestate-class/wavestate.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>

class SEELS
{
public:
    SEELS();
    SEELS(double energyMinimum,double energyMaximum,double dispersion);

    void setEnergyRange(double energyMinimum,double energyMaximum,double dispersion);
    void writeEnergyRangeToFile(std::string ENERGYFILE);

    bool readBASISFILE(std::string BASISFILE);
    bool readWAVEFILE(std::string WAVEFILE);

    void calculateSpectrum(vec3 q);
    void addSpectrumToFile(std::string SPECTRUMFILE);

    void addKPointToFile(std::string KPOINTFILE, vec3 q);

    double FermiDirac(double energy);

    double getEnergyMin() const;
    void setEnergyMin(double value);

    double getEnergyMax() const;
    void setEnergyMax(double value);

    double getEnergyDispersion() const;
    void setEnergyDispersion(double value);

    double getFermilevel() const;
    void setFermilevel(double value);

    vec3 getMomentumTransfer() const;
    void setMomentumTransfer(const vec3 &value);

    double getTemperature() const;
    void setTemperature(double temperature);

private:
    //EELS spectra details
    double m_energyMin;
    double m_energyMax;
    double m_energyDispersion;


    //WaveBasis
    int m_waveBasisLength;
    std::vector<vec3> m_waveBasis;

    //WaveStates
    std::vector<waveState> m_wavestates;
    int m_wavestatesLength;


    //Dielectric function calculation
    std::vector<double> m_matrixElements;
    int m_matrixElementsLenght;

    std::vector<double> m_energyTransitions;
    int m_energyTransitionsLength;

    std::vector<double> m_fermiWeights;
    int m_fermiWeightsLength;

    //Calculation parameters
    vec3 m_momentumTransfer;
    double m_temperature;
    double m_fermilevel;

    //Spectra storage
    double *m_energyRange;
    double *m_spectrum;
    int m_spectrumLength;

    
};

#endif // SEELS_H
