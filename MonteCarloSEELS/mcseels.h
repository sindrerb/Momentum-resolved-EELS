#ifndef MCSEELS_H
#define MCSEELS_H

#include <random>
#include <string>
#include <../vec3-class/vec3.h>
#include <../wavestate-class/wavestate.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>

class MCSEELS
{
public:
    MCSEELS();

    bool readBasisfile(std::string FILENAME);
    bool readWavefile(std::string FILENAME);

    void defineMomentumRange(float momentumMin, float momentumMax, float momentumBin);
    void defineEnergyRange(float energyMin, float energyMax, float energyBin);

    float fermiDirac(float energy, float energyFermi, float temperature);

    void calculateSpectrum(int cycles);
    void writeSpectrum(std::string FILENAME);

    float momentumMin() const;
    void setMomentumMin(float momentumMin);

    float momentumMax() const;
    void setMomentumMax(float momentumMax);

    float momentumBin() const;
    void setMomentumBin(float momentumBin);

    float energyMin() const;
    void setEnergyMin(float energyMin);

    float energyMax() const;
    void setEnergyMax(float energyMax);

    float energyBin() const;
    void setEnergyBin(float energyBin);

    bool setSpectrum();
    void setSpectrum(float momentumMin, float momentumMax, float momentumBin, float energyMin, float energyMax, float energyBin);
    float **spectrum() const;

    int momentumLength() const;
    void setMomentumLength(int momentumLength);

    int energyLength() const;
    void setEnergyLength(int energyLength);

    int waveBasisLength() const;
    void setWaveBasisLength(int waveBasisLength);

    int wavestatesLength() const;
    void setWavestatesLength(int wavestatesLength);

    float fermiEnergy() const;
    void setFermiEnergy(float fermiEnergy);

    float temperature() const;
    void setTemperature(float temperature);

private:
    //System spesifications
    float m_fermiEnergy;
    float m_temperature;

    //Momentum dispersion
    float m_momentumMin, m_momentumMax, m_momentumBin;

    //Energy dispersion
    float m_energyMin, m_energyMax, m_energyBin;

    //Spectrum matrix
    float **m_spectrum;
    int m_momentumLength;
    int m_energyLength;

    //WaveBasis
    std::vector<vec3> m_waveBasis;
    int m_waveBasisLength;

    //WaveStates
    std::vector<waveState> m_wavestates;
    int m_wavestatesLength;


};

#endif // MCSEELS_H
