#ifndef SEELS_H
#define SEELS_H

#include<vector>
#include<string>
#include<../wavestate-class/wavestate.h>
#include<../vec3-class/vec3.h>

class SEELS
{
public:
    SEELS();
    SEELS(double energyMinimum,double energyMaximum,double dispersion);

    void setEnergyRange(double energyMinimum,double energyMaximum,double dispersion);
    void writeEnergyRangeToFile();

    void readWavestates(std::string WAVEFILE);

    void calculateSpectrum(vec3 q);
    void addSpectrumToFile();


    double getEnergyMin() const;
    void setEnergyMin(double value);

    double getEnergyMax() const;
    void setEnergyMax(double value);

    double getEnergyDispersion() const;
    void setEnergyDispersion(double value);

private:
    //EELS spectra details
    double energyMin;
    double energyMax;
    double energyDispersion;

    //Spectra storage
    double *energyRange;
    double *spectrum;

    //Calculation parameters
    vec3 momentumTransfer;




};

#endif // SEELS_H
