#ifndef SEELS_H
#define SEELS_H

#include<vector>
#include<../wavestate-class/wavestate.h>
#include<../vec3-class/vec3.h>

class SEELS
{
public:
    SEELS();
    SEELS(double energyMinimum,double energyMaximum,double dispersion);

    void setEnergyRange(double energyMinimum,double energyMaximum,double dispersion);
    void writeEnergyRangeToFile();





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
    double *spectra;

    //Calculation parameters
    vec3 momentumTransfer;




};

#endif // SEELS_H
