#include "seels.h"

SEELS::SEELS() {

}

SEELS::SEELS(double energyMinimum, double energyMaximum, double dispersion) {
    setEnergyMin(energyMinimum);
    setEnergyMax(energyMaximum);
    setEnergyDispersion(dispersion);
}

double SEELS::getEnergyMin() const
{
    return energyMin;
}

void SEELS::setEnergyMin(double value)
{
    energyMin = value;
}

double SEELS::getEnergyMax() const
{
    return energyMax;
}

void SEELS::setEnergyMax(double value)
{
    energyMax = value;
}

double SEELS::getEnergyDispersion() const
{
    return energyDispersion;
}

void SEELS::setEnergyDispersion(double value)
{
    energyDispersion = value;
}


