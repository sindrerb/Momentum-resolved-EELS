#include "seels.h"

SEELS::SEELS() {

}

SEELS::SEELS(double energyMinimum, double energyMaximum, double dispersion) {
    setEnergyMin(energyMinimum);
    setEnergyMax(energyMaximum);
    setEnergyDispersion(dispersion);


    m_spectrumLength = (int) (m_energyMax-m_energyMin)/(m_energyDispersion);
    m_spectrum = new double[m_spectrumLength];
    m_energyRange = new double[m_spectrumLength];
    for(int i = 0; i<m_spectrumLength; i++){
        m_spectrum[i] = 0;
        m_energyRange[i] = m_energyMin+i*m_energyDispersion;
    }
}

void SEELS::writeEnergyRangeToFile(std::__cxx11::string ENERGYFILE) {
    std::ofstream ENERGY(ENERGYFILE);

    for(int i = 0; i<m_spectrumLength; i++) {
        ENERGY << m_energyRange[i] << "\t";
    }
    ENERGY.close();
}

bool SEELS::readBASISFILE(std::__cxx11::string BASISFILE) {
    m_waveBasis.clear();
    m_waveBasisLength = 0;
    std::fstream BASIS(BASISFILE, std::ios_base::in);

    if(!BASIS.good()){
        std::cout << "Could not find a BASISFILE" << std::endl;
        return false;
    }else{
        double x,y,z;
        while(BASIS >> x >> y >> z){
            m_waveBasis.push_back(vec3(x,y,z));
            m_waveBasisLength ++;
        }
        BASIS.close();
        return true;
    }
}


bool SEELS::readWAVEFILE(std::__cxx11::string WAVEFILE) {
    m_wavestates.clear();
    m_wavestatesLength = 0;

    std::fstream WAVES(WAVEFILE, std::ios_base::in);

    double energy;
    std::vector<double> weights;

    if(!WAVES.good()){ //If kPoint is not contained in file
        WAVES.close();
        std::cout << "Could not find a WAVEFILE" << std::endl;
        return false;
    }else{
        double x,y,z;
        double weight;
        vec3 kPointFromFile, effectiveG;
        while(WAVES) {
            WAVES >> x >> y >> z;
            WAVES >> energy;
            kPointFromFile = vec3(x,y,z);

            weights.clear();
            for(int i = 0; i<m_waveBasisLength; i++) {
                WAVES >> weight;
                weights.push_back(weight);
            }
            for(int j = 0; j<m_waveBasisLength; j++) {
                effectiveG += m_waveBasis[j]*weights[j];
            }
            m_wavestates.push_back(waveState(kPointFromFile,effectiveG,energy,weights));
            m_wavestatesLength ++;
        }

        //m_unperturbedStates.pop_back();
        if(m_wavestatesLength > 0){
            m_wavestatesLength --;
        }

        if(m_wavestatesLength == 0) {
            std::cout << "The file does not contain any reciprocal points" << std::endl;
            return false;
        }else{
            return true;
        }
    }
}

void SEELS::calculateSpectrum(vec3 q)
{
    m_matrixElements.clear();
    m_matrixElementsLenght = 0;

    m_energyTransitions.clear();
    m_energyTransitionsLength = 0;

    m_fermiWeights.clear();
    m_fermiWeightsLength = 0;

    vec3 momentumTransfer;
    double matrixElement, fermiWeight;
    for(int initial=0; initial<m_wavestatesLength; initial++) {
        for(int final=initial; final<m_wavestatesLength; final++) {

            momentumTransfer = m_wavestates[final].getK() -m_wavestates[initial].getK();
            if(m_wavestates[final].energy() <= 25.0 && m_wavestates[initial].energy() >= 6.0) {
                if( (momentumTransfer-q).length() < 0.001) {
                    matrixElement = 0;
                    for(int i = 0; i<m_waveBasisLength; i++) {
                        matrixElement += m_wavestates[final].weight(i)*m_wavestates[initial].weight(i)*(m_waveBasis[i].length()+m_wavestates[initial].getK()+0.5*q).length();
                    }
                    m_matrixElements.push_back(matrixElement*matrixElement);
                    m_matrixElementsLenght ++;

                    m_energyTransitions.push_back(m_wavestates[final].energy()-m_wavestates[initial].energy());
                    m_energyTransitionsLength ++;

                    fermiWeight = FermiDirac(m_wavestates[final].energy())-FermiDirac(m_wavestates[initial].energy());
                    m_fermiWeights.push_back(fermiWeight);
                    m_fermiWeightsLength ++;
    //                std::cout << matrixElement << "\t" << m_energyTransitions[m_energyTransitionsLength-1] << std::endl;
                }
            }
        }
    }
    std::cout << m_matrixElementsLenght << " # " << m_energyTransitionsLength << " # " << m_fermiWeightsLength << std::endl;

    std::complex<double> epsilon, spectrumIntensity;
    epsilon = 0.001*1i;
    double energy;
    energy = m_energyMin;
    for(int i = 0; i<m_spectrumLength; i++) { //
        spectrumIntensity = 0;
        for(int j = 0; j<m_matrixElementsLenght; j++) {
            spectrumIntensity += m_matrixElements[j]*m_fermiWeights[j]/(energy-m_energyTransitions[j]+epsilon);
        }
        m_spectrum[i] = spectrumIntensity.real();
        energy += m_energyDispersion;
    }
}

void SEELS::addSpectrumToFile(std::__cxx11::string SPECTRUMFILE) {
    std::ofstream SPECTRA(SPECTRUMFILE, std::ios_base::app);

    for(int i = 0; i<m_spectrumLength; i++) {
        SPECTRA << m_spectrum[i] << "\t";
    }
    SPECTRA << "\n";
    SPECTRA.close();
}

void SEELS::addKPointToFile(std::__cxx11::string KPOINTFILE, vec3 q) {
    std::ofstream KPOINT(KPOINTFILE, std::ios_base::app);
    KPOINT << q.x() << "\t" << q.y() << "\t" << q.z() << "\n";
    KPOINT.close();
}


double SEELS::FermiDirac(double energy) {
    return 1.0/(exp((energy-m_fermilevel)/(m_temperature))+1.0);
}
double SEELS::getEnergyMin() const
{
    return m_energyMin;
}

void SEELS::setEnergyMin(double value)
{
    m_energyMin = value;
}

double SEELS::getEnergyMax() const
{
    return m_energyMax;
}

void SEELS::setEnergyMax(double value)
{
    m_energyMax = value;
}

double SEELS::getEnergyDispersion() const
{
    return m_energyDispersion;
}

void SEELS::setEnergyDispersion(double value)
{
    m_energyDispersion = value;
}

double SEELS::getFermilevel() const
{
    return m_fermilevel;
}

void SEELS::setFermilevel(double value)
{
    m_fermilevel = value;
}

vec3 SEELS::getMomentumTransfer() const
{
    return m_momentumTransfer;
}

void SEELS::setMomentumTransfer(const vec3 &value)
{
    m_momentumTransfer = value;
}

double SEELS::getTemperature() const
{
    return m_temperature;
}

void SEELS::setTemperature(double temperature)
{
    m_temperature = temperature;
}


