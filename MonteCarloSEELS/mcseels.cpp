#include "mcseels.h"

MCSEELS::MCSEELS() {
    setEnergyMin(0.0);
    setEnergyMax(1.0);
    setEnergyBin(0.1);
    setMomentumMin(0.0);
    setMomentumMax(1.0);
    setMomentumBin(0.1);

    setMomentumLength(10);
    setEnergyLength(10);
}

bool MCSEELS::readBasisfile(std::__cxx11::string FILENAME) {
    m_waveBasis.clear();
    m_waveBasisLength = 0;
    std::fstream BASIS(FILENAME, std::ios_base::in);

    if(!BASIS.good()){
        std::cout << "Could not find a basisfile, check location and filename." << std::endl;
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

bool MCSEELS::readWavefile(std::__cxx11::string FILENAME) {
    m_wavestates.clear();
    m_wavestatesLength = 0;

    std::fstream WAVES(FILENAME, std::ios_base::in);

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
        WAVES.close();
        if(m_wavestatesLength == 0) {
            std::cout << "The file does not contain any reciprocal points" << std::endl;
            return false;
        }else{
            return true;
        }
    }
}

void MCSEELS::defineMomentumRange(float momentumMin, float momentumMax, float momentumBin) {
    int length;
    setMomentumMin(momentumMin);
    setMomentumMax(momentumMax);
    setMomentumBin(momentumBin);

    length = (int) ((momentumMax-momentumMin)/momentumBin);
    setMomentumLength(length);
}

void MCSEELS::defineEnergyRange(float energyMin, float energyMax, float energyBin) {
    int length;
    setEnergyMin(energyMin);
    setEnergyMax(energyMax);
    setEnergyBin(energyBin);

    length = (int) ((energyMax-energyMin)/energyBin);
    setEnergyLength(length);
}

float MCSEELS::fermiDirac(float energy, float energyFermi, float temperature) {
    float FD;

    FD = 1.0/(exp((energy-energyFermi)/(temperature*0.0259))+1);

    return FD;
}

void MCSEELS::calculateSpectrum(int cycles) {

    // Initialize the seed and call the Mersienne algo
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<float> randomGenerator(0.0,1.0);

    int initial, final;
    float energyInitial, energyFinal, energyTransfer, fermi, matrixElement;
    vec3 momentumInitial, momentumFinal, momentumTransfer;

    float metropolis;
    int momentumPlace, energyPlace;
    for(int cycle = 0; cycle<cycles; cycle++) {
        for(int count = 0; count < m_wavestatesLength; count ++){
            initial = (int) (randomGenerator(gen)*m_wavestatesLength) +1;
            final = (int) (randomGenerator(gen)*m_wavestatesLength) +1;
            energyInitial = m_wavestates[initial].energy();
            energyFinal = m_wavestates[final].energy();
            energyTransfer = energyFinal-energyInitial;
            momentumInitial = m_wavestates[initial].getK();
            momentumFinal = m_wavestates[final].getK();
            momentumTransfer = momentumFinal-momentumInitial;

            fermi = fermiDirac(energyInitial,m_fermiEnergy,m_temperature)-fermiDirac(energyFinal,m_fermiEnergy,m_temperature);
            matrixElement = 0;

            metropolis = randomGenerator(gen);
            if(fermi > metropolis) {
                for(int i = 0; i<m_waveBasisLength; i++) {
                    matrixElement += (m_wavestates[final].weight(i)*m_wavestates[initial].weight(i)*(m_waveBasis[i]+m_wavestates[initial].getK()+0.5*momentumTransfer).length());
                }
                momentumPlace = (int) (momentumTransfer.length()/m_momentumBin);
                energyPlace = (int) (energyTransfer/m_energyBin);

                if(momentumPlace < m_momentumLength && energyPlace < m_energyLength) {
                    m_spectrum[momentumPlace][energyPlace] += matrixElement;
                }
                std::cout << matrixElement << std::endl;

            }

        }

    }
}

void MCSEELS::writeSpectrum(std::__cxx11::string FILENAME) {
    std::ofstream FILE(FILENAME);

    for(int momentumPlace = 0; momentumPlace<m_momentumLength; momentumPlace++) {
        for(int energyPlace = 0; energyPlace<m_energyLength; energyPlace++) {
            FILE << m_spectrum[momentumPlace][energyPlace] << "\t";
        }
        FILE << "\n";
    }
    FILE.close();
}

bool MCSEELS::setSpectrum()
{
    if( m_momentumLength==0 || m_energyLength==0 ) {
        return false;
    }

    m_spectrum = new float*[m_momentumLength];
    for(int i = 0; i<m_momentumLength; i++ ) {
        m_spectrum[i] = new float[m_energyLength];
        for(int j = 0; j<m_momentumLength; j++ ) {
            m_spectrum[i][j] = 0.0;
        }
    }

    return true;
}

void MCSEELS::setSpectrum(float momentumMin, float momentumMax, float momentumBin, float energyMin, float energyMax, float energyBin) {
    defineEnergyRange(energyMin, energyMax, energyBin);
    defineMomentumRange(momentumMin, momentumMax, momentumBin);

    setSpectrum();
}

float **MCSEELS::spectrum() const
{
    return m_spectrum;
}


float MCSEELS::momentumMin() const
{
    return m_momentumMin;
}

void MCSEELS::setMomentumMin(float momentumMin)
{
    m_momentumMin = momentumMin;
}

float MCSEELS::momentumMax() const
{
    return m_momentumMax;
}

void MCSEELS::setMomentumMax(float momentumMax)
{
    m_momentumMax = momentumMax;
}

float MCSEELS::momentumBin() const
{
    return m_momentumBin;
}

void MCSEELS::setMomentumBin(float momentumBin)
{
    m_momentumBin = momentumBin;
}

float MCSEELS::energyMin() const
{
    return m_energyMin;
}

void MCSEELS::setEnergyMin(float energyMin)
{
    m_energyMin = energyMin;
}

float MCSEELS::energyMax() const
{
    return m_energyMax;
}

void MCSEELS::setEnergyMax(float energyMax)
{
    m_energyMax = energyMax;
}

float MCSEELS::energyBin() const
{
    return m_energyBin;
}

void MCSEELS::setEnergyBin(float energyBin)
{
    m_energyBin = energyBin;
}

int MCSEELS::momentumLength() const
{
    return m_momentumLength;
}

void MCSEELS::setMomentumLength(int momentumLength)
{
    m_momentumLength = momentumLength;
}

int MCSEELS::energyLength() const
{
    return m_energyLength;
}

void MCSEELS::setEnergyLength(int energyLength)
{
    m_energyLength = energyLength;
}

int MCSEELS::waveBasisLength() const
{
    return m_waveBasisLength;
}

void MCSEELS::setWaveBasisLength(int waveBasisLength)
{
    m_waveBasisLength = waveBasisLength;
}

int MCSEELS::wavestatesLength() const
{
    return m_wavestatesLength;
}

void MCSEELS::setWavestatesLength(int wavestatesLength)
{
    m_wavestatesLength = wavestatesLength;
}

float MCSEELS::fermiEnergy() const
{
    return m_fermiEnergy;
}

void MCSEELS::setFermiEnergy(float fermiEnergy)
{
    m_fermiEnergy = fermiEnergy;
}

float MCSEELS::temperature() const
{
    return m_temperature;
}

void MCSEELS::setTemperature(float temperature)
{
    m_temperature = temperature;
}
