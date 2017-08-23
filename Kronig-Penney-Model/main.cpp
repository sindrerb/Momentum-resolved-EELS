#include <iostream>
#include "kronigpenney.h"
#include "../wavestate-class/wavestate.h"
#include "../vec3-class/vec3.h"
#include <math.h>
#include <fstream>

#include "mpi.h"
#include <stdio.h>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{


// MPI VARIABLES
    double startUniversal,stopUniversal,step, totalSteps;

    int myStartInt, myStopInt;
    double myIntervalStart,myIntervalEnd;
    string myFileEnd;

    startUniversal = -0.25;
    stopUniversal = 0.25;
    step = 0.01;
    totalSteps = (stopUniversal-startUniversal)/step;


    // PROGRAM VARIABLES
        double energyCutOff = 40; //eV
        string cellFile = "CELLFILE";
        string potentialFile = "POTENTIALFILE";
        string basisFile = "BASISFILE";
        string waveFile = "WAVEFILE";
        string oldFile = waveFile+"_OLD";
        string resultFile = "RESULTFILE";


//INITIALIZE CELL AND PARAMETERS
    KronigPenney KP;

    //int result;

    KP.initializeCELL(cellFile);

    bool basis = KP.readBASISFILE(basisFile);
    if(!basis){
        KP.setWaveBasis(energyCutOff);
        KP.writeBASISFILE(basisFile);
    }

    string dummystring;
    fstream POTENTIAL(potentialFile, std::ios_base::in);

    //"Potential [eV], Periodicity (a,b,c)\n\n-10\t1.0000\t1.0000\t1.0000\n"
    POTENTIAL >> dummystring >> dummystring >> dummystring >> dummystring;


    double potential, x,y,z, volume;
    vec3 a,b,c;
    vec3 kPoint;

// Initialize MPI
    int numprocs, myRank;
    char hostname[256];

    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    gethostname(hostname,255);
    if(myRank ==  0 && argc <= 1) {
        cout << "Total: " << totalSteps << endl;
    }

// Controlling that every slave contributes
    printf("Process number %d of %d on host %s is now active.\n", myRank, numprocs, hostname);

//Evaluate where my interval begins and ends
    myStartInt = totalSteps*myRank/numprocs;
    myStopInt =  (totalSteps*(myRank+1))/numprocs-1;

    myIntervalStart =startUniversal+myStartInt*step;
    myIntervalEnd = startUniversal+myStopInt*step;


    while(POTENTIAL){ //LOOP OVER THE DEFINED POTENTIALS
        POTENTIAL >> potential >> x >> y >> z;
        if(potential == 0){
            MPI_Barrier(MPI_COMM_WORLD);
            if(myRank == 0){
                cout << "Combining!" << endl;
                ofstream combinedWave(waveFile, std::ios_base::app);
                ofstream combinedResult(resultFile, std::ios_base::app);

                for(int i = 0; i<numprocs; i++){
                    ifstream inWave(waveFile+to_string(i));
                    combinedWave << inWave.rdbuf();
                    inWave.close();
                    remove( (waveFile+to_string(i)).c_str() );

                    ifstream inResult(resultFile+to_string(i));
                    combinedResult << inResult.rdbuf();
                    inResult.close();
                    remove( (resultFile+to_string(i)).c_str() );
                }

                combinedResult.close();
                combinedWave.close();

//                rename("WAVEFILE","WAVEFILE_OLD");
//                rename("RESULTFILE","WAVEFILE");
//                rename("WAVEFILE_OLD","WAVEFILE_ORIGINAL");
            }
            POTENTIAL.close();
            MPI_Finalize();
            return 0;
        }

        if(myRank == 0) {
            cout << "Calculating with a potential of " << potential << " eV." << endl;
        }
        a = x*KP.aReal();
        b = y*KP.bReal();
        c = z*KP.cReal();
        volume = a.dot(b.cross(c));

        //loop over desired k-points
        for(double h = myIntervalStart; h<=myIntervalEnd; h+= step){
            for(double k = startUniversal; k<=stopUniversal; k+= step){
                //for(double l = startUniversal; l<=stopUniversal; l+= step){
                    kPoint = h*KP.aResiprocal()+k*KP.bResiprocal();//+l*KP.cResiprocal();

                    bool states = KP.readWAVEFILE(waveFile+to_string(myRank),kPoint);
                    if(!states){
                        KP.setWaveStates(kPoint);
                        KP.writeWAVEFILE(waveFile+to_string(myRank),kPoint);
                    }
                    KP.calculateEigenValues(kPoint,5,20.0, potential, volume);
                    KP.writeRESULTFILE(resultFile+to_string(myRank));
                //}
            }
        }
    }
    POTENTIAL.close();
    MPI_Finalize();

    return 0;
}
