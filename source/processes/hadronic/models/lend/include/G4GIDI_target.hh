//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef G4GIDI_target_h_included
#define G4GIDI_target_h_included 1

#include <vector>
#include <string>

//using namespace std;

#include <statusMessageReporting.h>

#include <MCGIDI.h>

typedef struct crossSectionData_s crossSectionData;
typedef struct G4GIDI_Product_s G4GIDI_Product;

struct crossSectionData_s {
    int start, end;
    std::vector<double> crossSection;
};

#define channelID std::string

struct G4GIDI_Product_s {
    int A, Z, m;
    double kineticEnergy, px, py, pz;
    double birthTimeSec;
};

class G4GIDI_target {

    public:
        void init( const char *fileName );
        std::string equalProbableBinSampleMethod;
        int nElasticIndices, nCaptureIndices, nFissionIndices, nOthersIndices;
        int *elasticIndices, *captureIndices, *fissionIndices, *othersIndices;

    public:
        GIDI::statusMessageReporting smr;
        int projectilesPOPID;
        std::string name;
        std::string sourceFilename;
        double mass;
        GIDI::MCGIDI_target *target;

        G4GIDI_target( const char *fileName );
        G4GIDI_target( std::string const &fileName );       
        ~G4GIDI_target( );

        std::string *getName( void );
        std::string *getFilename( void );
        int getZ( void );
        int getA( void );
        int getM( void );
        double getMass( void );
        int getTemperatures( double *temperatures );
        int readTemperature( int index );
        std::string getEqualProbableBinSampleMethod( void );
        int setEqualProbableBinSampleMethod( std::string method );

        int getNumberOfChannels( void );
        int getNumberOfProductionChannels( void );
        channelID getChannelsID( int channelIndex );
        std::vector<channelID> *getChannelIDs( void );
        std::vector<channelID> *getProductionChannelIDs( void );

        std::vector<double> *getEnergyGridAtTIndex( int index );

        double getTotalCrossSectionAtE( double e_in, double temperature );
        double getElasticCrossSectionAtE( double e_in, double temperature );
        double getCaptureCrossSectionAtE( double e_in, double temperature );
        double getFissionCrossSectionAtE( double e_in, double temperature );
        double getOthersCrossSectionAtE( double e_in, double temperature );
        double sumChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature );
        int sampleChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );

        double getElasticFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getCaptureFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getFissionFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getOthersFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getFinalState( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );

        double getReactionsThreshold( int index );
        double getReactionsDomain( int index, double *EMin, double *EMax );
};

#endif      // End of G4GIDI_target_h_included
