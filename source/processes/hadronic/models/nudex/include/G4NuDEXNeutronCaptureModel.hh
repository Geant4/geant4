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
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NuDEXNeutronCaptureModel
//
//      Author:        E.Mendoza & A.Ribon
// 
//      Creation date: 29 May 2024
//
//      Description:   This class (a proxy of the class G4NuDEX) uses
//                     the NuDEX model to produce gammas and internal
//                     conversion electrons from neutron capture.
//                     Whenever NuDEX is not applicable, G4PhotonEvaporation
//                     is used.
//                     The implementation of this class follows the code
//                     of the class G4NeutronRadCapture.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
// Class to use NuDEX model inside Geant4
// 

#ifndef G4NUDEXNEUTRONCAPTUREMODEL_HH
#define G4NUDEXNEUTRONCAPTUREMODEL_HH 1

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

class G4NuDEXStatisticalNucleus;
class G4VEvaporationChannel;


#define G4NUDEX_MAXZA 120000


class G4NuDEXNeutronCaptureModel : public G4HadronicInteraction {
  public:
    explicit G4NuDEXNeutronCaptureModel();
    virtual ~G4NuDEXNeutronCaptureModel();
  
    virtual G4HadFinalState* ApplyYourself( const G4HadProjectile &aTrack, G4Nucleus &targetNucleus ) final;
    virtual void InitialiseModel() final;

  private:
    G4NuDEXNeutronCaptureModel & operator=( const G4NuDEXNeutronCaptureModel &right ) = delete;
    G4NuDEXNeutronCaptureModel( const G4NuDEXNeutronCaptureModel& ) = delete;

    G4int GenerateNeutronCaptureCascade( G4int theZ, G4int theA, G4double NeutronEnergy, G4int InitialLevel,
  				         std::vector< char >& pType, std::vector< G4double >& pEnergy, std::vector< G4double >& pTime );

    // Initial level for neutron capture. If jspinx2v < 0 it is sampled according to the 2J+1 rule
    // l-spin = 0, 1, 2 --> s-wave, p-wave, d-wave ...
    G4int SelectInitialLevel( G4int theCompoundZ, G4int theCompoundA, G4double NeutronEnergy, G4int lspin, G4int jspinx2 );
    G4int SampleJ( G4int theCompoundZ, G4int theCompoundA, G4int lspin );
    G4int GetAllowedJx2values( G4int theCompoundZ, G4int theCompoundA, G4int lspin, G4int* jx2vals );

    const G4NuDEXStatisticalNucleus* GetStatisticalNucleus( G4int za ) { return theStatisticalNucleus[za]; }
    G4int Init( G4int theZA, unsigned int seed1 = 0, unsigned int seed2 = 0, unsigned int seed3 = 0 );
    void SetBandWidth( G4double bandWidth ) { BandWidth = bandWidth; }
    void SetBrOption( G4int brOption ) { BrOption = brOption; }  

    G4NuDEXStatisticalNucleus* theStatisticalNucleus[G4NUDEX_MAXZA];
    G4int HasData[G4NUDEX_MAXZA];  // -1:no; 0:don't know; 1:yes
    G4String NuDEXLibDirectory;
    G4int BrOption;
    G4double BandWidth;

    G4int secID;  // creator model ID for the other secondaries produced by this model
    G4double lowestEnergyLimit;
    G4double minExcitation; 
    G4VEvaporationChannel* photonEvaporation;  // Needed when NuDEX is not applicable
};

#endif
