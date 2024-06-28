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
//      Geant4 header file 
//
//      File name: G4ParticleHPInelasticURR.hh
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//	         Loic Thulliez (CEA France)
// 
//      Creation date: 4 June 2024
//
//      Description: Final state production model for a high precision
//                   (based on evaluated data libraries) description of
//                   neutron Inelastic scattering below 20 MeV.
//                   To be used in your physics list in case you need
//                   this physics.
//                   In this case you want to register an object of this
//                   class with the corresponding process.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//
#ifndef G4ParticleHPInelasticURR_h
#define G4ParticleHPInelasticURR_h 1

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include <vector>

class G4ParticleHPInelastic;


class G4ParticleHPInelasticURR : public G4HadronicInteraction {
  public:
    G4ParticleHPInelasticURR();
    ~G4ParticleHPInelasticURR();
  
    G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );

    virtual const std::pair< G4double, G4double > GetFatalEnergyCheckLevels() const;
    G4int GetVerboseLevel() const;
    void SetVerboseLevel( G4int );
    void BuildPhysicsTable( const G4ParticleDefinition& );
    virtual void ModelDescription( std::ostream& outFile ) const;

  private:
    G4ParticleHPInelastic* particleHPinelastic;
    std::vector< std::pair< G4double, G4double > >* URRlimits{ nullptr };
    G4bool doNOTusePTforInelastic{ true };
};

#endif
