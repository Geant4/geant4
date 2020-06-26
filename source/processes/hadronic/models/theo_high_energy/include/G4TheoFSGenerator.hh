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
//

// Class Description
// Final state production model for theoretical models of hadron inelastic
// scattering in geant4;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Note: This class is part of an implementation framework. You need to
// register corresponding high energy generators and transport codes to 
// fill it with life; decay of strong resonances is done directly,
// in case there is no residual nucleus. 
// Class Description - End

#ifndef G4TheoFSGenerator_h
#define G4TheoFSGenerator_h 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4VHighEnergyGenerator.hh"
#include "G4DecayStrongResonances.hh"
#include "G4HadFinalState.hh"
#include "G4QuasiElasticChannel.hh"

class G4TheoFSGenerator : public G4HadronicInteraction 

{
  public:
      explicit G4TheoFSGenerator(const G4String& name = "TheoFSGenerator");
      ~G4TheoFSGenerator() override;

      G4TheoFSGenerator(const G4TheoFSGenerator &right) = delete;
      const G4TheoFSGenerator & operator=(const G4TheoFSGenerator &right) = delete;
      G4bool operator==(const G4TheoFSGenerator &right) const = delete;
      G4bool operator!=(const G4TheoFSGenerator &right) const = delete;

  public:
      G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus) override;
      inline void SetTransport(G4VIntraNuclearTransportModel *const  value);
      inline void SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value);
      inline void SetQuasiElasticChannel(G4QuasiElasticChannel *const value);

      std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const override;
      void ModelDescription(std::ostream& outFile) const override;

  private: 
      G4VIntraNuclearTransportModel * theTransport;
      G4VHighEnergyGenerator * theHighEnergyGenerator;
      G4DecayStrongResonances theDecay;
      G4HadFinalState * theParticleChange;
      G4QuasiElasticChannel * theQuasielastic;
};

inline void G4TheoFSGenerator::SetTransport(G4VIntraNuclearTransportModel *const  value)
{
  theTransport = value;
}

inline void G4TheoFSGenerator::SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value)
{
  theHighEnergyGenerator= value;
}

inline void G4TheoFSGenerator::SetQuasiElasticChannel(G4QuasiElasticChannel *const value)
{
  theQuasielastic = value;
}

#endif


