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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, CEA (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4AblaInterface_hh
#define G4AblaInterface_hh 1

#ifdef ABLAXX_IN_GEANT4_MODE

#include "G4VPreCompoundModel.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4Abla.hh"

class G4AblaInterface : public G4VPreCompoundModel {
public:
  G4AblaInterface();
  virtual ~G4AblaInterface();

  virtual G4ReactionProductVector *DeExcite(G4Fragment &aFragment);

  virtual G4HadFinalState *ApplyYourself(G4HadProjectile const &, G4Nucleus &) {
    return NULL;
  }

  virtual void ModelDescription(std::ostream& outFile) const;
  virtual void DeExciteModelDescription(std::ostream& outFile) const;

private:
  G4VarNtp *ablaResult;
  G4Volant *volant;
  G4Abla *theABLAModel;
  G4long eventNumber;

  /// \brief Convert an Abla particle to a G4ReactionProduct
  G4ReactionProduct *toG4Particle(G4int A, G4int Z , G4double kinE, G4double px, G4double py, G4double pz) const;

  /// \brief Convert A and Z to a G4ParticleDefinition
  G4ParticleDefinition *toG4ParticleDefinition (G4int A, G4int Z) const;

};

#endif // ABLAXX_IN_GEANT4_MODE

#endif
