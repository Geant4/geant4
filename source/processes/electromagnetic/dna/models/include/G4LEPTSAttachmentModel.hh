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
#ifndef G4LEPTSAttachmentModel_h
#define G4LEPTSAttachmentModel_h 1

#include "G4VLEPTSModel.hh"
#include "G4ParticleChangeForGamma.hh"


class G4LEPTSAttachmentModel : public G4VLEPTSModel
{ 
public:
  G4LEPTSAttachmentModel(const G4String& processName ="G4LEPTSAttachmentModel");
  ~G4LEPTSAttachmentModel() override;

  void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin = 0.0,
                                 G4double tmax = DBL_MAX) override;

 // main method to compute cross section per Volume
  G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX) override;

private:
  G4ParticleChangeForGamma* fParticleChangeForGamma;

};


#endif
