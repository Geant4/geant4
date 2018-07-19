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
// $Id: G4ee2KChargedModel.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ee2KChargedModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 09.07.2008
//
// Modifications:
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4ee2KChargedModel_h
#define G4ee2KChargedModel_h 1

#include "G4Vee2hadrons.hh"
#include "globals.hh"
#include "G4eeCrossSections.hh"

class G4DynamicParticle;
class G4PhysicsVector;

class G4ee2KChargedModel : public G4Vee2hadrons
{

public:

  explicit G4ee2KChargedModel(G4eeCrossSections*,G4double,G4double);

  virtual ~G4ee2KChargedModel();

  virtual G4double PeakEnergy() const override;

  virtual G4double ComputeCrossSection(G4double) const override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
              G4double, const G4ThreeVector&) override;

private:

  // hide assignment operator
  G4ee2KChargedModel & operator=(const  G4ee2KChargedModel &right);
  G4ee2KChargedModel(const  G4ee2KChargedModel&);

  G4double massK;
  G4double massPhi;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
