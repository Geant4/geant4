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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4WeightWindowBiasing_h
#define G4WeightWindowBiasing_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4GeometrySampler.hh"
#include "G4WeightWindowAlgorithm.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4WeightWindowBiasing : public G4VPhysicsConstructor
{
public:

  G4WeightWindowBiasing(const G4String& name = "NoParallelWP");
  G4WeightWindowBiasing(G4GeometrySampler* mgs, G4VWeightWindowAlgorithm* wwAlg, G4PlaceOfAction placeOfAction, const G4String& name = "NoParallelWP");
  virtual ~G4WeightWindowBiasing();

public:

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

   // hide assignment operator
  G4WeightWindowBiasing & operator=(const G4WeightWindowBiasing &right);
  G4WeightWindowBiasing(const G4WeightWindowBiasing&);

  G4GeometrySampler* fGeomSampler;
  G4VWeightWindowAlgorithm* fWWalg;
  G4PlaceOfAction fPlaceOfAction;

  G4bool paraFlag;
  G4String paraName;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
