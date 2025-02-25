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
// --------------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet, S. Incerti, N. H. Tran, L. Morelli
//
// University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// --------------------------------------------------------------------------------

#include "SteppingAction.hh"

#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* runAction)
:G4UserSteppingAction(), fRunAction(runAction)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // ********************************************************************************
  // Avoid string comparison to extract material (1, 2 or 3) whic causes issues in MT
  // ********************************************************************************

  fMyPhantomParam = CellParameterisation::Instance();
  const G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4int preReplicaNumber = preStep->GetTouchableHandle()->GetReplicaNumber();
  G4int voxelMaterial =  fMyPhantomParam->GetMaterial(preReplicaNumber);

  // The absorbed energy is added to the "voxel energy" array in RunAction
  // Added protection to make sure Replica Number has been identified

  if (aStep->GetTotalEnergyDeposit()>0. && preReplicaNumber>0)
  {
    if (voxelMaterial == 1)
    {
      fRunAction->AddDoseBox(preReplicaNumber, aStep->GetTotalEnergyDeposit());
    }
    else if (voxelMaterial == 2)
    {
      fRunAction->AddDoseBox(preReplicaNumber, aStep->GetTotalEnergyDeposit());
    }
    else if (voxelMaterial == 3)
    {
      fRunAction->AddDoseBox(preReplicaNumber, aStep->GetTotalEnergyDeposit());
    }
  }
}
