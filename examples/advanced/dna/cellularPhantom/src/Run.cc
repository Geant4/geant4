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
// -----------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet (a), S. Incerti (a), N. H. Tran (a), L. Morelli (a,b)
//
// a) University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
// b) Politecnico di Milano, Italy
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// -----------------------------------------------------------------------------

#include "Run.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
:G4Run()
{
  G4int nbVoxel = CellParameterisation::Instance()->GetPhantomTotalPixels();
  fVoxelEdeposit = new G4double[nbVoxel];
  for (G4int i=0; i<nbVoxel; ++i) fVoxelEdeposit[i] = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
  delete[] fVoxelEdeposit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // Accumulate energy deposits per voxel
  G4int nbVoxel = CellParameterisation::Instance()->GetPhantomTotalPixels();
  for (G4int i=0; i<nbVoxel; ++i)
    fVoxelEdeposit[i] += localRun->fVoxelEdeposit[i];

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4double nrjRed=0;
  G4double doseRed=0;
  G4double nrjGreen=0;
  G4double doseGreen=0;
  G4double nrjBlue=0;
  G4double doseBlue=0;

  fMyPhantomParam = CellParameterisation::Instance();

  G4double redMassTot = fMyPhantomParam->GetRedMass();
  G4double greenMassTot = fMyPhantomParam->GetGreenMass();
  G4double blueMassTot = fMyPhantomParam->GetBlueMass();

  for (G4int i = 0; i < fMyPhantomParam->GetPhantomTotalPixels(); i++)
  {
    if (fVoxelEdeposit[i] > 0.)
    {
      if (fMyPhantomParam->GetMaterial(i) == 1)
      {
        nrjRed=nrjRed+fVoxelEdeposit[i];
        doseRed=doseRed+(fVoxelEdeposit[i]/redMassTot);
      }
      else if (fMyPhantomParam->GetMaterial(i) == 2)
      {
        nrjGreen=nrjGreen+fVoxelEdeposit[i];
        doseGreen=doseGreen+(fVoxelEdeposit[i]/greenMassTot);
      }
      else if (fMyPhantomParam->GetMaterial(i) == 3)
      {
        nrjBlue=nrjBlue+fVoxelEdeposit[i];
        doseBlue=doseBlue+(fVoxelEdeposit[i]/blueMassTot);
      }
    }
  }

  G4int numberVoxTot = fMyPhantomParam->GetPhantomTotalPixels();
  G4int numberVoxRed = fMyPhantomParam->GetRedTotalPixels();
  G4int numberVoxGreen = fMyPhantomParam->GetGreenTotalPixels();
  G4int numberVoxBlue = fMyPhantomParam->GetBlueTotalPixels();

  G4cout << G4endl;
  G4cout << "- Summary --------------------------------------------------" << G4endl;
  G4cout << G4endl;
  G4cout << "  Total number of voxels in phantom           = " << numberVoxTot << G4endl;
  G4cout << "  Total number of RED voxels in phantom       = " << numberVoxRed << G4endl;
  G4cout << "  Total number of GREEN voxels in phantom     = " << numberVoxGreen << G4endl;
  G4cout << "  Total number of BLUE voxels in phantom      = " << numberVoxBlue << G4endl;
  G4cout << G4endl;
  G4cout << "  Total absorbed energy in RED   voxels (MeV) = "  << nrjRed/MeV << G4endl;
  G4cout << "  Total absorbed energy in GREEN voxels (MeV) = "  << nrjGreen/MeV << G4endl;
  G4cout << "  Total absorbed energy in BLUE  voxels (MeV) = "  << nrjBlue/MeV << G4endl;
  G4cout << G4endl;
  G4cout << "  Total absorbed dose   in RED   voxels (Gy)  = "  << doseRed/(joule/kg) << G4endl;
  G4cout << "  Total absorbed dose   in GREEN voxels (Gy)  = "  << doseGreen/(joule/kg) << G4endl;
  G4cout << "  Total absorbed dose   in BLUE  voxels (Gy)  = "  << doseBlue/(joule/kg) << G4endl;
  G4cout << G4endl;
  G4cout << "------------------------------------------------------------" << G4endl;

}
