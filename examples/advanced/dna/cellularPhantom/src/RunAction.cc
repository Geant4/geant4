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

#include "RunAction.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
:G4UserRunAction()
{
  if (isMaster)
  {
    // Declare ntuples
    auto man = G4AnalysisManager::Instance();
    man->SetDefaultFileType("root");
    man->SetFirstNtupleId(1);

    // Create 1st ntuple (id = 1)
    man->CreateNtuple("ntuple1", "RED");
    man->CreateNtupleDColumn("x");
    man->CreateNtupleDColumn("y");
    man->CreateNtupleDColumn("z");
    man->CreateNtupleDColumn("energy");
    man->CreateNtupleDColumn("dose");
    man->CreateNtupleIColumn("voxelID");
    man->FinishNtuple();

    // Create 2nd ntuple (id = 2)
    man->CreateNtuple("ntuple2", "GREEN");
    man->CreateNtupleDColumn("x");
    man->CreateNtupleDColumn("y");
    man->CreateNtupleDColumn("z");
    man->CreateNtupleDColumn("energy");
    man->CreateNtupleDColumn("dose");
    man->CreateNtupleIColumn("voxelID");
    man->FinishNtuple();

    // Create 3rd ntuple (id = 3)
    man->CreateNtuple("ntuple3", "BLUE");
    man->CreateNtupleDColumn("x");
    man->CreateNtupleDColumn("y");
    man->CreateNtupleDColumn("z");
    man->CreateNtupleDColumn("energy");
    man->CreateNtupleDColumn("dose");
    man->CreateNtupleIColumn("voxelID");
    man->FinishNtuple();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  fRun = new Run();
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run *)
{
  if (isMaster)
  {
    // Analysis manager
    auto man = G4AnalysisManager::Instance();
    man->OpenFile("phantom");

    // Access phantom singleton
    fMyPhantomParam = CellParameterisation::Instance();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run * /*aRun*/)
{
  if (isMaster)
  {
    // Display results from merged local runs
    fRun->EndOfRun();

    // Fill ntuples
    auto man = G4AnalysisManager::Instance();

    G4double X, Y, Z;

    // Total mass of voxels
    G4double redMassTot=0.;
    G4double greenMassTot=0.;
    G4double blueMassTot=0.;

    redMassTot = fMyPhantomParam->GetRedMass();
    greenMassTot = fMyPhantomParam->GetGreenMass();
    blueMassTot = fMyPhantomParam->GetBlueMass();

    // (Optional) Numbers of voxel
    //G4double redVox=0;
    //G4double greenVox=0;
    //G4double blueVox=0;
    //redVox = fMyPhantomParam->GetRedTotalPixels();
    //greenVox = fMyPhantomParam->GetGreenTotalPixels();
    //blueVox = fMyPhantomParam->GetBlueTotalPixels();

    // (Optional) Single voxel mass
    //G4double redMass=0.;
    //G4double greenMass=0.;
    //G4double blueMass=0.;
    //redMass = redMassTot/redVox;
    //greenMass = greenMassTot/greenVox;
    //blueMass = blueMassTot/blueVox;

    // Save x, y, z and energy for every voxel having absorbed an energy above 0.
    // Energy is in keV
    // Dose is in Gy

    G4double cumulatedDeposit = 0;

    // Loop on voxels and collect energy and dose from merged local runs
    for (G4int i = 0; i < fMyPhantomParam->GetPhantomTotalPixels(); i++)
    {
      cumulatedDeposit = fRun->GetVoxelEdeposit(i);

      if (cumulatedDeposit > 0.)
      {
        X = (fMyPhantomParam->GetVoxelThreeVectorOriginal(i).x()) / um;
        Y = (fMyPhantomParam->GetVoxelThreeVectorOriginal(i).y()) / um;
        Z = (fMyPhantomParam->GetVoxelThreeVectorOriginal(i).z()) / um;

        if (fMyPhantomParam->GetMaterial(i) == 1)
        {
          man->FillNtupleDColumn(1,0,X);
          man->FillNtupleDColumn(1,1,Y);
          man->FillNtupleDColumn(1,2,Z);
          man->FillNtupleDColumn(1,3,cumulatedDeposit/keV);
          man->FillNtupleDColumn(1,4,((cumulatedDeposit/joule)/(redMassTot/kg)));
          man->FillNtupleIColumn(1,5,i);
          man->AddNtupleRow(1);
        }
        else if (fMyPhantomParam->GetMaterial(i) == 2)
        {
          man->FillNtupleDColumn(2,0,X);
          man->FillNtupleDColumn(2,1,Y);
          man->FillNtupleDColumn(2,2,Z);
          man->FillNtupleDColumn(2,3,cumulatedDeposit/keV);
          man->FillNtupleDColumn(2,4,((cumulatedDeposit/joule)/(greenMassTot/kg)));
          man->FillNtupleIColumn(2,5,i);
          man->AddNtupleRow(2);
        }
        else if (fMyPhantomParam->GetMaterial(i) == 3)
        {
          man->FillNtupleDColumn(3,0,X);
          man->FillNtupleDColumn(3,1,Y);
          man->FillNtupleDColumn(3,2,Z);
          man->FillNtupleDColumn(3,3,cumulatedDeposit/keV);
          man->FillNtupleDColumn(3,4,((cumulatedDeposit/joule)/(blueMassTot/kg)));
          man->FillNtupleIColumn(3,5,i);
          man->AddNtupleRow(3);
        }
      }
    }

    // Save histograms
    man->Write();
    man->CloseFile();
    man->Clear();
  }
}
