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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "MicrobeamRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamRunAction::MicrobeamRunAction(MicrobeamDetectorConstruction* det,
MicrobeamHistoManager * his)
:Detector(det),Histo(his)
{   
  saveRndm = 0;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamRunAction::~MicrobeamRunAction()
{
  delete[] dose3DDose;
  delete[] mapVoxels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamRunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{  
 
  // Histograms
  Histo->book();

  // save Rndm status
  if (saveRndm > 0)
    { 
      CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
    }
 
  numEvent = 0;
  nbOfHitsGas = 0;
    
  // ABSORBED DOSES INITIALIZATION
  DoseN = 0;
  DoseC = 0;
    
  massCytoplasm = Detector->GetMassCytoplasm();
  massNucleus   = Detector->GetMassNucleus();
  nbOfPixels = Detector->GetNbOfPixelsInPhantom();
  
  mapVoxels = new G4ThreeVector[nbOfPixels];
  dose3DDose = new G4float[nbOfPixels];

  for (G4int i=0; i<nbOfPixels; i++)
  {
  	mapVoxels [i]=myMicrobeamPhantomConfiguration.GetVoxelThreeVector(i);
  	dose3DDose[i]=0;
	G4ThreeVector v=mapVoxels[i];
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamRunAction::EndOfRunAction(const G4Run* /*aRun*/)
{     
  // save Rndm status
  if (saveRndm == 1)
  { 
    CLHEP::HepRandom::showEngineStatus();
    CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  }   
  
  for (G4int i=0; i<nbOfPixels; i++) 
  {  
    G4ThreeVector v;
    v = mapVoxels[i];
    if ( (GetNumEvent()+1) !=0) 
      {
	  Histo->FillNtuple(4,0,v.x());
	  Histo->FillNtuple(4,1,v.y());
	  Histo->FillNtuple(4,2,v.z());
	  Histo->FillNtuple(4,3,dose3DDose[i]/(GetNumEvent()+1));
	  Histo->AddRowNtuple(4);			     
      }
  }
   
  G4cout << "-> Total number of particles detected by the gas detector : " << GetNbOfHitsGas() << G4endl;  
  G4cout << G4endl;    
  
  //save histograms      
  Histo->save();

}
