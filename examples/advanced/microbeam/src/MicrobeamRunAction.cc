//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------------
// $Id: MicrobeamRunAction.cc,v 1.5 2006-06-28 13:43:35 gunter Exp $
// -------------------------------------------------------------------

#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "MicrobeamRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamRunAction::MicrobeamRunAction(MicrobeamDetectorConstruction* det)
:Detector(det)
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
  //drawing
  if (G4VVisManager::GetConcreteInstance()) G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
       
  // save Rndm status
  if (saveRndm == 1)
  { 
    CLHEP::HepRandom::showEngineStatus();
    CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  }   
  
  FILE *myFile;
  myFile=fopen("3DDose.txt","a");
  for (G4int i=0; i<nbOfPixels; i++) 
  {  
    G4ThreeVector v;
    v = mapVoxels[i];
    if ( (GetNumEvent()+1) !=0) 
      fprintf (myFile,"%f %f %f %f \n", v.x(), v.y(), v.z(), dose3DDose[i]/(GetNumEvent()+1));
  }
  fclose (myFile);				
    
  G4cout << "-> Total number of particles detected by the gas detector : " << GetNbOfHitsGas() << G4endl;  
  G4cout << G4endl;    
}
