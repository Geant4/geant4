// -------------------------------------------------------------------
// $Id: MicrobeamRunAction.cc,v 1.1 2006-04-06 15:32:44 sincerti Exp $
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
      HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
    }
 
  numEvent = 0;
  nbOfHitsGas = 0;
    
  // ABSORBED DOSES INITIALIZATION
  DoseN = 0;
  DoseC = 0;
    
  massPhantom   = Detector->GetMassPhantom();
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
    HepRandom::showEngineStatus();
    HepRandom::saveEngineStatus("endOfRun.rndm");
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
