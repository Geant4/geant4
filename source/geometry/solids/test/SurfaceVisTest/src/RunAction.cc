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

#include "RunAction.hh"
#include "RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4VSolid.hh"
#include "globals.hh"
#include <vector>
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

#include "G4ios.hh"
#include <iomanip>

#include "Randomize.hh"

//////////////////////////////////////////////////////////////////////////////

RunAction::RunAction()
  : saveRndm(0)
{
  runMessenger = new RunMessenger(this);
  saveNum=1000;
}

////////////////////////////////////////////////////////////////////////////

RunAction::~RunAction()
{
 delete runMessenger;
}

/////////////////////////////////////////////////////////////////////////////

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
  { 
      CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
  }  
   G4int check=CheckPoints(saveNum);
   if(check==0){G4cout<<"Check Point On Surface Ok! for "<<saveNum<<" points"<<G4endl;}
   else{G4cout<<"Check Point On Surface has "<<check<<" inconsistencies"<<G4endl;};
   DrawPoints(saveNum);
        


}

/////////////////////////////////////////////////////////////////////////////

void RunAction::EndOfRunAction(const G4Run*)
{
  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
  // DrawPoints(saveNum);
  // save Rndm status

  if (saveRndm == 1)
  { 
    CLHEP::HepRandom::showEngineStatus();
    CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  }     
}

//
//
////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void RunAction::DrawPoints(G4int num)
{
  //LogicalVolumeStore search for Solid 
         G4LogicalVolumeStore *theLVStore = G4LogicalVolumeStore::GetInstance();
         G4VSolid *solid;
         G4ThreeVector point;
         solid=theLVStore->GetVolume("aVolume_L",false)->GetSolid();
	
	  if(solid){
          //Loop for Points on the Surface
          G4int N=num;
          for (G4int j=0; j<N;j++){
           //Visualisation
            G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
            if(pVVisManager){
             point=solid->GetPointOnSurface();
	     G4Circle circle(point); // Instantiate a circle with its 3D
                            // position. The argument "position"
                            // is defined as G4Point3D instance
             circle.SetScreenDiameter (2.0); // Should be circle.SetScreenDiameter
                                 //  (1.0 * pixels) - to be implemented
             circle.SetFillStyle (G4Circle::filled); // Make it a filled circle
             G4Colour colour(1.,0.,0.);              // Define red color
             G4VisAttributes attribs(colour);        // Define a red visualization attribute
             circle.SetVisAttributes(attribs);       //
          // make a 3D data for visualization 
          pVVisManager->Draw(circle);               
          }
         }
	  }else{//no corresponding solid found
	    G4cout<<"RunAction::DrawPoints Sorry, No solid is found "<<G4endl;
   
	  }
}
// Check if the Randomly Choosen Points are On the Surface 
// return number of points not on the surface
G4int RunAction::CheckPoints(G4int num)
{
  //LogicalVolumeStore search for Solid 
         G4LogicalVolumeStore *theLVStore = G4LogicalVolumeStore::GetInstance();
         G4VSolid *solid;
         G4ThreeVector point;
         EInside surface;
         solid=theLVStore->GetVolume("aVolume_L",false)->GetSolid();
	  G4int count=0;
	  if(solid){
          //Loop for Points on the Surface
          G4int N=num;
          for (G4int j=0; j<N;j++){
          
             point=solid->GetPointOnSurface();
             surface=solid->Inside(point);
             if(surface != kSurface)count++;
	   }
	  }else{//no corresponding solid found
	    G4cout<<"RunAction::DrawPoints Sorry, No solid is found "<<G4endl;
            return count;
	  }
	  
 return count;
}
//
//
////////////////////////////////////////////////////////////////////////
 void RunAction:: SetNumberOfSurfacePoints(G4int val)
{saveNum=val;
 
}
//
//
////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
