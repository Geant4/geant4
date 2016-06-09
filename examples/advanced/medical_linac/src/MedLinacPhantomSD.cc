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
//// $Id: MedLinacPhantomSD.cc,v 1.7 2006/06/29 16:04:37 gunter Exp $
//
//
// Code developed by: M. Piergentili

//
//
#include "MedLinacPhantomSD.hh"
#include "MedLinacPhantomMessenger.hh"
#include "MedLinacPhantomHit.hh"
#include "MedLinacAnalysisManager.hh"
#include "MedLinacPhantomMessenger.hh"
#include "MedLinacDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

//....

MedLinacPhantomSD::MedLinacPhantomSD(G4String name):G4VSensitiveDetector(name)
{
 phantomMessenger = new MedLinacPhantomMessenger(this);
}

MedLinacPhantomSD::~MedLinacPhantomSD()
{
  delete phantomMessenger;
}

void MedLinacPhantomSD::Initialize(G4HCofThisEvent*)
{
}

void MedLinacPhantomSD::SetPhantomDimension (G4double val)
{
  phantomDimension = val;
  //G4cout <<"2==============================phantomDim  "<< phantomDimension/mm<<"mm"<<G4endl;
}

void MedLinacPhantomSD::SetNumberOfPhantomVoxels (G4int val)
{
  numberOfPhantomVoxels = val;
  //G4cout <<"2==============================numberOfVoxels  "<< numberOfPhantomVoxels<<G4endl;
}

G4bool MedLinacPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "Phantom_phys")
    return false;

  G4double energyDep = aStep->GetTotalEnergyDeposit();
  if(energyDep == 0.)
    return false;

    // Read Voxel indexes: i is the x index, k is the z index
  G4int k = ROhist->GetReplicaNumber(1);
  G4int i = ROhist->GetReplicaNumber(2);
  G4int j = ROhist->GetReplicaNumber();
  
  //G4int numberOfVoxelZ = 150;
  //G4double voxelWidthZ = 2. *mm;
  G4int numberOfVoxelZ = numberOfPhantomVoxels;
  G4double voxelWidthZ = phantomDimension/numberOfVoxelZ;
  //G4cout <<"########## phantomDimension  "<< phantomDimension/mm<<"mm"<<G4endl;
  //G4cout <<"########## numberOfVoxelZ "<<numberOfVoxelZ <<G4endl;
  //G4cout <<"########## voxelWidthZ "<<voxelWidthZ <<G4endl;

  G4double x = (-numberOfVoxelZ+1+2*i)*voxelWidthZ; 
  G4double y = (- numberOfVoxelZ+1+2*j)*voxelWidthZ;
  G4double z = (- numberOfVoxelZ+1+2*k)*voxelWidthZ;

 if(energyDep != 0)                       
	    {       
	     #ifdef G4ANALYSIS_USE	
                 MedLinacAnalysisManager* analysis = 
                                      MedLinacAnalysisManager::getInstance();   
//**** PDD in isocenter (Y and X Thickness = 5. mm) *********** 	   
	  if(energyDep != 0)                       
	    { 
	      if (y<=2.5*mm){if (y>= -2.5*mm)
		{
		  if(x<=2.5*mm){if (x>= -2.5*mm)
		{analysis->FillHistogram1WithEnergy(z,energyDep/MeV);}
		  }
		}
	      }
	    }
//***** flatness  along x ***** Depth=build-up (15 mm)*************************
  if(energyDep != 0)                       
    {  
       if (z<=137.5*mm){if (z>= 132.5*mm) 
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram2WithEnergy(x,energyDep/MeV);}
	 }
	 }
       }
    }

//***** flatness  along x ***** Depth=50mm 5*********************************
  if(energyDep != 0)                       
    {  
      if (z<=102.5*mm){if (z>= 97.5*mm)
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram3WithEnergy(x,energyDep/MeV);}
	 }}}}

//***** flatness  along x ***** Depth=100 mm********************************
  if(energyDep != 0)
    {
      if (z<=52.5*mm){if (z>= 47.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram4WithEnergy(x,energyDep/MeV);}                                                             
                        }}}}
//***** flatness  along x ***** Depth=200 mm********************************
  if(energyDep != 0)
    {
      if (z<=-47.5*mm){if (z>= -52.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram5WithEnergy(x,energyDep/MeV);} 
                        }}}}
//**** PDD in isocenter (Y and X Thickness = 5. mm) *********** 	   
	  if(energyDep != 0)                       
	    { 
	      if (y<=2.5*mm){if (y>= -2.5*mm)
		{
		  if(x<=2.5*mm){if (x>= -2.5*mm)
		{analysis->FillHistogram6WithEnergy(z,energyDep/MeV);}
		  }}}}

//***** flatness  along x ***** Depth=build-up (15 mm)*************************
  if(energyDep != 0)                       
    {  
       if (z<=137.5*mm){if (z>= 132.5*mm) 
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram7WithEnergy(x,energyDep/MeV);}
	                }}}}
//***** flatness  along x ***** Depth=50mm*********************************
  if(energyDep != 0)                       
    {  
      if (z<=102.5*mm){if (z>= 97.5*mm)
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram8WithEnergy(x,energyDep/MeV);}
	                }}}}
//***** flatness  along x ***** Depth=100 mm********************************
  if(energyDep != 0)                       
    {  
      if (z<=52.5*mm){if (z>= 47.5*mm) 
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram9WithEnergy(x,energyDep/MeV);}
	                }}}}
//***** flatness  along x ***** Depth=200 mm********************************
  if(energyDep != 0)                       
    {  
      if (z<=-47.5*mm){if (z>= -52.5*mm)
	 { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram10WithEnergy(x,energyDep/MeV);}
	                }
	 }
                     }
    }

//**** PDD in isocenter (Y and X Thickness = 5. mm) ***********
          if(energyDep != 0)
            {
              if (y<=2.5*mm){if (y>= -2.5*mm)
                {
                  if(x<=2.5*mm){if (x>= -2.5*mm)
                {analysis->FillHistogram11WithEnergy(z,energyDep/MeV);}
                  }
                }
              }
            }
//***** flatness  along x ***** Depth=build-up (15 mm)*************************
  if(energyDep != 0)
    {
       if (z<=137.5*mm){if (z>= 132.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram12WithEnergy(x,energyDep/MeV);}
         }
         }
       }
    }              
//***** flatness  along x ***** Depth=50mm 5*********************************
  if(energyDep != 0)
    {
      if (z<=102.5*mm){if (z>= 97.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram13WithEnergy(x,energyDep/MeV);}
         }}}}
//***** flatness  along x ***** Depth=100 mm********************************
  if(energyDep != 0)
    {
      if (z<=52.5*mm){if (z>= 47.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram14WithEnergy(x,energyDep/MeV);}
                        }}}}
//***** flatness  along x ***** Depth=200 mm********************************
  if(energyDep != 0)
    {
      if (z<=-47.5*mm){if (z>= -52.5*mm)
         { if (y<=2.5*mm){if (y>= -2.5*mm)
       {analysis->FillHistogram15WithEnergy(x,energyDep/MeV);}
     }}}}

#endif 	       
	    }
return true;
}

void MedLinacPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
}

void MedLinacPhantomSD::clear()
{
} 

void MedLinacPhantomSD::DrawAll()
{
}

void MedLinacPhantomSD::PrintAll()
{
}

