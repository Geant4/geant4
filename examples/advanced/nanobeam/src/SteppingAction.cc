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
// $Id: SteppingAction.cc,v 1.2 2008/01/25 20:49:24 sincerti Exp $
// -------------------------------------------------------------------

#include "SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det,PrimaryGeneratorAction* pri)
:Run(run),Detector(det),Primary(pri)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* s)
  
{ 


if  (       (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")

/*
// for doublet

	 && (s->GetTrack()->GetPosition().z()>-3230.2)
         && (s->GetTrack()->GetPosition().z()<-3229.8) 
*/

// for triplet and whole line

	 && (s->GetTrack()->GetPosition().z()>249.99999)
         && (s->GetTrack()->GetPosition().z()<250.00001) 
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Vol")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World")
     )
		
     {
     	      xIn = s->GetTrack()->GetPosition().x();
	      yIn = s->GetTrack()->GetPosition().y();
	      zIn = s->GetTrack()->GetPosition().z();
	      E   = s->GetTrack()->GetKineticEnergy();

              G4ThreeVector angleIn;
              angleIn = s->GetTrack()->GetMomentumDirection();

              thetaIn = std::asin(angleIn[0]/std::sqrt(angleIn[0]*angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));
              phiIn = std::asin(angleIn[1]/std::sqrt(angleIn[0]*angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));

              G4cout << "    =>IMAGE : X(microns)=" << xIn/micrometer <<" Y(microns)="<< yIn/micrometer << " THETA(mrad)=" << (thetaIn/mrad) << " PHI(mrad)=" << (phiIn/mrad) << G4endl;
	      G4cout << G4endl;
	      		
   	      FILE *myFile1;
	      myFile1=fopen ("results/x.txt","a");
              fprintf(myFile1,"%e \n",xIn*1000);
              fclose (myFile1);

              FILE *myFile2;
   	      myFile2=fopen ("results/y.txt","a");
   	      fprintf(myFile2,"%e \n",yIn*1000);
	      fclose (myFile2);

	      FILE *myFile3;
	      myFile3=fopen ("results/theta.txt","a");
              fprintf(myFile3,"%e \n",thetaIn*1000);
              fclose (myFile3);

              FILE *myFile4;
              myFile4=fopen ("results/phi.txt","a");
              fprintf(myFile4,"%e \n",phiIn*1000);
              fclose (myFile4);

              FILE *myFile5;
              myFile5=fopen ("results/image.txt","a");
              fprintf(myFile5,"%e %e\n",xIn*1000, yIn*1000);
              fclose (myFile5);

     }

if (Detector->GetProfile()==1) 
{

	if  (
	    (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Vol")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Vol") )
	{
     	      xIn = s->GetTrack()->GetPosition().x();
	      yIn = s->GetTrack()->GetPosition().y();
	      zIn = s->GetTrack()->GetPosition().z();
              FILE *myFile6;
              myFile6=fopen ("results/profile.txt","a");
              fprintf(myFile6,"%e %e %e\n",xIn*1000, yIn*1000, zIn);
              fclose (myFile6);
	}

}
	
if (Detector->GetGrid()==1) 
{

	if  (
	    (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ControlVol_GridShadow")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World") )
	{
     	      xIn = s->GetTrack()->GetPosition().x();
	      yIn = s->GetTrack()->GetPosition().y();
              E   = s->GetTrack()->GetKineticEnergy();
	      FILE *myFile7;
              myFile7=fopen ("results/grid.txt","a");
              fprintf(myFile7,"%e %e %e\n",xIn*1000, yIn*1000, E);
              fclose (myFile7);
	}

}

// end
}     
