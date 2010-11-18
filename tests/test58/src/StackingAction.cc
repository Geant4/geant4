
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
// StackingAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


StackingAction::StackingAction(RunAction* RA, EventAction* EA,
				DetectorConstruction* DE )
:runaction(RA), eventaction(EA), detector(DE)
{
  killSecondary  = 0;
  stackMessenger = new StackingMessenger(this);
  rec=0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




G4double  StackingAction::DamageEnergy(G4double T,G4double A, G4double Z)
{

//.................. T in  eV!!!!!!!!!!!!!
 G4double Z2= Z;
 G4double M2= A*mole/g;
 G4double k_d;
 G4double epsilon_d;
 G4double g_epsilon_d;
 G4double E_nu;

   k_d=0.1334*std::pow(Z2,(2./3.))*std::pow(M2,(-1./2.));
    epsilon_d=0.01014*std::pow(Z2,(-7./3.))*(T/eV);
    g_epsilon_d= epsilon_d+0.40244*std::pow(epsilon_d,(3./4.))+3.4008*std::pow(epsilon_d,(1./6.));



    E_nu=1./(1.+ k_d*g_epsilon_d);

return E_nu;//partition fraction!!!
}


//...................................................................


G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{

G4int IDp= aTrack->GetParentID();
  //keep primary particle
  if (IDp == 0) {return fUrgent;} 
      

 //
  //energy spectrum of secondaries
  //
  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();

		//NIEL 
  if (IDp==1 && charge  !=0){ 

 	G4Material*  material= detector->GetAbsorberMaterial();

	G4double A2  = material->GetA();
        G4double Z2  = material->GetZ();

 	        runaction->NumberRec(1);
                rec++;
	
		//Lindhard partition                
		G4double LT = DamageEnergy(energy,A2,Z2);
	      
		//total sum of T*L(T)
		runaction->SumTL(energy*LT); 

                //total sum of T
 
	//	G4cout<<"T "<<energy/eV<<G4endl;
               runaction->SumT(energy);




			}

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (killSecondary) 
    {if (killSecondary == 1) 

	//add secondary energy before killing
     	eventaction->AddEnergy(energy);
        eventaction->AddNonIonEnergy(energy);

     	status = fKill;}

  return status;
  

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
