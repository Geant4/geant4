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
// $Id: HadrontherapySteppingAction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyDetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyRunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include <iomanip.h>
#include "G4UImanager.hh"

// -----------------------------------------------------------------
HadrontherapySteppingAction::HadrontherapySteppingAction(HadrontherapyEventAction* EA)
{
  eventaction = EA;
}

// ----------------------------------------------------------------
HadrontherapySteppingAction::~HadrontherapySteppingAction()
{
}

// -----------------------------------------------------------------
void HadrontherapySteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  event_id = eventaction -> Trasporto();
  //the function Trasporto takes the event's number from EventAction Class
 
G4double x  = aStep -> GetPostStepPoint() -> GetPosition().x();
G4double y  = aStep -> GetPostStepPoint() -> GetPosition().y();
G4double z  = aStep -> GetPostStepPoint() -> GetPosition().z();
  
G4Track * theTrack = aStep -> GetTrack();
G4double TrackID = theTrack -> GetTrackID();
G4double KEnergy = theTrack -> GetKineticEnergy();
  
G4double thicknessPlane = 1 *mm;
G4double xPlane = 1048.59 *mm;


//------------------------------------------------------------
// The PLANE that simulate radiochromic film 
// Particles that interact within this "software plane" are registered;

 if (event_id == 0) {Controllo = 0;}
 
 if (x >= xPlane & x <= (xPlane+thicknessPlane) & Controllo != event_id & TrackID == 1)
   {
     Controllo = event_id;
// ----------------------------------------------------------
// WRITE ASCII FILES
// Properties of the particles interacting with the plane
// are registered in two different files. Exactly 
// Number of event and kinetic energy are registered in the disEnXX.dat file;
// Coordinates of the particles are registered in the disAngXX.dat file;
     
     std::ofstream pmtfile("energyDistribution.out", std::ios::app);
     if(pmtfile.is_open())
       
       {
	 pmtfile << KEnergy  << '\t' << "   " << event_id << '\t' << G4endl;
       }
     
     std::ofstream pmtfile2("angularDistribution.out", std::ios::app);
     if(pmtfile2.is_open())
       
       {
	 pmtfile2 << x/mm << '\t' <<"  " << y/mm << '\t' <<"  " << z/mm << '\t' << G4endl;
       }
   }
}





