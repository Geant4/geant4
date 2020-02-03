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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTStepMax.hh"
#include "IORTStepMaxMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
IORTStepMax::IORTStepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),MaxChargedStep(DBL_MAX)
{
  pMess = new IORTStepMaxMessenger(this);
}
 
/////////////////////////////////////////////////////////////////////////////
IORTStepMax::~IORTStepMax() { delete pMess; }

/////////////////////////////////////////////////////////////////////////////
G4bool IORTStepMax::IsApplicable(const G4ParticleDefinition& particle) 
{ 
  return (particle.GetPDGCharge() != 0.);
}

/////////////////////////////////////////////////////////////////////////////    
void IORTStepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}

/////////////////////////////////////////////////////////////////////////////
G4double IORTStepMax::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double ProposedStep = DBL_MAX;

  if((MaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() == "DetectorPhys"))
     ProposedStep = MaxChargedStep;
  /*                              IMPORTED FROM eliot_geant4.9.3p01_version
 if((MaxChargedStep > 0.) &&       // ATTENZIONE 1) CI POSSONO ESSERE + DI 2 DISCHI
     (aTrack.GetVolume() != 0) &&  // ATTENZIONE 2) ANCORA NON HO INSERITO I DISCHI CHE ANDREBBERO NEL PASSPROTONBL
     (aTrack.GetVolume()->GetName() == "DiscoIORTPhys"))
     ProposedStep = MaxChargedStep;

 if((MaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() == "DiscoIORTPhys1"))
     ProposedStep = MaxChargedStep;
    */
  return ProposedStep;
}

/////////////////////////////////////////////////////////////////////////////
G4VParticleChange* IORTStepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

