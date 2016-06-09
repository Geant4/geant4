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
//
// $Id: G4BlinePrimaryGeneratorAction.cc,v 1.2 2004/12/03 16:07:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
//
// --------------------------------------------------------------------
//
// G4BlinePrimaryGeneratorAction implementation
//
// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------

#include "G4BlinePrimaryGeneratorAction.hh"

#include "G4Types.hh"
#include "G4Event.hh"
#include "G4ChargedGeantino.hh"

///////////////////////////////////////////////////////////////////////////

G4BlinePrimaryGeneratorAction::G4BlinePrimaryGeneratorAction()
{
  fUserPrimaryAction = 0;
  FirstPartOfBline = true;
}

///////////////////////////////////////////////////////////////////////////

G4BlinePrimaryGeneratorAction::~G4BlinePrimaryGeneratorAction()
{ 
}

///////////////////////////////////////////////////////////////////////////

void G4BlinePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (!fUserPrimaryAction)
  {
    G4Exception("G4BlinePrimaryGeneratorAction::GeneratePrimaries()",
                "NullPointer", JustWarning,
                "Primary generator action not defined !");
    return;
  }

  // For the first part of a bline the start position and time are defined
  // by using the USER primary action while for the second part the previous
  // values are taken.
  
  if (FirstPartOfBline)
  {
    // set the position and time defined by using the USER primary action
 
    G4Event* tmpEvent = new G4Event();    
    fUserPrimaryAction->GeneratePrimaries(tmpEvent);
    BlineStartPosition = tmpEvent->GetPrimaryVertex()->GetPosition();
    T0 = tmpEvent->GetPrimaryVertex()->GetT0();
    delete tmpEvent;
  }
  FirstPartOfBline = false;   

  G4PrimaryVertex* primary_vertex = 
    new G4PrimaryVertex(BlineStartPosition, T0);

  // Define the particle to be tracked as Charged Geantino
    
  G4ChargedGeantino* pdef = G4ChargedGeantino::ChargedGeantino();
   
  G4double mass =  pdef->GetPDGMass();
  G4double energy = 10000.*MeV + mass;
  G4double pmom = std::sqrt(energy*energy-mass*mass);

  // The momentum direction and energy do not have an effect in tracing of 
  // bline but still need to be defined.
   
  G4double px = 0.;
  G4double py = 0.;
  G4double pz = pmom;

  G4PrimaryParticle* particle = new G4PrimaryParticle(pdef,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge(pdef->GetPDGCharge());
  primary_vertex->SetPrimary( particle );
 
  anEvent->AddPrimaryVertex( primary_vertex );
}
