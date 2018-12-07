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
/// \file field/BlineTracer/src/G4BlinePrimaryGeneratorAction.cc
/// \brief Implementation of the G4BlinePrimaryGeneratorAction class
//
//
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
#include "G4SystemOfUnits.hh"

///////////////////////////////////////////////////////////////////////////

G4BlinePrimaryGeneratorAction::G4BlinePrimaryGeneratorAction()
{
  fUserPrimaryAction = 0;
  fFirstPartOfBline = true;
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
  
  if (fFirstPartOfBline)
  {
    // set the position and time defined by using the USER primary action
 
    G4Event* tmpEvent = new G4Event();    
    fUserPrimaryAction->GeneratePrimaries(tmpEvent);
    fBlineStartPosition = tmpEvent->GetPrimaryVertex()->GetPosition();
    fT0 = tmpEvent->GetPrimaryVertex()->GetT0();
    delete tmpEvent;
  }
  fFirstPartOfBline = false;   

  G4PrimaryVertex* primary_vertex = 
    new G4PrimaryVertex(fBlineStartPosition, fT0);

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
