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
// $Id: G4SDProcessFilter.cc,v 1.1 2007-08-14 21:23:52 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VSensitiveDetector
#include "G4SDProcessFilter.hh"
#include "G4Step.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//  This class filters steps by partilce definition.
//
// Created: 2007-03-23  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

G4SDProcessFilter::G4SDProcessFilter(G4String name)
  :G4VSDFilter(name)
{
  theProcessDef.clear();
  theParticleDef.clear();
}

G4SDProcessFilter::G4SDProcessFilter(G4String name,
				     const G4VProcess* process, 
				     const G4String& particleName)
  :G4VSDFilter(name)
{
  theProcessDef.clear();
  theParticleDef.clear();
  G4ParticleDefinition* pd = 
      G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(!pd || !process)
  {
    G4String msg = "Process <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4SDProcessFilter::G4SDProcessFilter()",
                "DetUtil0000",FatalException,msg);
  }
  theParticleDef.push_back(pd);
  theProcessDef.push_back(process);
}

G4SDProcessFilter::~G4SDProcessFilter()
{ 
  theParticleDef.clear();
  theProcessDef.clear();
}

G4bool G4SDProcessFilter::Accept(const G4Step* aStep) const
{
  for ( size_t i = 0; i < theProcessDef.size(); i++){
    if ( theProcessDef[i] == aStep->GetPreStepPoint()->GetProcessDefinedStep() )
	if ( theParticleDef[i]   == aStep->GetTrack()->GetDefinition() ) {
	    return TRUE;
	}
  }
  return FALSE;
}

void G4SDProcessFilter::add(const G4VProcess* process, const G4String& particleName)
{
  G4ParticleDefinition* pd = 
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(!pd || !process )
  {
    G4String msg = "Process <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4SDProcessFilter::add()",
                "DetUtil0000",FatalException,msg);
  }

  for ( size_t i = 0; i < theProcessDef.size(); i++){
    if ( theProcessDef[i] == process && theParticleDef[i] == pd ) return;
  }
  theProcessDef.push_back(process);
  theParticleDef.push_back(pd);
}

void G4SDProcessFilter::show(){
  G4cout << "----G4SDProcessFilter particle list------"<<G4endl;
  for ( size_t i = 0; i < theProcessDef.size(); i++){
      G4cout << " Process : " << theProcessDef[i] 
	     << " Particle : "<< theParticleDef[i]->GetParticleName()<< G4endl;
  }
  G4cout << "-------------------------------------------"<<G4endl;
}



