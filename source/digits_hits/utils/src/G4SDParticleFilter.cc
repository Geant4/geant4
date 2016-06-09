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
// $Id: G4SDParticleFilter.cc,v 1.2 2005/11/22 21:41:55 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4VSensitiveDetector
#include "G4SDParticleFilter.hh"
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
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

G4SDParticleFilter::G4SDParticleFilter(G4String name)
  :G4VSDFilter(name)
{
  thePdef.clear();
}

G4SDParticleFilter::G4SDParticleFilter(G4String name,
				       const G4String& particleName)
  :G4VSDFilter(name)
{
  thePdef.clear();
  G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(!pd)
  {
    G4String msg = "Particle <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4SDParticleFilter::G4SDParticleFilter()",
                "DetUtil0000",FatalException,msg);
  }
  thePdef.push_back(pd);
}

G4SDParticleFilter::G4SDParticleFilter(G4String name, 
			       const std::vector<G4String>& particleNames)
  :G4VSDFilter(name)
{
  thePdef.clear();
  for ( size_t i = 0; i < particleNames.size(); i++){
   G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle(particleNames[i]);
   if(!pd)
   {
     G4String msg = "Particle <";
     msg += particleNames[i];
     msg += "> not found.";
     G4Exception("G4SDParticleFilter::G4SDParticleFilter()",
                "DetUtil0000",FatalException,msg);
   }
   thePdef.push_back(pd);
  }
}

G4SDParticleFilter::G4SDParticleFilter(G4String name, 
		       const std::vector<G4ParticleDefinition*>& particleDef)
  :G4VSDFilter(name), thePdef(particleDef)
{
  for ( size_t i = 0; i < particleDef.size(); i++){
    if(!particleDef[i]) G4Exception("G4SDParticleFilter::G4SDParticleFilter()","DetUtil0001",FatalException,
       "NULL pointer is found in the given particleDef vector.");
  }
}

G4SDParticleFilter::~G4SDParticleFilter()
{ 
  thePdef.clear();
}

G4bool G4SDParticleFilter::Accept(const G4Step* aStep) const
{
  for ( size_t i = 0; i < thePdef.size(); i++){
    if ( thePdef[i] == aStep->GetTrack()->GetDefinition() ) return TRUE;
  }
  return FALSE;
}

void G4SDParticleFilter::add(const G4String& particleName)
{
  G4ParticleDefinition* pd = 
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(!pd)
  {
     G4String msg = "Particle <";
     msg += particleName;
     msg += "> not found.";
     G4Exception("G4SDParticleFilter::add()",
                "DetUtil0002",FatalException,msg);
  }
  for ( size_t i = 0; i < thePdef.size(); i++){
    if ( thePdef[i] == pd ) return;
  }
  thePdef.push_back(pd);
}

void G4SDParticleFilter::show(){
  G4cout << "----G4SDParticleFileter particle list------"<<G4endl;
  for ( size_t i = 0; i < thePdef.size(); i++){
    G4cout << thePdef[i]->GetParticleName() << G4endl;
  }
  G4cout << "-------------------------------------------"<<G4endl;
}



