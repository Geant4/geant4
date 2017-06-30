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
// $Id:$
// GEANT4 tag $Name:$
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------

#include "G4VPhysicsConstructor.hh"
#include "G4PhysicsBuilderInterface.hh"
#include <algorithm>
// This field helps to use the class G4VPCManager
//
G4VPCManager G4VPhysicsConstructor::subInstanceManager;

void G4VPCData::initialize()
{
    _aParticleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
    _builders = new PhysicsBuilders_V;
}


G4VPhysicsConstructor::G4VPhysicsConstructor(const G4String& name)
   : verboseLevel(0), namePhysics(name), typePhysics(0)
{
  g4vpcInstanceID = subInstanceManager.CreateSubInstance();
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  //aParticleIterator = theParticleTable->GetIterator();

  // PhysicsListHelper
  //aPLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
}

G4VPhysicsConstructor::G4VPhysicsConstructor(const G4String& name, G4int type)
    : verboseLevel(0), namePhysics(name), typePhysics(type)
{
    g4vpcInstanceID = subInstanceManager.CreateSubInstance();
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  //aParticleIterator = theParticleTable->GetIterator();

  if (type<0) typePhysics = 0;

  // PhysicsListHelper
  //aPLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
}

G4VPhysicsConstructor::~G4VPhysicsConstructor()
{
  //Master/Sequential needs to cleanup too
  G4VPhysicsConstructor::TerminateWorker();
}

G4ParticleTable::G4PTblDicIterator* G4VPhysicsConstructor::GetParticleIterator() const
{
	return (subInstanceManager.offset[g4vpcInstanceID])._aParticleIterator;
}

G4VPhysicsConstructor::PhysicsBuilder_V G4VPhysicsConstructor::GetBuilders() const
{
  const auto& tls = *((subInstanceManager.offset[g4vpcInstanceID])._builders);
  PhysicsBuilder_V copy(tls.size());
  int i = 0;
  for ( const auto& el : tls ) { copy[i++] = el; }
  return copy;
}

void G4VPhysicsConstructor::AddBuilder(G4PhysicsBuilderInterface* bld)
{
  (subInstanceManager.offset[g4vpcInstanceID])._builders->push_back(bld);
}

void G4VPhysicsConstructor::TerminateWorker()
{
  if ( subInstanceManager.offset[g4vpcInstanceID]._builders != nullptr ) {
  std::for_each( subInstanceManager.offset[g4vpcInstanceID]._builders->begin() ,
		 subInstanceManager.offset[g4vpcInstanceID]._builders->end() ,
		 [](PhysicsBuilder_V::value_type bld) { delete bld;});
  subInstanceManager.offset[g4vpcInstanceID]._builders->clear();
  }

}
