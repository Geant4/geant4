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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2ParticleContainer.cc,v 1.3 2001-07-11 10:02:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "tst2ParticleContainer.hh"

 tst2ParticleContainer::tst2ParticleContainer()
{
  pVector =  new  tst2ParticleVector;
  pTable = G4ParticleTable::GetParticleTable();
}

 tst2ParticleContainer::~tst2ParticleContainer()
{
  pVector->clearAndDestroy();
  delete pVector;
}    

  tst2ParticleContainer::tst2ParticleContainer(const tst2ParticleContainer & right)
{
	pVector = new tst2ParticleVector;
	for (G4int i=0; i< right.entries(); i++) {
	  tst2ContainerElement* element =  new tst2ContainerElement( right.GetParticle(i), right.GetEncoding(i) );
      pVector->insert(element);

    }
}

 tst2ParticleContainer & tst2ParticleContainer::operator=(const tst2ParticleContainer & right)
{
  if (this != &right) { 
    pVector->clearAndDestroy();
    delete pVector;
	pVector = new tst2ParticleVector;
	for (G4int i=0; i< right.entries(); i++) {
	  tst2ContainerElement* element =  new tst2ContainerElement( right.GetParticle(i), right.GetEncoding(i) );
      pVector->insert(element);
    }
  }
  return *this;
}

G4int  tst2ParticleContainer::Insert( G4ParticleDefinition * particle){
  if (particle == 0) return 0;

  G4int encoding = particle->GetPDGEncoding();
  
  // anti_particle is not filled
  if (encoding <0) return 0;

  tst2ContainerElement* element =  new tst2ContainerElement(particle, encoding);
  pVector->insert(element);
 
  return 1;
}

 G4ParticleDefinition*  tst2ParticleContainer::GetAntiParticle(G4int index) const
{
  G4ParticleDefinition* p = GetParticle(index);

  if (p==0) return 0;
  if (p->GetAntiPDGEncoding() == 0) return 0;

  if (p->GetAntiPDGEncoding() == p->GetPDGEncoding()) return p;

  return pTable->FindParticle(p->GetAntiPDGEncoding());
}

