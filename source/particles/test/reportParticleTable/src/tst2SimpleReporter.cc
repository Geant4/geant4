// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2SimpleReporter.cc,v 1.2 1999-12-15 14:51:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "tst2SimpleReporter.hh"
#include "G4ios.hh"
#include "globals.hh"

 tst2SimpleReporter::tst2SimpleReporter():tst2VParticleReporter()
{
 
}

 tst2SimpleReporter::~tst2SimpleReporter()
{
}    

 void tst2SimpleReporter::Print(const tst2ParticleContainer& container, 
                               const G4String& option)
{
  pList = &container;
  G4cout << " Encoding    " << "name " << G4endl;
  for (G4int i=0; i< entries(); i++){
	G4cout << GetEncoding(i) << " :   " << GetParticle(i)->GetParticleName() << G4endl;
  }
}    
