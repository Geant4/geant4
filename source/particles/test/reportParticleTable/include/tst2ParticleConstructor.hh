// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2ParticleConstructor.hh,v 1.3 1999-12-15 14:51:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2ParticleConstructor_h
#define tst2ParticleConstructor_h 1

#include "globals.hh"

class tst2ParticleConstructor
{
  public:
	void ConstructParticle();

  protected:
	void ConstructAllBosons();
	void ConstructAllLeptons();
	void ConstructAllMesons();
	void ConstructAllBaryons();
	void ConstructAllIons();
	void ConstructAllShortLiveds();
};
#endif
