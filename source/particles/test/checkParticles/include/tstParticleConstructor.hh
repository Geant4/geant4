// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tstParticleConstructor.hh,v 1.2 1999-10-03 09:13:18 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tstParticleConstructor_h
#define tstParticleConstructor_h 1

#include "globals.hh"

class tstParticleConstructor
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
