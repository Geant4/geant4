// $Id: ParticleGunAction.cc,v 1.1 2006-02-27 09:52:54 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   ParticleGunAction.cc
//
//                                         2005 Q
// ====================================================================
#include "ParticleGunAction.hh"
#include "G4ParticleGun.hh"
 
// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////
ParticleGunAction::ParticleGunAction()
//////////////////////////////////////
{
  particleGun= new G4ParticleGun;
}

///////////////////////////////////////
ParticleGunAction::~ParticleGunAction()
///////////////////////////////////////
{
  delete particleGun;
}

///////////////////////////////////////////////////////////
void ParticleGunAction::GeneratePrimaries(G4Event* anEvent)
///////////////////////////////////////////////////////////
{
  particleGun-> GeneratePrimaryVertex(anEvent);
}

