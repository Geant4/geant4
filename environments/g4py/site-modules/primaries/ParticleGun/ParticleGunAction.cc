// $Id: ParticleGunAction.cc,v 1.2 2006-04-25 10:29:52 kmura Exp $
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

