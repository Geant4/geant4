// $Id: QPrimaryGeneratorAction.cc,v 1.1 2006-05-11 03:00:08 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QPrimaryGeneratorAction.cc
//
//                                         2005 Q
// ====================================================================
#include "QPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
 
// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////////////
QPrimaryGeneratorAction::QPrimaryGeneratorAction()
//////////////////////////////////////////////////
{
  particleGun= new G4ParticleGun;
}

///////////////////////////////////////////////////
QPrimaryGeneratorAction::~QPrimaryGeneratorAction()
///////////////////////////////////////////////////
{
  delete particleGun;
}

/////////////////////////////////////////////////////////////////
void QPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
/////////////////////////////////////////////////////////////////
{
  particleGun-> GeneratePrimaryVertex(anEvent);
}

