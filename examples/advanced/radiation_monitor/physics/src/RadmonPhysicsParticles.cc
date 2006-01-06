//
// File name:     RadmonPhysicsParticles.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsParticles.cc,v 1.1 2006-01-06 12:52:32 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsParticles.hh"

#include "G4ProcessManager.hh"

#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsParticles :: New(void) const
{
 return new RadmonPhysicsParticles;
}



void                                            RadmonPhysicsParticles :: ConstructParticle(void)
{
  G4LeptonConstructor lepton;
  lepton.ConstructParticle();
 
  G4BosonConstructor boson;
  boson.ConstructParticle();

  G4MesonConstructor meson;
  meson.ConstructParticle();

  G4BaryonConstructor baryon;
  baryon.ConstructParticle();

  G4ShortLivedConstructor shortLived;
  shortLived.ConstructParticle();

  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}



void                                            RadmonPhysicsParticles :: ConstructProcess(void)
{
}



void                                            RadmonPhysicsParticles :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsParticles :: Provides(void) const
{
 return infoList;
}
