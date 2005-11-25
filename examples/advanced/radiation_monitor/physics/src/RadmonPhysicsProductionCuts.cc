//
// File name:     RadmonPhysicsProductionCuts.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsProductionCuts.cc,v 1.1 2005-11-25 11:52:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsProductionCuts.hh"

#include "G4ProcessManager.hh"

#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4ParticleTable.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsProductionCuts :: New(void) const
{
 return new RadmonPhysicsProductionCuts;
}



void                                            RadmonPhysicsProductionCuts :: ConstructParticle(void)
{
}



void                                            RadmonPhysicsProductionCuts :: ConstructProcess(void)
{
}



void                                            RadmonPhysicsProductionCuts :: SetCuts(void)
{
 G4double cuts(GetAttributeAsMeasure("Cuts", "Length", -1.));

 if (cuts<0)
 {
  G4cout << "RadmonPhysicsProductionCuts::SetCuts: \"Cuts\" attribute is not defined." << G4endl;
  return;
 }

 SetProductionCut(cuts, "gamma");
 SetProductionCut(cuts, "e-");
 SetProductionCut(cuts, "e+");
}





const RadmonPhysicsInfoList &                   RadmonPhysicsProductionCuts :: Provides(void) const
{
 return infoList;
}





void                                            RadmonPhysicsProductionCuts :: SetProductionCut(G4double cut, const G4String& particleName) const
{
 G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
 G4ParticleDefinition * particleDefinition(particleTable->FindParticle(particleName));

 if (particleDefinition==0)
  return;
  
 if (particleDefinition->IsShortLived())
  return;

 G4ProductionCutsTable * cutsTable(G4ProductionCutsTable::GetProductionCutsTable());
 if (!cutsTable)
  return;
  
 G4ProductionCuts * productionCuts(cutsTable->GetDefaultProductionCuts());
 if (!productionCuts)
  return;

 productionCuts->SetProductionCut(cut, particleDefinition);
}


