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
// File name:     RadmonPhysicsProductionCuts.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsProductionCuts.cc,v 1.3 2006/06/29 16:19:52 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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


