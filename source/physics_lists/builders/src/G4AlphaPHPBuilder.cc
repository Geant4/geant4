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
#include "G4AlphaPHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4SystemOfUnits.hh"

G4AlphaPHPBuilder::
G4AlphaPHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
  theParticlePHPModel = 0;
}

void G4AlphaPHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4AlphaPHPBuilder::Build() not adding elastic" << G4endl;
}

void G4AlphaPHPBuilder::
Build(G4AlphaInelasticProcess * aP)
{
  G4cout << " G4AlphaPHPBuilder " << G4endl;
  G4ParticleHPInelasticData* theAlphaHPInelasticData=new G4ParticleHPInelasticData(G4Alpha::Alpha());
  theAlphaHPInelasticData->SetMinKinEnergy(theMin);
  theAlphaHPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theAlphaHPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4Alpha::Alpha(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);

}

