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
#include "G4HadSecondary.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicException.hh"

G4HadSecondary::G4HadSecondary(G4DynamicParticle * aT, G4double aWeight, 
			       G4int mod) :
  theP(aT), theWeight(aWeight), theTime(-1), theCreatorModel(mod),
  theParentResonanceDef(nullptr), theParentResonanceID(0)
{
  if(aT->GetKineticEnergy()<0)
  {
    throw G4HadronicException(__FILE__, __LINE__, 
    "ATTEMPTING TO CREATE A SECONDARY WITH NEGATIVE KINETIC ENERGY.");
  }
}

G4HadSecondary::~G4HadSecondary()
{}
