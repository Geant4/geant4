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
//
// Hits filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#include "G4ModelCommandsT.hh"
#include "G4HitFilterFactories.hh"
#include "G4AttributeFilterT.hh"

// Attribute filter
G4HitAttributeFilterFactory::G4HitAttributeFilterFactory()
  :G4VModelFactory< G4VFilter<G4VHit> >("attributeFilter") 
{}

G4HitAttributeFilterFactory::~G4HitAttributeFilterFactory() {}

G4HitAttributeFilterFactory::ModelAndMessengers
G4HitAttributeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  typedef G4AttributeFilterT<G4VHit> G4HitAttributeFilter;
  // Create model
  G4HitAttributeFilter* model = new G4HitAttributeFilter(name);
  
  // Create associated messengers
  Messengers messengers;
  
  messengers.push_back(new G4ModelCmdSetString<G4HitAttributeFilter>(model, placement, "setAttribute"));
  messengers.push_back(new G4ModelCmdInvert<G4HitAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4HitAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4HitAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4HitAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdAddInterval<G4HitAttributeFilter>(model, placement, "addInterval"));
  messengers.push_back(new G4ModelCmdAddValue<G4HitAttributeFilter>(model, placement, "addValue"));
  
  return ModelAndMessengers(model, messengers);
}

