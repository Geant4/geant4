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
/// $Id: G4DigiFilterFactories.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//
// Digi filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#include "G4ModelCommandsT.hh"
#include "G4DigiFilterFactories.hh"
#include "G4AttributeFilterT.hh"

// Attribute filter
G4DigiAttributeFilterFactory::G4DigiAttributeFilterFactory()
  :G4VModelFactory< G4VFilter<G4VDigi> >("attributeFilter") 
{}

G4DigiAttributeFilterFactory::~G4DigiAttributeFilterFactory() {}

G4DigiAttributeFilterFactory::ModelAndMessengers
G4DigiAttributeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  typedef G4AttributeFilterT<G4VDigi> G4DigiAttributeFilter;
  // Create model
  G4DigiAttributeFilter* model = new G4DigiAttributeFilter(name);
  
  // Create associated messengers
  Messengers messengers;
  
  messengers.push_back(new G4ModelCmdSetString<G4DigiAttributeFilter>(model, placement, "setAttribute"));
  messengers.push_back(new G4ModelCmdInvert<G4DigiAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4DigiAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4DigiAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4DigiAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdAddInterval<G4DigiAttributeFilter>(model, placement, "addInterval"));
  messengers.push_back(new G4ModelCmdAddValue<G4DigiAttributeFilter>(model, placement, "addValue"));
  
  return ModelAndMessengers(model, messengers);
}

