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
// File name:     RadmonMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessenger.cc,v 1.3 2006/06/29 16:14:47 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Include files
#include "RadmonMessenger.hh"
#include "RadmonTokenizer.hh"
#include "G4UnitsTable.hh"
#include "G4String.hh"

#include <fstream>

G4bool                                          RadmonMessenger :: ProcessArguments(const G4String & rawArguments, G4int nArgs, G4String * arguments)
{
 RadmonTokenizer args(rawArguments);
 
 for (G4int i(0); i<nArgs; i++)
 {
  if (args.eos())
  {
   G4cout << "RadmonMessenger::ProcessArguments: " << (nArgs-i) <<  " arguments missing." << G4endl;
   return false;
  }

  arguments[i]=args();
 }

 if (!args.eos())
 {
  G4cout << "RadmonMessenger::ProcessArguments: Unexpected arguments after \"" << arguments[nArgs-1] << "\"." << G4endl;
  return false;
 }
 
 return true;
}



G4double                                        RadmonMessenger :: GetUnit(const G4String & unitStr, const char * category)
{
 G4double dbl(G4UnitDefinition::GetValueOf(unitStr));

 if (dbl<=0.)
 {
  G4cout << "RadmonMessenger::GetUnit(): Unknown unit of measure \"" << unitStr << "\"." << G4endl;
  return -1.;
 }
 
 if (G4UnitDefinition::GetCategory(unitStr)!=category)
 {
  G4cout << "RadmonMessenger::GetUnit(): Unit of measure \"" << unitStr << "\" is not a " << category << '.' << G4endl;
  return -1.;
 }
 
 return dbl;
}





                                                RadmonMessenger :: RadmonMessenger(const char * path, const char * guidance)
{
 directory=new G4UIdirectory(path);
 
 if (directory==0)
 {
  G4String text("RadmonMaterialsMessenger::RadmonMaterialsMessenger: \"");
  text+=path;
  text+="\" directory not allocated.";
  
  G4Exception(text);
  
  return;
 }
 
 if (guidance)
  directory->SetGuidance(guidance); 
}





std::istream *                                  RadmonMessenger :: OpenForInput(const G4String & fileName) const
{
 std::ifstream * in(new std::ifstream(fileName.data()));
 
 if (!in)
 {
  G4cout << "RadmonMessenger::OpenForInput(): Cannot allocate class std::ifstream." << G4endl;
  return 0;
 }

 if (!in->good())
 {
  G4cout << "RadmonMessenger::OpenForInput(): Cannot access to file \"" << fileName << "\" for read." << G4endl;
  delete in;
  return 0;
 }
 
 return in;
}



std::ostream *                                  RadmonMessenger :: OpenForOutput(const G4String & fileName) const
{
 std::ofstream * out(new std::ofstream(fileName.data()));
 
 if (!out)
 {
  G4cout << "RadmonMessenger::OpenForOutput(): Cannot allocate class std::ofstream." << G4endl;
  return 0;
 }

 if (!out->good())
 {
  G4cout << "RadmonMessenger::OpenForOutput(): Cannot access to file \"" << fileName << "\" for write." << G4endl;
  delete out;
  return 0;
 }
 
 return out;
}
