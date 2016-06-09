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
// $Id$
//
// 23-Jan-2009 V.Ivanchenko make the class to be a singleton
// 17-Aug-2012 V.Ivanchenko added hadronic model factories

// Class Description
// This is a singleton class to store all hadronic interactions
// Class Description - End

#ifndef G4HadronicInteractionRegistry_h
#define G4HadronicInteractionRegistry_h 1

#include <vector>
#include "globals.hh"

class G4HadronicInteraction;

class G4HadronicInteractionRegistry
{
public:

  static G4HadronicInteractionRegistry* Instance();
  // access 
  
  ~G4HadronicInteractionRegistry();

  void Clean();
  // delete models
  
  void RegisterMe(G4HadronicInteraction * aModel);
  // register new model

  void RemoveMe(G4HadronicInteraction * aModel);
  // deregister model

  G4HadronicInteraction* FindModel(const G4String& name);
  // find existing hadronic interaction by name

private:

  G4HadronicInteractionRegistry();

  static G4HadronicInteractionRegistry* theInstance;
  
  std::vector <G4HadronicInteraction*>  allModels;

};

#endif
