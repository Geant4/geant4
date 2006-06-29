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
#ifndef G4HadronicInteractionRegistry_h
#define G4HadronicInteractionRegistry_h 1

#include <vector>
#include "globals.hh"
class G4HadronicInteraction;

class G4HadronicInteractionRegistry
{
  public:
  
  ~G4HadronicInteractionRegistry();
  
  static void RegisterMe(G4HadronicInteraction * aModel);
  static void RemoveMe(G4HadronicInteraction * aModel);
    
  protected:

  G4HadronicInteractionRegistry(G4String );
  static G4HadronicInteractionRegistry & theRegistry();

  private:

  //  !!!  can not use "copy constructor" nor "default constructor" !!!!
       G4HadronicInteractionRegistry(const G4HadronicInteractionRegistry &right) 
       { nModels = right.nModels; }
       G4HadronicInteractionRegistry() {nModels = 0;}

  //  !!!  Assignment operation is forbidden !!!
      const G4HadronicInteractionRegistry & operator=(const G4HadronicInteractionRegistry &) 
      { return *this;}

  void AddModel(G4HadronicInteraction * aModel);
  
  G4int nModels;
  std::vector <G4HadronicInteraction *> allModels;

};

#endif
