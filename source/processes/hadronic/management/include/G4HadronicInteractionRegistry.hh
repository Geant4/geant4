//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4HadronicInteractionRegistry_h
#define G4HadronicInteractionRegistry_h 1

#include "g4std/vector"
#include "globals.hh"
class G4HadronicInteraction;

class G4HadronicInteractionRegistry
{
  public:
  
  ~G4HadronicInteractionRegistry();
  
  static void RegisterMe(G4HadronicInteraction * aModel);
  static void RemoveMe(G4HadronicInteraction * aModel);
  
  protected:

  G4HadronicInteractionRegistry(G4String aString) 
  { G4Exception("G4HadronicInteractionRegistry meant as a singleton; please do not inherit");}

  private:

  //  !!!  can not use "copy constructor" nor "default constructor" !!!!
       G4HadronicInteractionRegistry(const G4HadronicInteractionRegistry &right) 
       { nModels = right.nModels; }
       G4HadronicInteractionRegistry() {nModels = 0;}

  //  !!!  Assignment operation is forbidden !!!
      const G4HadronicInteractionRegistry & operator=(const G4HadronicInteractionRegistry &right) 
      { return *this;}

  void AddModel(G4HadronicInteraction * aModel);
  
  G4int nModels;
  G4std::vector <G4HadronicInteraction *> allModels;
  static G4HadronicInteractionRegistry theRegistry;

};

#endif
