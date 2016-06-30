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
// $Id: G4AtomicShell.hh 96626 2016-04-27 08:36:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  
//  16 Sept 2001 Modified according to a design iteration in the 
//              LowEnergy category
//  24 July 2009 Moved to utils subdirectory and make inline methods (VI) 
//
// -------------------------------------------------------------------

// Class description:
// A container of atomic shell data   

// -------------------------------------------------------------------


#ifndef G4AtomicShell_h 
#define G4AtomicShell_h 1
#include "globals.hh"

class G4AtomicShell {

public:

  // The data and the methods of this class are relative to
  // a given shell
 
  explicit G4AtomicShell(G4int id, G4double e):identifier(id),bindingEnergy(e){};
 
  // Returns the binding energy of the shell
  inline G4double BindingEnergy() const; 

  // Returns the id of the shell
  inline G4int ShellId() const;

private:
  
  G4AtomicShell & operator=(const  G4AtomicShell &right) = delete;
  G4AtomicShell(const  G4AtomicShell&) = delete;

  G4int identifier;
  G4double bindingEnergy;

};

inline G4double G4AtomicShell::BindingEnergy() const {

  return bindingEnergy;
}

inline G4int G4AtomicShell::ShellId() const{

  return identifier;
}

#endif


