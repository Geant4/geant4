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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AtomicShell.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  
//  16 Sept 2001 Modofied according to a design iteration in the 
//              LowEnergy category
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, a data container
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4AtomicShell_h 
#define G4AtomicShell_h 1
#include "globals.hh"

class G4AtomicShell {

public:

  // The data and the methods of this class are relative to
  // a given shell
 
  G4AtomicShell(G4int,G4double);
 
  ~G4AtomicShell();

  // Returns the binding energy of the shell
  G4double BindingEnergy() const; 

  // Returns the id of the shell
  G4int ShellId() const;

private:
  
  G4int identifier;
  G4double bindingEnergy;

};

#endif


