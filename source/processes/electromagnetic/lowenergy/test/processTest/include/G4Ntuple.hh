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
// $Id: G4Ntuple.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: MGP
//
// History:
// -----------
//  

//
// -------------------------------------------------------------------

// Class description:

// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4NTUPLE_HH
#define G4NTUPLE_HH 1

#include "globals.hh"
#include "g4std/vector"

// For old CLHEP stuff
#include "G4DataVector.hh"

// For NtupleTag from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
using namespace Lizard;

class G4Ntuple {


public: 

  G4Ntuple();
 
  ~G4Ntuple();

  void Book(const G4String& ntupleName);

  G4bool AddAndBind(const G4String& id, float q);

  // First one for obsolete CLHEP stuff, second one for new code
  void AddRow(const G4DataVector& row);
  //  void AddRow();

  const G4String& GetName() { return name; }

  G4int NumberOfAttributes();

private:
  
  G4String name;
  
  // Used for old CLHEP ntuples
  G4std::vector<G4String> attributes;

  NTuple* ntuple;

};

#endif
