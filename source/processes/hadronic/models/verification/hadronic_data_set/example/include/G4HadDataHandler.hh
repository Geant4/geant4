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
// R&D: Vladimir.Grichine@cern.ch

#ifndef G4HadDataHandler_HH
#define G4HadDataHandler_HH 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"
#include "G4AugerTransition.hh"

class G4DataVector;

class G4HadDataHandler
{
public:

  G4HadDataHandler();

  ~G4HadDataHandler();

private:

  typedef G4std::map<G4int,G4std::
                 vector<G4AugerTransition>,G4std::less<G4int> > trans_Table;

  trans_Table augerTransitionTable;

  G4std::vector<G4int> nXSC;
  G4std::vector<G4int> numberOfSecondaries;
  
};

#endif





