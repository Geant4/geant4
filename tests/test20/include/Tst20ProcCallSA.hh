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

#ifndef Tst20ProcCallSA_H
#define Tst20ProcCallSA_H 1
//
// stepping Action:
//  counts number of times a process is called, per material, per particle type.
// 
//  veronique.lefebure@cern.ch  20.01.2000
//

//#define histo

#include "globals.hh"
#include <g4std/map>

class G4Step;
class Tst20ProcCallSA{
public:
  Tst20ProcCallSA();
  ~Tst20ProcCallSA();
  
  void execute(const G4Step*);
  
private:
  void print();
  
private:
  typedef G4std::map<G4String, G4int, G4std::less<G4String> > intMap;
  typedef G4std::map<G4String, G4int, G4std::less<G4String> >::iterator intMapIter;
  intMap calls;  
  
};

#endif

