// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst14ProcCallSA_H
#define Tst14ProcCallSA_H 1
//
// stepping Action:
//  counts number of times a process is called, per material, per particle type.
// 
//  veronique.lefebure@cern.ch  20.01.2000
//
#include "globals.hh"
#include <g4std/map>
class G4Step;
class Tst14ProcCallSA{
  public:
    Tst14ProcCallSA();
     ~Tst14ProcCallSA();

     void execute(const G4Step*);

  private:
     void print();
          
  private:
     typedef map<G4String, G4int, less<G4String> >::iterator iter;
     G4std::map<G4String, G4int, less<G4String> > particles;   
     G4std::map<G4String, G4int, less<G4String> > processes;  
     G4std::map<G4String, G4int, less<G4String> > materials;  
     G4std::map<G4int, G4int, less<G4int> > calls;  
     
     G4int  nparticles, nprocesses, nmaterials;
};

#endif

