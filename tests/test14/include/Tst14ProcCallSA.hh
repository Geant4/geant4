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

#ifndef Tst14ProcCallSA_H
#define Tst14ProcCallSA_H 1
//
// stepping Action:
//  counts number of times a process is called, per material, per particle type.
// 
//  veronique.lefebure@cern.ch  20.01.2000
//

//#define histo

#include "globals.hh"
#include <map>

#ifdef histo
class HepTupleManager;
class HepHistogram;
#endif
class G4Step;
class Tst14ProcCallSA{
  public:
    Tst14ProcCallSA();
     ~Tst14ProcCallSA();

     void execute(const G4Step*);

  private:
     void print();
          
  private:
     typedef std::map<G4String, G4int, std::less<G4String> > intMap;
     typedef std::map<G4String, G4int, std::less<G4String> >::iterator intMapIter;
     intMap calls;  
    
     
#ifdef histo
private:

  typedef std::map<G4String, HepHistogram*, std::less<G4String> > histoMap;
  typedef std::map<G4String, HepHistogram*, std::less<G4String> >::iterator histoMapIter;
  HepTupleManager* hbookManager;
  histoMap       hist,start;
#endif
};

#endif

