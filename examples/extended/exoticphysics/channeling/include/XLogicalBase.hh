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

#ifndef XLogicalBase_h
#define XLogicalBase_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"
#include "XLogicalAtomicLattice.hh"
#include "G4Element.hh"

using namespace std;

class XLogicalBase{

private:
    XLogicalAtomicLattice* fLattice;
    G4Element* fElement;
    
public:
    // Retrieval methods
    XLogicalAtomicLattice* GetLattice();
    G4Element* GetElement();
    
    // Set methods
    void SetLattice(XLogicalAtomicLattice*);
    void SetElement(G4Element*);
    
    // Calculation methods
    // ints == Miller indexes
    virtual G4double ComputeAtomicFormFactor();
    //Kittel - chapter 2 Eq. (42) for single atomic kind
    G4complex ComputeStructureFactorSingleAtomicKind(G4int,G4int,G4int);
    //Kittel - chapter 2 Eq. (46) for single atomic kind
  
    XLogicalBase(G4Element*,XLogicalAtomicLattice*);
    XLogicalBase();
    ~XLogicalBase();
};

#endif
