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

#ifndef XLogicalAtomicLattice_h
#define XLogicalAtomicLattice_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"

#ifndef MAXLATTICEATOMS
#define MAXLATTICEATOMS 64
#endif

using namespace std;

class XLogicalAtomicLattice{

private:
    // position of the atoms are saved in unit cell system, i.e MIN 0. & MAX 1.
    G4ThreeVector fLatticeAtomPosition[MAXLATTICEATOMS];
    G4int fLatticeAtomNumber;
    
public:    
    void InitializeXLogicalAtomicLattice();

    // Get methods
    G4ThreeVector GetAtomPosition(G4int i);
    G4int GetLatticeNumberOfAtoms();
    
    // Set methods
    void AddAtom(G4ThreeVector);
    void DeleteAtom(G4ThreeVector);
    

    // Calculation methods
    // ints == Miller indexes
    G4complex ComputeGeometricalStructureFactorSingleKind(G4int,G4int,G4int);

    // Definition methods
    XLogicalAtomicLattice();
    ~XLogicalAtomicLattice();
};

#endif
