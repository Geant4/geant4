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

#ifndef XUnitCell_h
#define XUnitCell_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"

#define MAXBASENUMBER 32

using namespace std;

class XUnitCell{

private:
    G4int fNumberOfBases;
    XLogicalBase* fBase[MAXBASENUMBER];
    
    G4ThreeVector fSize;
    G4ThreeVector fAngle;
    
public:
    //Retrieval methods
    G4ThreeVector GetSize();
    G4ThreeVector GetAngle();
    XLogicalBase* GetBase(G4int);
    
    
    //Set methods
    void SetSize(G4ThreeVector);
    void SetAngle(G4ThreeVector);
    void SetBase(G4int,XLogicalBase*);
    void AddBase(XLogicalBase*);

    //Calculation methods
    G4double ComputeVolume();

    G4double ComputeMillerOverSizeSquared(G4int,G4int,G4int);
    G4double ComputeMillerPerSizeSquared(G4int,G4int,G4int);

    G4double ComputeReciprocalVectorSquared(G4int,G4int,G4int);
    G4double ComputeReciprocalVector(G4int,G4int,G4int);

    G4double ComputeDirectVectorSquared(G4int,G4int,G4int);
    G4double ComputeDirectVector(G4int,G4int,G4int);
    
    G4double ComputeDirectPeriodSquared(G4int,G4int,G4int);
    G4double ComputeDirectPeriod(G4int,G4int,G4int);

    
    G4double ComputeAtomVolumeDensity();
    G4complex ComputeStructureFactor(G4int,G4int,G4int);
    //Kittel - chapter 2 Eq. (46)

    //Check method
    G4bool IsOrthogonal();
    G4bool IsCubic();
    
    //Contructors
    XUnitCell();
    ~XUnitCell();
};

#endif
