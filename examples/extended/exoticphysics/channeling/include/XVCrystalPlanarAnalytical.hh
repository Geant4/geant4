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

#ifndef XVCrystalPlanarAnalytical_h
#define XVCrystalPlanarAnalytical_h

#include "XVCrystalCharacteristic.hh"

class XVCrystalPlanarAnalytical:public XVCrystalCharacteristic {

private:
    G4int fNumberOfPlanes;

public:
    //set function
    void SetNumberOfPlanes(G4int);

    //retrieval function
    G4int GetNumberOfPlanes();

    virtual G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*) = 0;
    // G4double = position in the channel

    //virtual function of XVCrystalCharacteristic
    G4ThreeVector ComputeEC(G4ThreeVector,XPhysicalLattice*);
    G4ThreeVector ComputeECFromVector(G4ThreeVector);
    G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*);
    //G4double = position in the channel; G4double& = interplanar distance
    
    virtual G4double ComputeMaximum(XPhysicalLattice*);
    virtual G4double ComputeMinimum(XPhysicalLattice*);

    virtual void PrintOnFile(const G4String&,XPhysicalLattice*,G4double);
    virtual void ReadFromFile(const G4String&,XPhysicalLattice*,G4double = 1);
    virtual void ReadFromECHARM(const G4String&,G4double = 1);

    void InitializeVector();
    //Contructors
    XVCrystalPlanarAnalytical();
    ~XVCrystalPlanarAnalytical();
};

#endif
