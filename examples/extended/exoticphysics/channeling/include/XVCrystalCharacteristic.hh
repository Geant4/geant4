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

#ifndef XVCrystalCharacteristic_h
#define XVCrystalCharacteristic_h

#include "XLatticeManager3.hh"
#include "G4Track.hh"
#include "G4PhysicsVector.hh"

class XVCrystalCharacteristic {

private:
    XLatticeManager3* fLatticeManager;

protected:
    G4double fMaximum;
    G4double fMinimum;
    XPhysicalLattice *fPhysicalLattice;
    G4PhysicsVector* fVectorEC;
    
public:
    //retrieval functions
    XPhysicalLattice* GetXPhysicalLattice(G4VPhysicalVolume*);
    XUnitCell* GetXUnitCell(G4VPhysicalVolume*);
    XLogicalLattice* GetLogicalLattice(G4VPhysicalVolume*);
    void InitializePhysicalLattice(XPhysicalLattice*);
    
    //virtual function to compute value starting from the point in the xtal
    //reference frame and the physical volume of the xtal
    G4ThreeVector GetEC(G4ThreeVector,XPhysicalLattice*);
    virtual G4ThreeVector ComputeEC(G4ThreeVector,XPhysicalLattice*) = 0;
    virtual G4ThreeVector ComputeECFromVector(G4ThreeVector) = 0;
    virtual G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,
                                                    XPhysicalLattice*);
    
    virtual G4double ComputeTFScreeningRadius(XPhysicalLattice*);

    virtual G4double GetMaximum(XPhysicalLattice*);
    virtual G4double GetMinimum(XPhysicalLattice*);

    virtual G4double ComputeMaximum(XPhysicalLattice*);
    virtual G4double ComputeMinimum(XPhysicalLattice*);
    
    virtual void PrintOnFile(const G4String&,XPhysicalLattice*,
                             G4double = 1) = 0;
    virtual void ReadFromFile(const G4String&,XPhysicalLattice*,
                              G4double = 1) = 0;
    virtual void ReadFromECHARM(const G4String&,
                                G4double = 1) = 0;

    G4bool IsInitialized(XPhysicalLattice*);
    virtual void InitializeVector() = 0;
    
    //Contructors
    XVCrystalCharacteristic();
    ~XVCrystalCharacteristic();
};

#endif
