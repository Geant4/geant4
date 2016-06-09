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
// $Id: G4PVDivisionFactory.hh,v 1.1 2004/05/13 14:56:30 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//
// class G4PVDivisionFactory
//
// Class description:
//
// Implementation of the interfaces for creating volumes divisions
// (defined in G4VPVDivisionFactory) for G4PVDivision type.   

// Author: Ivana Hrivnacova, 4.5.2004  (Ivana.Hrivnacova@cern.ch)
// ------------------------------------------------------------------------
#ifndef G4PVDIVISION_FACTORY_HH
#define G4PVDIVISION_FACTORY_HH

#include "geomdefs.hh"
#include "G4VPVDivisionFactory.hh"

class G4LogicalVolume;

class G4PVDivisionFactory : public G4VPVDivisionFactory
{
  public:  // with description

    virtual ~G4PVDivisionFactory();
    
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double width,
                             const G4double offset );
      // Create division - with number of divisions and width.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double offset );
      // Create division - with number of divisions.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4double width,
                             const G4double offset );
      // Create division - with width.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const G4VPVParameterisation* param);
      // Create division - with parameterisation.

    virtual G4bool IsPVDivision(const G4VPhysicalVolume* pv) const; 
      // Returns true if pv is division.
    
    static G4PVDivisionFactory* GetInstance(); 
      // Create the unique instance of the singleton.

  protected:

    G4PVDivisionFactory();
     
};

#endif
