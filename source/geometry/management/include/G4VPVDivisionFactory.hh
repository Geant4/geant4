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
// $Id: G4VPVDivisionFactory.hh,v 1.1 2004/05/13 14:53:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//
// class G4VPVDivisionFactory
//
// Class description:
//
// Abstract factory that defines interfaces for creating 
// volumes divisions. 

// Author: Ivana Hrivnacova, 4.5.2004  (Ivana.Hrivnacova@cern.ch)
// ------------------------------------------------------------------------
#ifndef G4VPVDIVISION_FACTORY_HH
#define G4VPVDIVISION_FACTORY_HH

#include "geomdefs.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VPVParameterisation;

class G4VPVDivisionFactory
{
  public:  // with description

    virtual ~G4VPVDivisionFactory();
    
    static  G4VPVDivisionFactory* Instance();
    
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double width,
                             const G4double offset ) = 0;
      // Create division - with number of divisions and width.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double offset ) = 0;
      // Create division - with number of divisions.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4double width,
                             const G4double offset ) = 0;
      // Create division - with width.

    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMotherLogical,
                             const G4VPVParameterisation* param) = 0;
      // Create division - with parameterisation.

    virtual G4bool IsPVDivision(const G4VPhysicalVolume* pv) const = 0;
      // Returns true if pv is division.

  protected:

   G4VPVDivisionFactory();

  protected:

   static G4VPVDivisionFactory* fgInstance;
     
};

#endif
