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
// $Id: G4PVDivisionFactory.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
