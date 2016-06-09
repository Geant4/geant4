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
// $Id: G4PVDivisionFactory.cc,v 1.1 2004/05/13 14:56:31 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// class G4PVDivisionFactory Implementation file
//
// Author: Ivana Hrivnacova, 4.5.2004  (Ivana.Hrivnacova@cern.ch)
// --------------------------------------------------------------------

#include "G4PVDivisionFactory.hh"
#include "G4PVDivision.hh"
#include "G4VDivisionParameterisation.hh"

//_____________________________________________________________________________

G4PVDivisionFactory::G4PVDivisionFactory()
  : G4VPVDivisionFactory()
{
  // Protected singleton constructor.
  // ---
}

//_____________________________________________________________________________

G4PVDivisionFactory::~G4PVDivisionFactory()
{
}

//_____________________________________________________________________________

G4PVDivisionFactory* G4PVDivisionFactory::GetInstance()
{
  static G4PVDivisionFactory theFactory;
  if (!fgInstance)
  {
    fgInstance = &theFactory;
  }
  return &theFactory;
} 

//_____________________________________________________________________________

G4VPhysicalVolume* 
G4PVDivisionFactory::CreatePVDivision(const G4String& pName,
                             G4LogicalVolume* pLogical,
                             G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double width,
                             const G4double offset )
{     
  // Create division - with number of divisions and width
  // ---

  return new G4PVDivision(pName, pLogical, pMotherLogical, 
                          pAxis, nReplicas, width, offset);
}    

//_____________________________________________________________________________

G4VPhysicalVolume* 
G4PVDivisionFactory::CreatePVDivision(const G4String& pName,
                             G4LogicalVolume* pLogical,
                             G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double offset )
{     
  // Create division - with number of divisions 
  // ---

  return new G4PVDivision(pName, pLogical, pMotherLogical, 
                          pAxis, nReplicas, offset);
}    

//_____________________________________________________________________________

G4VPhysicalVolume* 
G4PVDivisionFactory::CreatePVDivision(const G4String& pName,
                             G4LogicalVolume* pLogical,
                             G4LogicalVolume* pMotherLogical,
                             const EAxis pAxis,
                             const G4double width,
                             const G4double offset )
{     
  // Create division - with width
  // ---

  return new G4PVDivision(pName, pLogical, pMotherLogical, 
                          pAxis, width, offset);
}    

//_____________________________________________________________________________

G4VPhysicalVolume* 
G4PVDivisionFactory::CreatePVDivision(const G4String& pName,
                             G4LogicalVolume* pLogical,
                             G4LogicalVolume* pMotherLogical,
                             const G4VPVParameterisation* param)
{     
  // Create division - with parameterisation
  // ---

  // Get parameterisation data
  //
  const G4VDivisionParameterisation* divParam
    = dynamic_cast<const G4VDivisionParameterisation*>(param);

  if (!divParam)
  {
    G4Exception("G4PVDivisionFactory::CreatePVDivision()",
                "WrongType", FatalException,
                "Unexpected parameterisation type !");
  }

  EAxis axis = divParam->GetAxis();
  G4int nofDivisions = divParam->GetNoDiv();
  G4double width = divParam->GetWidth();
  G4double offset = divParam->GetOffset();

  return new G4PVDivision(pName, pLogical, pMotherLogical, 
                          axis, nofDivisions, width, offset);
}    

//_____________________________________________________________________________

G4bool G4PVDivisionFactory::IsPVDivision(const G4VPhysicalVolume* pv) const
{ 
  // Returns true if pv is division
  // ---

  if (dynamic_cast<const G4PVDivision*>(pv))
    return true;
  else
    return false;  
}

