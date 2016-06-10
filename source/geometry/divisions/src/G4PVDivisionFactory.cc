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
// $Id: G4PVDivisionFactory.cc 68718 2013-04-05 09:16:24Z gcosmo $
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
  if (!fgInstance)
  {
    fgInstance =  new G4PVDivisionFactory;
  }
  return dynamic_cast<G4PVDivisionFactory*>(fgInstance);
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
                "GeomDiv0001", FatalException,
                "Unexpected parameterisation type!");
    return 0;
  }
  else
  {
    EAxis axis = divParam->GetAxis();
    G4int nofDivisions = divParam->GetNoDiv();
    G4double width = divParam->GetWidth();
    G4double offset = divParam->GetOffset();

    return new G4PVDivision(pName, pLogical, pMotherLogical, 
                            axis, nofDivisions, width, offset);
  }
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

