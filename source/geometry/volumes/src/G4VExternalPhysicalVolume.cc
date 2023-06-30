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
// G4VExternalPhysicalVolume Implementation
//
// Author: J.Apostolakis, October 2019
// ----------------------------------------------------------------------

#include "G4VExternalPhysicalVolume.hh"
#include "G4AffineTransform.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

// ----------------------------------------------------------------------
// Constructor
//
G4VExternalPhysicalVolume::
G4VExternalPhysicalVolume( G4RotationMatrix* pRot,
                     const G4ThreeVector& tlate,
                           G4LogicalVolume* pLogical,
                     const G4String& pName,
                           G4VPhysicalVolume* pMother )
   : G4VPhysicalVolume(pRot, tlate, pName, pLogical, pMother)
{
  if (pMother != nullptr)
  {
    G4LogicalVolume* motherLogical = pMother->GetLogicalVolume();
    if (pLogical == motherLogical)
    {
      G4Exception("G4VExternalPhysicalVolume::G4VExternalPhysicalVolume()",
                  "GeomVol0002", FatalException,
                  "Cannot place a volume inside itself!");
    }
    SetMotherLogical(motherLogical);
    motherLogical->AddDaughter(this);
  }
}

// ----------------------------------------------------------------------
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VExternalPhysicalVolume::G4VExternalPhysicalVolume( __void__& a )
  : G4VPhysicalVolume(a)
{
}

// ----------------------------------------------------------------------
// Destructor
//
G4VExternalPhysicalVolume::~G4VExternalPhysicalVolume() = default;

// ----------------------------------------------------------------------
// IsMany
//
G4bool G4VExternalPhysicalVolume::IsMany() const
{
  return fMany;
}

// ----------------------------------------------------------------------
// SetMany
//
void G4VExternalPhysicalVolume::SetMany(G4bool overlap)
{
  fMany = overlap;
}

// ----------------------------------------------------------------------
// IsReplicated
//
G4bool G4VExternalPhysicalVolume::IsReplicated() const
{
  return false;
}

// ----------------------------------------------------------------------
// IsParameterised
//
G4bool G4VExternalPhysicalVolume::IsParameterised() const
{
  return false;
}

// ----------------------------------------------------------------------
// GetParameterisation
//
G4VPVParameterisation* G4VExternalPhysicalVolume::GetParameterisation() const
{
  return nullptr;
}

// ----------------------------------------------------------------------
// GetReplicationData
//
void G4VExternalPhysicalVolume::
GetReplicationData( EAxis&, G4int&, G4double&, G4double&, G4bool& ) const
{
  // No-operations
}

// ----------------------------------------------------------------------
// IsRegularRepeatedStructure
//
// This is for specialised repeated volumes (replicas, parameterised vol.)
//
G4bool G4VExternalPhysicalVolume::IsRegularStructure() const
{
  return false;
}           

// ----------------------------------------------------------------------
// IsRegularRepeatedStructure
//
// This is for specialised repeated volumes (replicas, parameterised vol.)
//
G4int G4VExternalPhysicalVolume::GetRegularStructureId() const
{
  return 0;  
}           

// ----------------------------------------------------------------------
// VolumeType
//
// Information to help identify sub-navigator which will be used.
//
EVolume G4VExternalPhysicalVolume::VolumeType() const
{
  return kExternal;
}
