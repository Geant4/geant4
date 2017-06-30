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
// $Id: G4VPhysicalVolume.cc 103096 2017-03-15 15:21:33Z gcosmo $
//
// 
// class G4VPhysicalVolume Implementation
//
// --------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

// This new field helps to use the class G4PVManager
//
G4PVManager G4VPhysicalVolume::subInstanceManager;

// These macros change the references to fields that are now encapsulated
// in the class G4PVData.
//
#define G4MT_rot ((subInstanceManager.offset[instanceID]).frot)
#define G4MT_trans ((subInstanceManager.offset[instanceID]).ftrans)
#define G4MT_pvdata (subInstanceManager.offset[instanceID])

// Constructor: init parameters and register in Store
//
G4VPhysicalVolume::G4VPhysicalVolume( G4RotationMatrix *pRot,
                                const G4ThreeVector &tlate,
                                const G4String& pName,
                                      G4LogicalVolume* pLogical,
                                      G4VPhysicalVolume* )
  : flogical(pLogical),
    fname(pName), flmother(0)
{
  instanceID = subInstanceManager.CreateSubInstance();

  this->SetRotation( pRot );       // G4MT_rot = pRot;
  this->SetTranslation( tlate );   // G4MT_trans = tlate;

  // Initialize 'Shadow' data structure - for use by object persistency
  pvdata = new G4PVData();
  pvdata->frot = pRot;
  pvdata->ftrans = G4ThreeVector(tlate);

  G4PhysicalVolumeStore::Register(this);
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VPhysicalVolume::G4VPhysicalVolume( __void__& )
  : flogical(0), fname(""), flmother(0), pvdata(0)
{
  // Register to store
  //
  instanceID = subInstanceManager.CreateSubInstance();

  G4PhysicalVolumeStore::Register(this);
}

// Destructor -  remove from Store
//
G4VPhysicalVolume::~G4VPhysicalVolume() 
{
  delete pvdata;
  G4PhysicalVolumeStore::DeRegister(this);
}

// This method is similar to the constructor. It is used by each worker
// thread to achieve the same effect as that of the master thread exept
// to register the new created instance. This method is invoked explicitly.
// It does not create a new G4VPhysicalVolume instance.
// It only assign the value for the fields encapsulated by the class G4PVData.
//
void G4VPhysicalVolume::
InitialiseWorker( G4VPhysicalVolume* /*pMasterObject*/,
                  G4RotationMatrix *pRot,
                  const G4ThreeVector &tlate)
{
  subInstanceManager.SlaveCopySubInstanceArray();

  this->SetRotation( pRot );      // G4MT_rot   = pRot;
  this->SetTranslation( tlate );  // G4MT_trans = tlate;
  //  G4PhysicalVolumeStore::Register(this);
}

// Release memory allocated for offset
//
void G4VPhysicalVolume::Clean()
{
  subInstanceManager.FreeSlave();
}

// This method is similar to the destructor. It is used by each worker
// thread to achieve the partial effect as that of the master thread.
// For G4VPhysicalVolume instances, nothing more to do here.
//
void G4VPhysicalVolume::TerminateWorker( G4VPhysicalVolume* /*pMasterObject*/)
{
}

// Returns the private data instance manager.
//
const G4PVManager& G4VPhysicalVolume::GetSubInstanceManager()
{
  return subInstanceManager;
}

G4int G4VPhysicalVolume::GetMultiplicity() const
{
  return 1;
}

const G4ThreeVector& G4VPhysicalVolume::GetTranslation() const
{
  return G4MT_trans;
}

void G4VPhysicalVolume::SetTranslation(const G4ThreeVector &vec)
{
  G4MT_trans=vec;
}

const G4RotationMatrix* G4VPhysicalVolume::GetRotation() const
{
  return G4MT_rot;
}

G4RotationMatrix* G4VPhysicalVolume::GetRotation()
{
  return G4MT_rot;
}

void G4VPhysicalVolume::SetRotation(G4RotationMatrix *pRot)
{
  G4MT_rot=pRot;
}

G4RotationMatrix* G4VPhysicalVolume::GetObjectRotation() const
{
  static G4RotationMatrix aRotM;
  static G4RotationMatrix IdentityRM;

  G4RotationMatrix* retval = &IdentityRM;

  // Insure against frot being a null pointer
  if(this->GetRotation())
  {
     aRotM = GetRotation()->inverse();
     retval= &aRotM;
  }
  return retval;
}

G4RotationMatrix G4VPhysicalVolume::GetObjectRotationValue() const
{
  G4RotationMatrix  aRotM;   // Initialised to identity

  // Insure against G4MT_rot being a null pointer
  if(G4MT_rot)
  {
     aRotM= G4MT_rot->inverse();
  }
  return aRotM;
}

G4ThreeVector  G4VPhysicalVolume::GetObjectTranslation() const
{
  return G4MT_trans;
}

const G4RotationMatrix* G4VPhysicalVolume::GetFrameRotation() const
{
  return G4MT_rot;
}

G4ThreeVector  G4VPhysicalVolume::GetFrameTranslation() const
{
  return -G4MT_trans;
}

// Only implemented for placed and parameterised volumes.
// Not required for replicas.
//
G4bool G4VPhysicalVolume::CheckOverlaps(G4int, G4double, G4bool, G4int)
{
  return false;
}
