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
// $Id$
//
// 
// class G4VPhysicalVolume Implementation
//
// --------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//This static member is thread local. For each thread, it points to the
//array of PhysicalVolumePrivateSubclass instances.
template <class PhysicalVolumePrivateSubclass> __thread PhysicalVolumePrivateSubclass* G4MTPrivateSubInstanceManager<PhysicalVolumePrivateSubclass>::offset = 0;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//This new field helps to use the class G4VPhysicalVolumeSubInstanceManager
//introduced in the "G4VPhysicalVolume.hh" file.
G4VPhysicalVolumeSubInstanceManager G4VPhysicalVolume::g4vphysicalVolumeSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//This method is similar to the constructor. It is used by each worker
//thread to achieve the same effect as that of the master thread exept
//to register the new created instance. This method is invoked explicitly.
//It does not create a new G4VPhysicalVolume instance. It only assign the value
//for the fields encapsulated by the class PhysicalVolumePrivateSubclass.
void G4VPhysicalVolume::SlaveG4VPhysicalVolume( G4VPhysicalVolume* /*pMasterObject*/, G4RotationMatrix *pRot,
				      const G4ThreeVector &tlate)
{
  g4vphysicalVolumeSubInstanceManager.SlaveCopySubInstanceArray();

  frotG4MTThreadPrivate = pRot;
  ftransG4MTThreadPrivate = tlate;
  //  G4PhysicalVolumeStore::Register(this);
}

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//This method is similar to the destructor. It is used by each worker
//thread to achieve the partial effect as that of the master thread.
//For G4VPhysicalVolume instances, nothing more to do here.
void G4VPhysicalVolume::DestroySlaveG4VPhysicalVolume( G4VPhysicalVolume* /*pMasterObject*/)
{
}

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
  g4vphysicalVolumeInstanceID = g4vphysicalVolumeSubInstanceManager.CreateSubInstance();
  frotG4MTThreadPrivate = pRot;
  ftransG4MTThreadPrivate = tlate; 
  G4PhysicalVolumeStore::Register(this);
}



// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VPhysicalVolume::G4VPhysicalVolume( __void__& )
  : flogical(0), fname(""), flmother(0)
{
  // Register to store
  //
  g4vphysicalVolumeInstanceID = g4vphysicalVolumeSubInstanceManager.CreateSubInstance();

  frotG4MTThreadPrivate = 0; 
  G4PhysicalVolumeStore::Register(this);
}

// Destructor -  remove from Store
//
G4VPhysicalVolume::~G4VPhysicalVolume() 
{
  G4PhysicalVolumeStore::DeRegister(this);
}

G4int G4VPhysicalVolume::GetMultiplicity() const
{
  return 1;
}

G4RotationMatrix* G4VPhysicalVolume::GetObjectRotation() const
{
  static __thread G4RotationMatrix  *aRotM_G4MT_TLS_ = 0 ; if (!aRotM_G4MT_TLS_) aRotM_G4MT_TLS_ = new  G4RotationMatrix   ;  G4RotationMatrix  &aRotM = *aRotM_G4MT_TLS_; 
  static __thread G4RotationMatrix  *IdentityRM_G4MT_TLS_ = 0 ; if (!IdentityRM_G4MT_TLS_) IdentityRM_G4MT_TLS_ = new  G4RotationMatrix   ;  G4RotationMatrix  &IdentityRM = *IdentityRM_G4MT_TLS_;  // Never changed (from "1")
  G4RotationMatrix* retval; 

  // Insure against frot being a null pointer
  if(frotG4MTThreadPrivate)
  {
    aRotM= frotG4MTThreadPrivate->inverse();
    retval= &aRotM;
  }
  else
  {
    retval= &IdentityRM;
  }
  return retval;
}

// Only implemented for placed and parameterised volumes.
// Not required for replicas.
//
G4bool G4VPhysicalVolume::CheckOverlaps(G4int, G4double, G4bool)
{
  return false;
}
