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
#ifndef PAR04HIT_HH
#define PAR04HIT_HH

#include <stddef.h>                      // for size_t
#include <G4Types.hh>                    // for G4int, G4double
#include <map>                           // for map
#include <tls.hh>                        // for G4ThreadLocal
#include <vector>                        // for vector
#include "G4Allocator.hh"                // for G4Allocator
#include "G4RotationMatrix.hh"           // for G4RotationMatrix
#include "G4THitsCollection.hh"          // for G4THitsCollection
#include "G4ThreeVector.hh"              // for G4ThreeVector
#include "G4VHit.hh"                     // for G4VHit
class G4AttDef;
class G4AttValue;
class G4LogicalVolume;
class G4String;

/**
 * @brief Hit class to store energy deposited in the sensitive detector.
 *
 * Hit class registers position and energy deposited within the sensitive
 * detector. Cell ID is stored using identifiers of readout segmentation (z,
 * phi, rho). Additionally, pointer to cell logical volume, its position and
 * rotation are saved for visualisation purposes. Time allows to filter hits in
 * visualisation. Type of hit allows to distinguish between hits originating
 * from full simulation (type 0) and fast simulation (type 1).
 *
 */

class Par04Hit : public G4VHit
{
 public:
  Par04Hit();
  Par04Hit(const Par04Hit& aRight);
  virtual ~Par04Hit();

  const Par04Hit& operator=(const Par04Hit& aRight);
  int operator==(const Par04Hit& aRight) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);
  /// Visualise hits. If pointer to the logical volume was set, cell shape is
  /// drawn taking into account proper radial position (taken from fRhoId)
  virtual void Draw() final;
  /// Retrieve atributes' names in order to allow filtering
  virtual const std::map<G4String, G4AttDef>* GetAttDefs() const final;
  /// Create attributes for the visualisation.
  virtual std::vector<G4AttValue>* CreateAttValues() const final;
  /// Print hit properties.
  virtual void Print() final;
  /// Set position
  inline void SetPos(G4ThreeVector aXYZ) { fPos = aXYZ; }
  /// Get position
  inline G4ThreeVector GetPos() const { return fPos; }
  /// Set rotation
  inline void SetRot(G4RotationMatrix aXYZ) { fRot = aXYZ; }
  /// Get rotation
  inline G4RotationMatrix GetRot() const { return fRot; }
  /// Set energy
  inline void SetEdep(G4double aEdep) { fEdep = aEdep; }
  /// Add energy to previous value
  inline void AddEdep(G4double aEdep) { fEdep += aEdep; }
  /// Get energy
  inline G4double GetEdep() const { return fEdep; }
  /// Set Z id of the cell in the readout segmentation
  inline void SetZid(G4int aZ) { fZId = aZ; }
  /// Get Z id of the cell in the readout segmentation
  inline G4int GetZid() const { return fZId; }
  /// Set Rho id of the cell in the readout segmentation
  inline void SetRhoId(G4int aRho) { fRhoId = aRho; }
  /// Get rho id of the cell in the readout segmentation
  inline G4int GetRhoId() const { return fRhoId; }
  /// Set phi id of the cell in the readout segmentation
  inline void SetPhiId(G4int aPhi) { fPhiId = aPhi; }
  /// Get phi id of the cell in the readout segmentation
  inline G4int GetPhiId() const { return fPhiId; }
  /// Set time
  inline void SetTime(G4double aTime) { fTime = aTime; }
  /// Get time
  inline G4double GetTime() const { return fTime; }
  /// Set type (0 = full sim, 1 = fast sim)
  inline void SetType(G4int aType) { fType = aType; }
  /// Get type (0 = full sim, 1 = fast sim)
  inline G4int GetType() const { return fType; }
  // Set pointer to cell logical volume
  inline void SetLogV(G4LogicalVolume* aLogVol) { fLogVol = aLogVol; }
  // Get pointer to cell logical volume
  inline const G4LogicalVolume* GetLogVol() { return fLogVol; }

 public:
  /// Energy deposit
  G4double fEdep = 0;
  /// Z ID of readout cell
  G4int fZId = -1;
  /// Rho ID of readout cell
  G4int fRhoId = -1;
  /// Phi ID of readout cell
  G4int fPhiId = -1;
  /// Position
  G4ThreeVector fPos = { -1, -1, -1 };
  /// Rotation
  G4RotationMatrix fRot;
  /// Time
  G4double fTime = -1;
  /// Type: 0 = full sim, 1 = fast sim
  G4int fType = -1;
  /// Pointer to logical volume for visualisation
  G4LogicalVolume* fLogVol = nullptr;
};

typedef G4THitsCollection<Par04Hit> Par04HitsCollection;

extern G4ThreadLocal G4Allocator<Par04Hit>* Par04HitAllocator;

inline void* Par04Hit::operator new(size_t)
{
  if(!Par04HitAllocator)
    Par04HitAllocator = new G4Allocator<Par04Hit>;
  return (void*) Par04HitAllocator->MallocSingle();
}

inline void Par04Hit::operator delete(void* aHit)
{
  Par04HitAllocator->FreeSingle((Par04Hit*) aHit);
}

#endif /* PAR04HIT_HH */
