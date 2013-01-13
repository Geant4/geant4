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
// class G4TouchableHistory Implementation
//
// ----------------------------------------------------------------------

#include "G4TouchableHistory.hh"

__thread G4Allocator<G4TouchableHistory> *aTouchableHistoryAllocator_G4MT_TLS_ = 0;

G4TouchableHistory::G4TouchableHistory()
  : frot(G4RotationMatrix()),
    ftlate(G4ThreeVector(0.,0.,0.)),
    fhistory()
{  ;;;   if (!aTouchableHistoryAllocator_G4MT_TLS_) aTouchableHistoryAllocator_G4MT_TLS_ = new G4Allocator<G4TouchableHistory>  ; G4Allocator<G4TouchableHistory> &aTouchableHistoryAllocator = *aTouchableHistoryAllocator_G4MT_TLS_;  ;;;  
   G4VPhysicalVolume* pPhysVol=0;
   fhistory.SetFirstEntry(pPhysVol);
}

G4TouchableHistory::G4TouchableHistory( const G4NavigationHistory &history )
  : fhistory(history)
{  ;;;   if (!aTouchableHistoryAllocator_G4MT_TLS_) aTouchableHistoryAllocator_G4MT_TLS_ = new G4Allocator<G4TouchableHistory>  ; G4Allocator<G4TouchableHistory> &aTouchableHistoryAllocator = *aTouchableHistoryAllocator_G4MT_TLS_;  ;;;  
  G4AffineTransform tf(fhistory.GetTopTransform().Inverse());
  ftlate = tf.NetTranslation();
  frot = tf.NetRotation();
}

G4TouchableHistory::~G4TouchableHistory()
{  ;;;   if (!aTouchableHistoryAllocator_G4MT_TLS_) aTouchableHistoryAllocator_G4MT_TLS_ = new G4Allocator<G4TouchableHistory>  ; G4Allocator<G4TouchableHistory> &aTouchableHistoryAllocator = *aTouchableHistoryAllocator_G4MT_TLS_;  ;;;  
}

const G4ThreeVector&
G4TouchableHistory::GetTranslation(G4int depth) const
{  ;;;   if (!aTouchableHistoryAllocator_G4MT_TLS_) aTouchableHistoryAllocator_G4MT_TLS_ = new G4Allocator<G4TouchableHistory>  ; G4Allocator<G4TouchableHistory> &aTouchableHistoryAllocator = *aTouchableHistoryAllocator_G4MT_TLS_;  ;;;  
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  static __thread G4ThreeVector *currTranslation_G4MT_TLS_ = 0 ; if (!currTranslation_G4MT_TLS_) currTranslation_G4MT_TLS_ = new  G4ThreeVector  ;  G4ThreeVector &currTranslation = *currTranslation_G4MT_TLS_;
  if(depth==0.0)
  {
    return ftlate;
  }
  else
  {
    currTranslation =
      fhistory.GetTransform(CalculateHistoryIndex(depth)).NetTranslation();
    return currTranslation;
  }
}

const G4RotationMatrix*
G4TouchableHistory::GetRotation(G4int depth) const
{  ;;;   if (!aTouchableHistoryAllocator_G4MT_TLS_) aTouchableHistoryAllocator_G4MT_TLS_ = new G4Allocator<G4TouchableHistory>  ; G4Allocator<G4TouchableHistory> &aTouchableHistoryAllocator = *aTouchableHistoryAllocator_G4MT_TLS_;  ;;;  
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  static __thread G4RotationMatrix *rotM_G4MT_TLS_ = 0 ; if (!rotM_G4MT_TLS_) rotM_G4MT_TLS_ = new  G4RotationMatrix  ;  G4RotationMatrix &rotM = *rotM_G4MT_TLS_;

  if(depth==0.0)
  {
    return &frot;
  }
  else
  {
    rotM = fhistory.GetTransform(CalculateHistoryIndex(depth)).NetRotation();
    return &rotM;
  }
}
