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
// class G4TouchableHistory Implementation
//
// Created: Paul Kent, August 1996
// ----------------------------------------------------------------------

#include "G4TouchableHistory.hh"

G4Allocator<G4TouchableHistory>*& aTouchableHistoryAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4TouchableHistory>* _instance = nullptr;
  return _instance;
}

G4TouchableHistory::G4TouchableHistory()
  : frot(G4RotationMatrix()),
    ftlate(G4ThreeVector(0.,0.,0.)),
    fhistory()
{ 
  fhistory.SetFirstEntry(nullptr);
}

G4TouchableHistory::G4TouchableHistory( const G4NavigationHistory& history )
  : fhistory(history)
{
  const G4AffineTransform& tf = fhistory.GetTopTransform();
  ftlate = tf.InverseNetTranslation();
  frot = tf.InverseNetRotation();
}

G4TouchableHistory::~G4TouchableHistory()
{ 
}

const G4ThreeVector&
G4TouchableHistory::GetTranslation(G4int depth) const
{ 
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  static G4ThreadLocal G4ThreeVector* ctrans = nullptr;
  if ( ctrans == nullptr )  { ctrans = new G4ThreeVector; }
  if(depth==0.0)
  {
    return ftlate;
  }
  else
  {
    *ctrans =
      fhistory.GetTransform(CalculateHistoryIndex(depth)).NetTranslation();
    return *ctrans;
  }
}

const G4RotationMatrix*
G4TouchableHistory::GetRotation(G4int depth) const
{ 
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  static G4ThreadLocal G4RotationMatrix* rotM = nullptr;
  if ( rotM == nullptr )  { rotM = new G4RotationMatrix; }

  if(depth==0.0)
  {
    return &frot;
  }
  else
  {
    *rotM = fhistory.GetTransform(CalculateHistoryIndex(depth)).NetRotation();
    return rotM;
  }
}
