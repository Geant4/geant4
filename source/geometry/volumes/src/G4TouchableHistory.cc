// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TouchableHistory.cc,v 1.4 2000-11-01 16:51:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory Implementation

#include "G4TouchableHistory.hh"

G4TouchableHistory::~G4TouchableHistory()
{
}

const G4ThreeVector&
G4TouchableHistory::GetTranslation(G4int depth) const
{
   // The value returned will change at the next call
   //  Copy it if you want to use it!
   //   
   static G4ThreeVector currTranslation;
   if(depth==0.0) {
      return ftlate;
   }else{
      currTranslation =
         fhistory.GetTransform(CalculateHistoryIndex(depth)).NetTranslation();
      return currTranslation;
   }
}

const G4RotationMatrix*
G4TouchableHistory::GetRotation(G4int depth) const
{
   // The value returned will change at the next call
   //  Copy it if you want to use it!
   //
   static G4RotationMatrix rotM;

   if(depth==0.0) {
      return &frot;
   }else{
      rotM = fhistory.GetTransform(CalculateHistoryIndex(depth)).NetRotation();
      return &rotM;
   }
}
