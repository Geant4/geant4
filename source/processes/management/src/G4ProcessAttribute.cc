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
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1997, H.Kurashige
//   ----------------  G4ProcessAttribute -----------------
// History:
//   adds copy constructor            27 June 1998 H.Kurashige
// ------------------------------------------------------------

#include "G4ProcessAttribute.hh"

// default constructor //////////////////////////
G4ProcessAttribute::G4ProcessAttribute():
         pProcess(nullptr),
	 isActive(true),
	 idxProcessList(-1)
{
  // clear 
  for (G4int idx=0; idx<G4ProcessManager::SizeOfProcVectorArray; idx++){
    idxProcVector[idx] = -1;
    ordProcVector[idx] = -1;
  }
}

// copy constructor //////////////////////////
G4ProcessAttribute::G4ProcessAttribute(const G4ProcessAttribute &right):
         pProcess(right.pProcess),
	 isActive(right.isActive),
	 idxProcessList(right.idxProcessList)
{
  // copy all contents in idxProcVector[] and ordProcVector[]
  //   deep copy 
  for (G4int idx=0; idx<G4ProcessManager::SizeOfProcVectorArray; idx++){
    idxProcVector[idx] = right.idxProcVector[idx];
    ordProcVector[idx] = right.ordProcVector[idx];
  }
}

// destructor   ///////////////////////
G4ProcessAttribute::~G4ProcessAttribute()
{
   // do nothing
}

// assignment operator //////////
G4ProcessAttribute & G4ProcessAttribute::operator=(const G4ProcessAttribute &right)
{
  if (this != &right) {
    pProcess       = right.pProcess;
    idxProcessList = right.idxProcessList;
    isActive = right.isActive;
    // copy all contents in idxProcVector[] and ordProcVector[]
    //   deep copy 
    for (G4int idx=0; idx<G4ProcessManager::SizeOfProcVectorArray; idx++){
      idxProcVector[idx] = right.idxProcVector[idx];
      ordProcVector[idx] = right.ordProcVector[idx];
    }
  }
  return *this;
}



