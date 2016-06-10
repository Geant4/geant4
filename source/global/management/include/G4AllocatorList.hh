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
// $Id: G4AllocatorList.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 
// ------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// A class to store all G4Allocator objects in a thread for the sake
// of cleanly deleting them.
//
// ------------------------------------------------------------

#ifndef G4AllocatorList_h
#define G4AllocatorList_h 1

#include <vector>
#include "globals.hh"

class G4AllocatorBase;

class G4AllocatorList
{
  public:  // with description

    static G4AllocatorList* GetAllocatorList();
    static G4AllocatorList* GetAllocatorListIfExist();

  public:

    ~G4AllocatorList();
    void Register(G4AllocatorBase*);
    void Destroy(G4int nStat=0, G4int verboseLevel=0);
    G4int Size() const;

  private:

    G4AllocatorList();

  private:

    static G4ThreadLocal G4AllocatorList* fAllocatorList;
    std::vector<G4AllocatorBase*> fList;
};

#endif
