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
// class G4tgbParallelGeomMgr
//
// Class description:
//
// Class to manage the building of parallel worlds. It is a singleton.

// History:
// - Created.                                 P.Arce, CIEMAT (February 2009)
// -------------------------------------------------------------------------

#ifndef G4tgbParallelGeomMgr_h
#define G4tgbParallelGeomMgr_h

#include "globals.hh"
class G4VUserParallelWorld;
#include <set>
#include <vector>
class G4tgbVParallelWorldCreator;

//----------------------------------------------------------------------------  
class G4tgbParallelGeomMgr 
{
  protected:
    G4tgbParallelGeomMgr();

  public:  // with description  
    ~G4tgbParallelGeomMgr();

    static G4tgbParallelGeomMgr* GetInstance();  
      // Get the only instance 

    void AddParallelWorldIndex( G4int index );

    std::vector<G4VUserParallelWorld*> CreateParalleWorlds();

    void BuildPhysicsProcess( const G4String& volName, const G4int index );

    std::vector<G4VUserParallelWorld*> GetParallelWorlds() const { 
      return theParallelWorlds; 
    }

  G4tgbVParallelWorldCreator* GetWorldCreator() const {
    return theWorldCreator; 
  }

  void SetWorldCreator(G4tgbVParallelWorldCreator* wc ){
    theWorldCreator = wc;
  }

  private:

    static G4tgbParallelGeomMgr* theInstance;
    std::set<G4int> theIndices;
    std::vector<G4VUserParallelWorld*> theParallelWorlds;

    G4tgbVParallelWorldCreator* theWorldCreator;
};

#endif
