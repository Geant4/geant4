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
// class G4tgbVParallelWorldCreator
//
// Class description:
//
// Abstract class to tell G4tgbParallelGeomMgr which G4VUserParallelWorld's are instantiated
// History:
// - Created.                                 P.Arce, CIEMAT (March 2009)
// -------------------------------------------------------------------------

#ifndef G4tgbVParallelWorldCreator_h
#define G4tgbVParallelWorldCreator_h

#include "globals.hh"
class G4VUserParallelWorld;
#include <set>
#include <vector>

//----------------------------------------------------------------------------  
class G4tgbVParallelWorldCreator 
{
  public: 
    G4tgbVParallelWorldCreator(){};

    ~G4tgbVParallelWorldCreator(){};

    virtual G4VUserParallelWorld* CreateParallelWorld( const G4String& parallelWorldName, const G4int index ) = 0;

  private:
};

#endif
