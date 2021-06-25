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
// G4PhysicsListOrderingParameter
// 
// Class description:
//
// This class defins a parameter ordering used by G4PhysicsListHelper.

// Author: H.Kurashige, 29 April 2011
// --------------------------------------------------------------------
#ifndef G4PhysicsListOrderingParameter_hh
#define G4PhysicsListOrderingParameter_hh 1

#include "G4ios.hh"
#include "globals.hh"

class G4PhysicsListHelper;
class G4PhysicsListOrderingParameter
{
  friend class G4PhysicsListHelper;

  public:

    G4PhysicsListOrderingParameter();
    virtual ~G4PhysicsListOrderingParameter();

    inline const G4String& GetTypeName() const { return processTypeName; }
    inline G4int GetType() const { return processType; }
    inline G4int GetSubType() const { return processSubType; }
    inline G4int GetOrdering(G4int idx) const
    {
      return ((idx < -1) || (idx > 2)) ? -1 : ordering[idx];
    }
    inline G4bool GetDuplicable() const { return isDuplicable; }

  private:

    G4String processTypeName = "NONE";
    G4int processType = -1;
    G4int processSubType = -1;
    G4int ordering[3];
    G4bool isDuplicable = false;
};

#endif
