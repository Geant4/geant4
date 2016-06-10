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
// $Id: G4tgrRotationMatrixFactory.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrRotationMatrixFactory
//
// Class description:
//
// Singleton class to manage the building of transient rotation matrix.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrRotationMatrixFactory_h
#define G4tgrRotationMatrixFactory_h

#include "globals.hh"
#include "G4tgrRotationMatrix.hh"

#include <vector>
#include <map>

typedef std::map< G4String, G4tgrRotationMatrix* > G4mstgrrotm;

class G4tgrRotationMatrixFactory 
{
  public:  // with desctiption

    static G4tgrRotationMatrixFactory* GetInstance();
      // Get the only instance (it it does not exists, create it)
  
    G4tgrRotationMatrix* AddRotMatrix( const std::vector<G4String>& wl );
      // Build a G4tgrRotationMatrix and add it to theTgrRotMats

    G4tgrRotationMatrix* FindRotMatrix(const G4String& rotm);
     // Look for an G4tgrRotationMatrix and if not found return 0

    const G4mstgrrotm& GetRotMatMap() const
      { return theTgrRotMats; }
    std::vector<G4tgrRotationMatrix*> GetRotMatList() const
      { return theTgrRotMatList; }

    void DumpRotmList();
      // Dump list of rotation matrices

  private:

    G4tgrRotationMatrixFactory();
   ~G4tgrRotationMatrixFactory();

  private:

    static G4ThreadLocal G4tgrRotationMatrixFactory* theInstance;
   
    std::vector<G4tgrRotationMatrix*> theTgrRotMatList;
    G4mstgrrotm theTgrRotMats;
};

#endif
