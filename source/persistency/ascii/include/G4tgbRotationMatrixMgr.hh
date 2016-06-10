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
// $Id: G4tgbRotationMatrixMgr.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgbRotationMatrixMgr
//
// Class description:
//
// Singleton class to manage the building of transient rotation matrix,
// as well as the construction of the corresponding G4RotationMatrix's.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbRotationMatrixMgr_h
#define G4tgbRotationMatrixMgr_h

#include "globals.hh"

#include <iostream>
#include <map>

#include "G4tgbRotationMatrix.hh"

typedef std::map< G4String, G4tgbRotationMatrix*,
                  std::less<G4String> > G4mstgbrotm;
typedef std::map< G4String, G4RotationMatrix*,
                  std::less<G4String> > G4msg4rotm;

class G4tgbRotationMatrixMgr 
{
  public:  // with description

    ~G4tgbRotationMatrixMgr();
  
    static G4tgbRotationMatrixMgr* GetInstance();
      // Get only instance (if it does not exists, create it)

    G4RotationMatrix* FindOrBuildG4RotMatrix(const G4String& name);
      // Look for a G4RotationMatrix and if not found create it
       // from the corresponding G4tgbRotationMatrix
    G4RotationMatrix* FindG4RotMatrix(const G4String& name);
      // Look for a G4RotationMatrix and if not found return 0

    G4tgbRotationMatrix* FindOrBuildTgbRotMatrix(const G4String& name);
      // Look for an G4tgbRotationMatrix and if not found exit
    G4tgbRotationMatrix* FindTgbRotMatrix(const G4String& name);
      // Look for an G4tgbRotationMatrix and if not found return 0

  public:  // without description

    const G4mstgbrotm GetTgbRotMatList() const { return theTgbRotMats; }
    const G4msg4rotm& GetG4RotMatList() const { return theG4RotMats; }

  private:

    G4tgbRotationMatrixMgr();
    void CopyRotMats();

    static G4ThreadLocal G4tgbRotationMatrixMgr* theInstance;
   
  private:

    G4mstgbrotm theTgbRotMats;
    G4msg4rotm theG4RotMats; 
};

std::ostream& operator<<(std::ostream&, const G4RotationMatrix &);

#endif
