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
// $Id: G4VReadOutGeometry.hh 81087 2014-05-20 15:44:27Z gcosmo $
//
// ------------------------------------------------------------

#ifndef G4VReadOutGeometry_h
#define G4VReadOutGeometry_h

#include "G4SensitiveVolumeList.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"

class G4Navigator;

class G4VReadOutGeometry 
{
  protected:
      virtual G4VPhysicalVolume* Build() = 0;  // must return the world of the ROGeometry;

  public:
      G4VReadOutGeometry();
      G4VReadOutGeometry(G4String);
      virtual ~G4VReadOutGeometry();

      G4int operator==(const G4VReadOutGeometry &right) const;
      G4int operator!=(const G4VReadOutGeometry &right) const;

      // buildROGeomety must be invoked to Build (ie Build() method)
      // the ROGeometry. It sets up in addition the needed
      // the G4Navigator used to navigate inside this ROGeometry.
      void BuildROGeometry();
      virtual G4bool CheckROVolume(G4Step*,G4TouchableHistory*&);

  protected:
      G4VReadOutGeometry(const G4VReadOutGeometry &right);
      G4VReadOutGeometry & operator=(const G4VReadOutGeometry &right);

      virtual G4bool FindROTouchable(G4Step*);

  protected:
      G4VPhysicalVolume* ROworld;
      G4SensitiveVolumeList* fincludeList;
      G4SensitiveVolumeList* fexcludeList;
      G4String name;

      G4Navigator*        ROnavigator;
      G4TouchableHistory* touchableHistory;

  public:
      inline const G4SensitiveVolumeList* GetIncludeList() const
      { return fincludeList; }
      inline void SetIncludeList(G4SensitiveVolumeList* value)
      { fincludeList = value; }
      inline const G4SensitiveVolumeList* GetExcludeList() const
      { return fexcludeList; }
      inline void SetExcludeList(G4SensitiveVolumeList* value)
      { fexcludeList = value; }
      inline G4String GetName() const
      { return name; }
      inline void SetName(G4String value)
      { name = value; }
      // ADDED:
      inline G4VPhysicalVolume* GetROWorld() const
      { return ROworld;}
};


#endif

