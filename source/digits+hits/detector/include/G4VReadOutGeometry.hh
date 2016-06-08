// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VReadOutGeometry.hh,v 1.1 1999/01/07 16:06:25 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
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
      G4VReadOutGeometry(const G4VReadOutGeometry &right);
      G4VReadOutGeometry(G4String);
      virtual ~G4VReadOutGeometry();

      const G4VReadOutGeometry & operator=(const G4VReadOutGeometry &right);

      G4int operator==(const G4VReadOutGeometry &right) const;
      G4int operator!=(const G4VReadOutGeometry &right) const;

      // buildROGeomety must be invoked to Build (ie Build() method)
      // the ROGeometry. It sets up in addition the needed
      // the G4Navigator used to navigate inside this ROGeometry.
      void BuildROGeometry();
      G4bool CheckROVolume(G4Step*,G4TouchableHistory*&);

  protected:
      G4bool FindROTouchable(G4Step*);

      G4VPhysicalVolume* ROworld;
      G4SensitiveVolumeList* fincludeList;
      G4SensitiveVolumeList* fexcludeList;
      G4String name;

  private:
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

