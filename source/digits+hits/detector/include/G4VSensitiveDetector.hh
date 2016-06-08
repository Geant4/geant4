// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSensitiveDetector.hh,v 1.1 1999/01/07 16:06:25 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4VSensitiveDetector_h
#define G4VSensitiveDetector_h 1

#include "G4VHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include <rw/tvordvec.h>

class G4VSensitiveDetector 
{

  public:
      G4VSensitiveDetector(G4String name);
      G4VSensitiveDetector(const G4VSensitiveDetector &right);

      virtual ~G4VSensitiveDetector();

      const G4VSensitiveDetector & operator=(const G4VSensitiveDetector &right);

      int operator==(const G4VSensitiveDetector &right) const;
      int operator!=(const G4VSensitiveDetector &right) const;

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);
      virtual void clear();

      virtual void DrawAll();
      virtual void PrintAll();

  protected:
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) = 0;

      G4int GetCollectionID(G4int i);
      G4String SensitiveDetectorName; // detector name
      G4String thePathName;           // directory path
      G4String fullPathName;          // path + detector name
      G4int verboseLevel;
      G4bool active;
      RWTValOrderedVector<G4String> collectionName;

  private:
      G4VReadOutGeometry * ROgeometry;

  public:
      inline G4bool Hit(G4Step*aStep)
      {
        G4bool ack = true; 
        G4TouchableHistory* ROhis = NULL;
        if(!isActive()) 
        { ack = false; }
        else if(ROgeometry)
        { ack = ROgeometry->CheckROVolume(aStep,ROhis); }
        if(ack)
        { ack = ProcessHits(aStep,ROhis); }
        return ack;
      }
      inline G4int GetNumberOfCollections()
      { return collectionName.entries(); };
      inline G4String GetCollectionName(G4int id)
      { return collectionName[id]; };
      inline void SetVerboseLevel(G4int vl)
      { verboseLevel = vl; };
      inline void Activate(G4bool activeFlag)
      { active = activeFlag; };
      inline G4bool isActive()
      { return active; };
      inline G4String GetName()
      { return SensitiveDetectorName; };
      inline G4String GetPathName()
      { return thePathName; };
      inline G4String GetFullPathName()
      { return fullPathName; };
      inline void SetROgeometry(G4VReadOutGeometry*value)
      { ROgeometry = value; };
      inline G4VReadOutGeometry* GetROgeometry()
      { return ROgeometry; };
};




#endif

