// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSensitiveDetector.hh,v 1.2.2.1.2.1 1999/12/07 20:47:43 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4VSensitiveDetector_h
#define G4VSensitiveDetector_h 1

#include "G4VHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "g4rw/tvordvec.h"

// class description:
//
//  This is the abstract base class of the sensitive detector. The user's
// sensitive detector which generates hits must be derived from this
// class.
//  In the derived class constructor, name(s) of hits collection(s) which
// are made by the sensitive detector must be set to "collectionName" string
// vector.

class G4VSensitiveDetector 
{

  public: // with description
      G4VSensitiveDetector(G4String name);
      G4VSensitiveDetector(const G4VSensitiveDetector &right);
      // Constructors. The user's concrete class must use one of these constructors
      // by the constructor initializer of the derived class. The name of
      // the sensitive detector must be unique.

  public:
      virtual ~G4VSensitiveDetector();

      const G4VSensitiveDetector & operator=(const G4VSensitiveDetector &right);

      int operator==(const G4VSensitiveDetector &right) const;
      int operator!=(const G4VSensitiveDetector &right) const;

  public: // with description
      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);
      //  These two methods are invoked at the begining and at the end of each
      // event. The hits collection(s) created by this sensitive detector must
      // be set to the G4HCofThisEvent object at one of these two methods.
      virtual void clear();
      //  This method is invoked if the event abortion is occured. Hits collections
      // created but not beibg set to G4HCofThisEvent at the event should be deleted.
      // Collection(s) which have already set to G4HCofThisEvent will be deleted 
      // automatically.

  public:
      virtual void DrawAll();
      virtual void PrintAll();

  protected: // with description
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) = 0;
      //  The user MUST implement this method for generating hit(s) using the 
      // information of G4Step object. Note that the volume and the position
      // information is kept in PreStepPoint of G4Step.
      //  Be aware that this method is a protected method and it sill be invoked 
      // by Hit() method of Base class after Readout geometry associated to the
      // sensitive detector is handled.
      //  "ROhist" will be given only is a Readout geometry is defined to this
      // sensitive detector. The G4TouchableHistory object of the tracking geometry
      // is stored in the PreStepPoint object of G4Step.
      G4int GetCollectionID(G4int i);
      //  This is a utility method which returns the hits collection ID of the
      // "i"-th collection. "i" is the order (starting with zero) of the collection
      // whose name is stored to the collectionName protected vector.
      G4RWTValOrderedVector<G4String> collectionName;
      //  This protected name vector must be filled at the constructor of the user's
      // concrete class for registering the name(s) of hits collection(s) being
      // created by this particular sensitive detector.

  protected:
      G4String SensitiveDetectorName; // detector name
      G4String thePathName;           // directory path
      G4String fullPathName;          // path + detector name
      G4int verboseLevel;
      G4bool active;

  private:
      G4VReadOutGeometry * ROgeometry;

  public: // with description
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
      //  This is the public method invoked by G4SteppingManager for generating
      // hit(s). The actual user's implementation for generating hit(s) must be
      // implemented in GenerateHits() virtual protected method. This method
      // MUST NOT be overrided.
      inline void SetROgeometry(G4VReadOutGeometry*value)
      { ROgeometry = value; }
      //  Register the Readout geometry.

  public:
      inline G4int GetNumberOfCollections()
      { return collectionName.entries(); }
      inline G4String GetCollectionName(G4int id)
      { return collectionName[id]; }
      inline void SetVerboseLevel(G4int vl)
      { verboseLevel = vl; }
      inline void Activate(G4bool activeFlag)
      { active = activeFlag; }
      inline G4bool isActive()
      { return active; }
      inline G4String GetName()
      { return SensitiveDetectorName; }
      inline G4String GetPathName()
      { return thePathName; }
      inline G4String GetFullPathName()
      { return fullPathName; }
      inline G4VReadOutGeometry* GetROgeometry()
      { return ROgeometry; }
};




#endif

