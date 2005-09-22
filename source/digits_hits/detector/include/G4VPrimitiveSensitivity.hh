//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPrimitiveSensitivity.hh,v 1.2 2005-09-22 22:21:36 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VPrimitiveSensitivity_h
#define G4VPrimitiveSensitivity_h 1

class G4Step;
class G4HCofThisEvent;
class G4MultiFunctionalDetector;
class G4TouchableHistory;
#include "globals.hh"
#include "G4VSDFilter.hh"

// class description:
//
//  This is the base class of the sensitive detector which owns
// only one hits collection.
//  A concrete class object derived from this base class can be
// used either as a sensitive detector or to be registered to 
// G4MultiFunctionalDetector to define multiple functionalities.

class G4VPrimitiveSensitivity
{
  friend class G4MultiFunctionalDetector;

  public: // with description
      G4VPrimitiveSensitivity(G4String);
      virtual ~G4VPrimitiveSensitivity();

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*)=0;
      // This is the method must be implemented in each concrete class.

      virtual G4int GetIndex(G4Step*);
      // This is a function mapping from copy number(s) to an index of 
      // the hit collection. In the default implementation, just the
      // copy number of the physical volume is taken.

  public: // with description
      G4int GetCollectionID(G4int);
      // This method returns the ID of its hitsCollection. This mehod
      // gives valid value only after it is registered to G4MultiFunctionalDetector
      // and the G4MultiFunctionalDetector is registered to G4SDManager.

      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();
      // These five methods are exactly identical to those in G4VSensitiveDetector.
      // These methods are invoked by G4SDManager through G4MultiFunctionalDetector.

  private:
      G4String primitiveName;
      G4MultiFunctionalDetector* detector;
      G4VSDFilter* filter;
      G4int verboseLevel;

  public: // with description
      // Set/Get methods
      inline void SetMultiFunctionalDetector(G4MultiFunctionalDetector* d)
      { detector = d; }
      inline G4MultiFunctionalDetector* GetMultiFunctionalDetector() const
      { return detector; }
      inline G4String GetName() const
      { return primitiveName; }
      inline void SetFilter(G4VSDFilter* f)
      { filter = f; }
      inline G4VSDFilter* GetFilter() const
      { return filter; }
      inline void SetVerboseLevel(G4int vl)
      { verboseLevel = vl; }

  private:
      inline G4bool HitPrimitive(G4Step*aStep,G4TouchableHistory*ROhis)
      {
        if(filter)
        { if(!(filter->Accept(aStep))) return false; }
        return ProcessHits(aStep,ROhis);
      }
 
};



#endif

