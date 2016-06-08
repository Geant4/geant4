// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SDStructure.hh,v 1.2.2.1.2.1 1999/12/07 20:47:41 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef G4SDStructure_h
#define G4SDStructure_h 1

// Globals
#include "globals.hh"
// G4VSensitiveDetector
#include "G4VSensitiveDetector.hh"
// G4RWTPtrOrderedVector
#include "g4rw/tpordvec.h"

class G4HCofThisEvent;

// class description:
//
//  This class is exclusively used by G4SDManager for handling the tree
// structure of the user's sensitive detector names.
//

class G4SDStructure 
{
  public:
      G4SDStructure(G4String aPath);
      ~G4SDStructure();

      int operator==(const G4SDStructure &right) const;

      void AddNewDetector(G4VSensitiveDetector*aSD, G4String treeStructure);
      void Activate(G4String aName, G4bool sensitiveFlag);
      void Initialize(G4HCofThisEvent*HCE);
      void Terminate(G4HCofThisEvent*HCE);
      G4VSensitiveDetector* FindSensitiveDetector(G4String aName);
      G4VSensitiveDetector* GetSD(G4String aName);
      void ListTree();

  private:
      G4SDStructure* FindSubDirectory(G4String subD);
      G4String ExtractDirName(G4String aPath);

  private:
      G4RWTPtrOrderedVector<G4SDStructure> structure;
      G4RWTPtrOrderedVector<G4VSensitiveDetector> detector;
      G4String pathName;
      G4String dirName;
      G4int verboseLevel;

  public:
      inline void SetVerboseLevel(G4int vl) 
      {
        verboseLevel = vl;
        for(int i=0; i<structure.entries(); i++)
        { structure(i)->SetVerboseLevel(vl); }
        for(int j=0; j<detector.entries(); j++)
        { detector(j)->SetVerboseLevel(vl); }
      };

};




#endif

