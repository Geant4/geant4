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
// $Id: G4SDStructure.hh,v 1.2 2004/05/03 08:14:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//

#ifndef G4SDStructure_h
#define G4SDStructure_h 1

// Globals
#include "globals.hh"
// G4VSensitiveDetector
#include "G4VSensitiveDetector.hh"
#include <vector>

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

      G4int operator==(const G4SDStructure &right) const;

      void AddNewDetector(G4VSensitiveDetector*aSD, G4String treeStructure);
      void Activate(G4String aName, G4bool sensitiveFlag);
      void Initialize(G4HCofThisEvent*HCE);
      void Terminate(G4HCofThisEvent*HCE);
      G4VSensitiveDetector* FindSensitiveDetector(G4String aName, G4bool warning = true);
      G4VSensitiveDetector* GetSD(G4String aName);
      void ListTree();

  private:
      G4SDStructure* FindSubDirectory(G4String subD);
      G4String ExtractDirName(G4String aPath);

  private:
      std::vector<G4SDStructure*> structure;
      std::vector<G4VSensitiveDetector*> detector;
      G4String pathName;
      G4String dirName;
      G4int verboseLevel;

  public:
      inline void SetVerboseLevel(G4int vl) 
      {
        verboseLevel = vl;
        for(size_t i=0; i<structure.size(); i++)
        { structure[i]->SetVerboseLevel(vl); }
        for(size_t j=0; j<detector.size(); j++)
        { detector[j]->SetVerboseLevel(vl); }
      };

};




#endif

