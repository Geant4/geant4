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
// $Id: G4SDStructure.hh 74048 2013-09-20 09:34:57Z gcosmo $
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
      void RemoveSD(G4VSensitiveDetector*);

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

