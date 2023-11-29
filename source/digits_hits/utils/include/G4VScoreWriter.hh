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
//

#ifndef G4VScoreWriter_h
#define G4VScoreWriter_h 1

#include "globals.hh"
class G4VScoringMesh;

// class description:
//
//  This class represents storing the scored quantity into a file.
//

class G4VScoreWriter
{
 public:
  G4VScoreWriter() = default;
  virtual ~G4VScoreWriter() = default;

 public:
  // store a quantity into a file
  virtual void DumpQuantityToFile(const G4String& psName,
                                  const G4String& fileName,
                                  const G4String& option);
  // store all quantities into a file
  virtual void DumpAllQuantitiesToFile(const G4String& fileName,
                                       const G4String& option);

  // set a socring mesh to retrieve its quantities
  void SetScoringMesh(G4VScoringMesh* sm);
  // set a verbose level
  inline void SetVerboseLevel(G4int vl) { verboseLevel = vl; }
  inline void SetFactor(G4double val = 1.0) { fact = val; }
  inline G4double GetFactor() const { return fact; }

 protected:
  // get an index from (x,y,z)
  G4int GetIndex(G4int x, G4int y, G4int z) const;

 protected:
  G4int fNMeshSegments[3] = {0, 0, 0};  // number of segments of the mesh
  G4VScoringMesh* fScoringMesh = nullptr;
  G4int verboseLevel = 0;
  G4double fact = 1.0;
};

#endif
