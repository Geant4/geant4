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

#ifndef G4TScoreHistFiller_h
#define G4TScoreHistFiller_h 1

#include "G4VScoreHistFiller.hh"
#include "globals.hh"

#include <memory>

// class description:
//
// This class implements filling histogram
// In order to avoid introducing dependency on the analysis category,
// the analysis manager type is defined via template.
//
// Created : M. Asai (Sept. 2020)
//

template <typename T>
class G4TScoreHistFiller : public G4VScoreHistFiller
{
 public:
  G4TScoreHistFiller();
  virtual ~G4TScoreHistFiller();

  // methods
  virtual void FillH1(G4int id, G4double value, G4double weight = 1.0);
  virtual void FillH2(G4int id, G4double xvalue, G4double yvalue,
                      G4double weight = 1.0);
  virtual void FillH3(G4int id, G4double xvalue, G4double yvalue,
                      G4double zvalue, G4double weight = 1.0);
  virtual void FillP1(G4int id, G4double xvalue, G4double yvalue,
                      G4double weight = 1.0);
  virtual void FillP2(G4int id, G4double xvalue, G4double yvalue,
                      G4double zvalue, G4double weight = 1.0);

  virtual G4bool CheckH1(G4int id);
  virtual G4bool CheckH2(G4int id);
  virtual G4bool CheckH3(G4int id);
  virtual G4bool CheckP1(G4int id);
  virtual G4bool CheckP2(G4int id);

  void SetVerboseLevel(G4int value);
  G4int GetVerboseLevel() const { return fVerboseLevel; }

 protected:
  // methods
  virtual G4VScoreHistFiller* CreateInstance() const;

 private:
  // methods
  void CreateAnalysisManager();

  // data members
  T* fAnalysisManager        = nullptr;
  G4int fVerboseLevel        = 0;
  G4bool fIsInitialized      = false;
};

#include "G4TScoreHistFiller.icc"

#endif
