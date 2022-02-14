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
// Class G4VScoreHistFiller
//
// Class description:
//  Abstract base class that defines functions to fill histogram and
//  profile plot. This class avoids the direct dependency to G4Analysis
//  category. Template class G4TScoreHistFiller should be used by the
//  user to specify the type of the file format.
//
//  Created: M. Asai (Sept. 2020)
//

#ifndef G4VScoreHistFiller_h
#define G4VScoreHistFiller_h 1

#include "G4Threading.hh"
#include "globals.hh"

// class description:
//
// This class implements the interface for filling histograms from scorer
// type vith Geant4 analysis tools.

class G4VScoreHistFiller
{
 public:
  virtual ~G4VScoreHistFiller();

  // static method
  static G4VScoreHistFiller* Instance();

  // methods
  virtual void FillH1(G4int id, G4double value, G4double weight = 1.0) = 0;
  virtual void FillH2(G4int id, G4double xvalue, G4double yvalue,
                      G4double weight = 1.0)                           = 0;
  virtual void FillH3(G4int id, G4double xvalue, G4double yvalue,
                      G4double zvalue, G4double weight = 1.0)          = 0;
  virtual void FillP1(G4int id, G4double xvalue, G4double yvalue,
                      G4double weight = 1.0)                           = 0;
  virtual void FillP2(G4int id, G4double xvalue, G4double yvalue,
                      G4double zvalue, G4double weight = 1.0)          = 0;

  virtual G4bool CheckH1(G4int id) = 0;
  virtual G4bool CheckH2(G4int id) = 0;
  virtual G4bool CheckH3(G4int id) = 0;
  virtual G4bool CheckP1(G4int id) = 0;
  virtual G4bool CheckP2(G4int id) = 0;

 protected:
  G4VScoreHistFiller();
  virtual G4VScoreHistFiller* CreateInstance() const = 0;

  // static data members
  static G4VScoreHistFiller* fgMasterInstance;
  static G4ThreadLocal G4VScoreHistFiller* fgInstance;
};

#endif
