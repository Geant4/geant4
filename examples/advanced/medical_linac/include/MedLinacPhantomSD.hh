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
//// $Id: MedLinacPhantomSD.hh,v 1.7 2006/06/29 16:03:55 gunter Exp $
//
//

// Code developed by: M. Piergentili
//
#ifndef MedLinacPhantomSD_h
#define MedLinacPhantomSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MedLinacPhantomHit.hh"
class MedLinacPhantomMessenger;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class MedLinacPhantomSD : public G4VSensitiveDetector
{
public:
  MedLinacPhantomSD(G4String name);
  ~MedLinacPhantomSD();


  void SetPhantomDimension (G4double);
  void SetNumberOfPhantomVoxels (G4int);

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

public:
  G4double GetPhantomDimension()  {return phantomDimension;};
  G4int GetNumberOfPhantomVoxels()  {return numberOfPhantomVoxels;};

  private:
    G4double phantomDimension;
  G4int numberOfPhantomVoxels;
  MedLinacPhantomMessenger* phantomMessenger;
};
#endif


