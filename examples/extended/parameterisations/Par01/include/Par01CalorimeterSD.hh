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
/// \file Par01/include/Par01CalorimeterSD.hh
/// \brief Definition of the Par01CalorimeterSD class
//
//
//

#ifndef Par01CalorimeterSD_h
#define Par01CalorimeterSD_h 1

#include "Par01CalorimeterHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class Par01CalorimeterSD : public G4VSensitiveDetector
{
  
public:
  Par01CalorimeterSD(G4String name, G4int nCells, G4String colName);
  ~Par01CalorimeterSD();
  
  virtual void   Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  virtual void   EndOfEvent(G4HCofThisEvent*HCE);
  virtual void   clear();
  virtual void   DrawAll();
  virtual void   PrintAll();
  
private:
  Par01CalorimeterHitsCollection *fCalCollection;
  
  int* fCellID;
  int  fNumberOfCells;
  int  fHCID;
};


#endif

