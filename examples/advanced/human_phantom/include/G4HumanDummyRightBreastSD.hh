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
//    ********************************
//    *                              *
//    *       G4HumanDummySD.hh       *
//    *                              *
//    ********************************
//
// Authors: S. Guatelli , M. G. Pia, INFN Genova and F. Ambroglini INFN Perugia, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
// Dummy sensitive used only to flag sensitivity in cells of RO geometry.

#ifndef G4HumanDummyRightBreastSD_h
#define G4HumanDummyRightBreastSD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class G4HumanDummyRightBreastSD : public G4VSensitiveDetector
{
public:
  G4HumanDummyRightBreastSD();
  ~G4HumanDummyRightBreastSD() {};
  
  void Initialize(G4HCofThisEvent* ) {};
  G4bool ProcessHits(G4Step* ,G4TouchableHistory*) {return false;}
  void EndOfEvent(G4HCofThisEvent*) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};
};

G4HumanDummyRightBreastSD::G4HumanDummyRightBreastSD() : G4VSensitiveDetector("dummyRightBreastSD")
{}
#endif
