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
// $Id: GammaRayTelDummySD.hh,v 1.5 2003/05/30 15:08:47 flongo Exp $
// GEANT4 tag $Name: geant4-06-00 $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDummySD  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************
//
// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef GammaRayTelDummySD_h
#define GammaRayTelDummySD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;

class GammaRayTelDummySD : public G4VSensitiveDetector
{
public:
  GammaRayTelDummySD(G4String name):G4VSensitiveDetector(name){};
  ~GammaRayTelDummySD() {};
  
  void Initialize(G4HCofThisEvent* ) {};
  G4bool ProcessHits(G4Step*,G4TouchableHistory*) {return false;}
  void EndOfEvent(G4HCofThisEvent*) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};
};

#endif
