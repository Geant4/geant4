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
// $Id: test201DetectorMessenger.hh,v 1.4 2006-06-29 21:47:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction Messenger for visualization testing.
// John Allison 25th April 1997

#ifndef test201DetectorMessenger_hh
#define test201DetectorMessenger_hh

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"

class G4UIcommand;
class test201DetectorConstruction;

class test201DetectorMessenger: public G4UImessenger
{
public:
  test201DetectorMessenger (test201DetectorConstruction* test201Det);
  ~test201DetectorMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  test201DetectorConstruction* test201Detector;
  G4UIcommand* fpTest201DetCommandDirectory;
  G4UIcommand* fpDetectorCommand;
};

#endif
