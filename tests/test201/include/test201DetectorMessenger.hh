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
// $Id: test201DetectorMessenger.hh,v 1.3 2001-07-11 10:10:22 gunter Exp $
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
