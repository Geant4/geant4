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
// $Id: MyDetectorMessenger.hh,v 1.3 2001-07-11 10:09:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyDetectorMessenger_h
#define MyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"

class MyDetectorConstruction;

class MyDetectorMessenger: public G4UImessenger
{
public:
  MyDetectorMessenger(MyDetectorConstruction * myDet);
  ~MyDetectorMessenger ();
  void SetNewValue(G4UIcommand * command,G4String newValues);
private:
  MyDetectorConstruction * myDetector;
  G4UIcommand* fpMyDetCommandDirectory;
  G4UIcommand* fpCalMaterialCommand;
};

#endif

