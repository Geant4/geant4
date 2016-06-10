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
// $Id: G4ASCIITreeMessenger.hh 66870 2013-01-14 23:38:59Z adotti $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy to standard output as
//   ASCII stream.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#ifndef G4ASCIITREEMESSENGER_HH
#define G4ASCIITREEMESSENGER_HH

#include "G4UImessenger.hh"

#include <vector>

class G4ASCIITree;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class G4ASCIITreeMessenger: public G4UImessenger {
public:
  G4ASCIITreeMessenger(G4ASCIITree*);
  virtual ~G4ASCIITreeMessenger();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
  static std::vector<G4String> fVerbosityGuidance;
private:
  G4ASCIITree* fpASCIITree;
  G4UIcommand* fpDirectory;
  G4UIcommand* fpDirectorySet;
  G4UIcmdWithAnInteger* fpCommandVerbose;
  G4UIcmdWithAString* fpCommandSetOutFile;
};

#endif
