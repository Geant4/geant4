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
// ====================================================================
//
//   G4HepMCAsciiReaderMessenger.hh
//   $Id: G4HepMCAsciiReaderMessenger.hh,v 1.2 2006/06/29 17:08:01 gunter Exp $
//
// ====================================================================
#ifndef G4_HEPMC_ASCII_READER_MESSENGER_H
#define G4_HEPMC_ASCII_READER_MESSENGER_H

#include "G4UImessenger.hh"

class G4HepMCAsciiReader;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class G4HepMCAsciiReaderMessenger : public G4UImessenger {
private:
  G4HepMCAsciiReader* gen;

  G4UIdirectory* dir;
  G4UIcmdWithAnInteger* verbose;
  G4UIcmdWithAString* open;

public:
  G4HepMCAsciiReaderMessenger(G4HepMCAsciiReader* agen);
  ~G4HepMCAsciiReaderMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
};

#endif
