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
//   G4HepMCAsciiReader.hh
//   $Id: G4HepMCAsciiReader.hh,v 1.3 2006/06/29 17:07:58 gunter Exp $
//
// ====================================================================
#ifndef G4_HEPMC_ASCII_READER_H
#define G4_HEPMC_ASCII_READER_H

#include <fstream>
#include "G4HepMCInterface.hh"
#include <CLHEP/HepMC/ReadHepMC.h>

class G4HepMCAsciiReaderMessenger;

class G4HepMCAsciiReader : public G4HepMCInterface
{
protected:
  G4String filename;
  std::ifstream from;

  G4int verbose;
  G4HepMCAsciiReaderMessenger* messenger;

  virtual HepMC::GenEvent* GenerateHepMCEvent();

public:
  G4HepMCAsciiReader();
  ~G4HepMCAsciiReader();

  // set/get methods
  void SetFileName(G4String name);
  G4String GetFileName() const;

  void SetVerboseLevel(G4int i);
  G4int GetVerboseLevel() const; 

  // methods...
  void Initialize();
};

// ====================================================================
// inline functions
// ====================================================================

inline void G4HepMCAsciiReader::SetFileName(G4String name)
{
  filename= name;
}

inline G4String G4HepMCAsciiReader::GetFileName() const
{
  return filename;
}

inline void G4HepMCAsciiReader::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int G4HepMCAsciiReader::GetVerboseLevel() const
{
  return verbose;
}

#endif
