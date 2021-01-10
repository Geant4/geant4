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
/// \file eventgenerator/HepMC/HepMCEx01/include/HepMC3G4AsciiReader.hh
/// \brief Definition of the HepMC3G4AsciiReader class
//
//

#ifndef HEPMC3_G4_ASCII_READER_H
#define HEPMC3_G4_ASCII_READER_H

#include "HepMC3G4Interface.hh"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/Print.h"

class HepMC3G4AsciiReaderMessenger;

class HepMC3G4AsciiReader : public HepMC3G4Interface {
protected:
  G4String filename;
  HepMC3::ReaderAsciiHepMC2* asciiInput;

  G4int verbose;
  HepMC3G4AsciiReaderMessenger* messenger;

  virtual HepMC3::GenEvent* GenerateHepMCEvent();

public:
  HepMC3G4AsciiReader();
  ~HepMC3G4AsciiReader();

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

inline void HepMC3G4AsciiReader::SetFileName(G4String name)
{
  filename= name;
}

inline G4String HepMC3G4AsciiReader::GetFileName() const
{
  return filename;
}

inline void HepMC3G4AsciiReader::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int HepMC3G4AsciiReader::GetVerboseLevel() const
{
  return verbose;
}

#endif
