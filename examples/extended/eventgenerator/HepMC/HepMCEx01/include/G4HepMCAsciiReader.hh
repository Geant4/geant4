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

// ====================================================================
//
//   G4HepMCAsciiReader.hh
//   $Id: G4HepMCAsciiReader.hh,v 1.1 2002-04-29 20:43:51 asaim Exp $
//
// ====================================================================
#ifndef G4_HEPMC_ASCII_READER_H
#define G4_HEPMC_ASCII_READER_H

#include "G4HepMCInterface.hh"
#include "HepMC/IO_Ascii.h"

class G4HepMCAsciiReaderMessenger;

class G4HepMCAsciiReader : public G4HepMCInterface {
protected:
  G4String filename;
  HepMC::IO_Ascii* asciiInput;

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
