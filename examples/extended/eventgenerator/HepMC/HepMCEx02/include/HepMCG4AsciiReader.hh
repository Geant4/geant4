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
//   HepMCG4AsciiReader.hh
//   $Id: HepMCG4AsciiReader.hh,v 1.1 2002-05-28 14:10:53 murakami Exp $
//
// ====================================================================
#ifndef HEPMC_G4_ASCII_READER_H
#define HEPMC_G4_ASCII_READER_H

#include "HepMCG4Interface.hh"
#include "HepMC/IO_Ascii.h"

class HepMCG4AsciiReaderMessenger;

class HepMCG4AsciiReader : public HepMCG4Interface {
protected:
  G4String filename;
  HepMC::IO_Ascii* asciiInput;

  G4int verbose;
  HepMCG4AsciiReaderMessenger* messenger;

  virtual HepMC::GenEvent* GenerateHepMCEvent();

public:
  HepMCG4AsciiReader();
  ~HepMCG4AsciiReader();

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

inline void HepMCG4AsciiReader::SetFileName(G4String name)
{
  filename= name;
}

inline G4String HepMCG4AsciiReader::GetFileName() const
{
  return filename;
}

inline void HepMCG4AsciiReader::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int HepMCG4AsciiReader::GetVerboseLevel() const
{
  return verbose;
}

#endif
