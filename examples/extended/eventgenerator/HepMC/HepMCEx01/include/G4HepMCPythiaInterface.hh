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
//   G4HepMCPythiaInterface.hh
//   $Id: G4HepMCPythiaInterface.hh,v 1.1 2002-04-29 20:43:53 asaim Exp $
//
//   A generic interface class with Pythia event generator via HepMC.
//
// ====================================================================
#ifndef G4_HEPMC_PYTHIA_INTERFACE_H
#define G4_HEPMC_PYTHIA_INTERFACE_H

#include "G4HepMCInterface.hh"
#include "HepMC/IO_HEPEVT.h"

class G4HepMCPythiaMessenger;

class G4HepMCPythiaInterface : public G4HepMCInterface {
protected:
  G4int verbose;
  G4int mpylist;
  HepMC::IO_HEPEVT hepevtio;

  G4HepMCPythiaMessenger* messenger;

  // In default, this is automatic conversion, Pythia->HEPEVT->HepMC, 
  // by HepMC utilities.
  virtual HepMC::GenEvent* GenerateHepMCEvent();

public:
  G4HepMCPythiaInterface();
  ~G4HepMCPythiaInterface();

  // set/get methods
  void SetVerboseLevel(G4int i);
  G4int GetVerboseLevel() const; 

  void SetPylist(G4int i);
  G4int GetPylist() const; 

  // call pyxxx
  void CallPyinit(G4String frame, G4String beam, G4String target, 
		  G4double win);
  void CallPystat(G4int istat);

  // random numbers operations
  void SetRandomSeed(G4int iseed);
  void CallPygive(G4String par);
  void CallPyrget(G4int lun, G4int move);
  void CallPyrset(G4int lun, G4int move);
  void PrintRandomStatus(G4std::ostream& ostr=G4std::cout) const;

  // setup user parameters (empty in default).
  // Implement your parameters in a delived class if you want.
  virtual void SetUserParameters();

  virtual void Print() const;
};

// ====================================================================
// inline functions
// ====================================================================

inline void G4HepMCPythiaInterface::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int G4HepMCPythiaInterface::GetVerboseLevel() const
{
  return verbose;
}

inline void G4HepMCPythiaInterface::SetPylist(G4int i)
{
  mpylist= i;
}

inline G4int G4HepMCPythiaInterface::GetPylist() const
{
  return mpylist;
}

#endif
