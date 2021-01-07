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
/// \file eventgenerator/HepMC/HepMCEx01/include/HepMC3G4PythiaInterface.hh
/// \brief Definition of the HepMC3G4PythiaInterface class
//
//

#ifndef HEPMC3_G4_PYTHIA_INTERFACE_H
#define HEPMC3_G4_PYTHIA_INTERFACE_H

#include "HepMC3G4Interface.hh"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/HEPEVT_Wrapper.h"

class HepMC3G4PythiaMessenger;

/// A generic interface class with Pythia event generator via HepMC.

class HepMC3G4PythiaInterface : public HepMC3G4Interface {
protected:
  G4int verbose;
  G4int mpylist;
  std::shared_ptr<HepMC3::GenRunInfo> fGenRunInfo;

  HepMC3G4PythiaMessenger* messenger;

  // In default, this is automatic conversion, Pythia->HEPEVT->HepMC,
  // by HepMC utilities.
  virtual HepMC3::GenEvent* GenerateHepMCEvent();

public:
  HepMC3G4PythiaInterface();
  ~HepMC3G4PythiaInterface();

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
  void PrintRandomStatus(std::ostream& ostr=G4cout) const;

  // setup user parameters (empty in default).
  // Implement your parameters in a delived class if you want.
  virtual void SetUserParameters();

  virtual void Print() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void HepMC3G4PythiaInterface::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int HepMC3G4PythiaInterface::GetVerboseLevel() const
{
  return verbose;
}

inline void HepMC3G4PythiaInterface::SetPylist(G4int i)
{
  mpylist= i;
}

inline G4int HepMC3G4PythiaInterface::GetPylist() const
{
  return mpylist;
}

#endif
