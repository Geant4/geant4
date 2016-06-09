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
//   HepMCG4PythiaInterface.cc
//   $Id: HepMCG4PythiaInterface.cc,v 1.7 2010-05-24 05:29:44 kmura Exp $
//
// ====================================================================

#ifdef G4LIB_USE_PYTHIA

#include "HepMCG4PythiaInterface.hh"
#include "HepMCG4PythiaMessenger.hh"

#include "HepMC/GenEvent.h" 
#include "HepMC/PythiaWrapper6_4.h" 

// additional pythia calls
#define pygive pygive_
#define pyrget pyrget_
#define pyrset pyrset_

extern "C" {
  void pygive(const char*, int);
  void pyrget(int*, int*);
  void pyrset(int*, int*);
}

void call_pygive(G4String s) { pygive(s.c_str(), s.length()); }
void call_pyrget(int a, int b) { pyrget(&a, &b); }
void call_pyrset(int a, int b) { pyrset(&a, &b); }

////////////////////////////////////////////////
HepMCG4PythiaInterface::HepMCG4PythiaInterface()
  : verbose(0), mpylist(0)
////////////////////////////////////////////////
{
#ifdef NEED_INITPYDATA  
  initpydata();
  // Some platforms may require the initialization of pythia PYDATA block 
  // data as external - if you get pythia initialization errors try 
  // commenting in/out the below call to initpydata().
#endif

  messenger= new HepMCG4PythiaMessenger(this);
}

/////////////////////////////////////////////////
HepMCG4PythiaInterface::~HepMCG4PythiaInterface()
/////////////////////////////////////////////////
{
  delete messenger;
}

/////////////////////////////////////////////////////
void HepMCG4PythiaInterface::CallPygive(G4String par)
/////////////////////////////////////////////////////
{
  call_pygive(par);
}

//////////////////////////////////////////////////////////////////////
void HepMCG4PythiaInterface::CallPyinit(G4String frame, G4String beam, 
					G4String target, G4double win)
//////////////////////////////////////////////////////////////////////
{
  call_pyinit(frame.c_str(), beam.c_str(), target.c_str(), win);
}

////////////////////////////////////////////////////
void HepMCG4PythiaInterface::CallPystat(G4int istat)
////////////////////////////////////////////////////
{
  call_pystat(istat);
}

///////////////////////////////////////////////////////
void HepMCG4PythiaInterface::SetRandomSeed(G4int iseed)
///////////////////////////////////////////////////////
{
  pydatr.mrpy[1-1]= iseed;
}

//////////////////////////////////////////////////////////////
void HepMCG4PythiaInterface::CallPyrget(G4int lun, G4int move)
//////////////////////////////////////////////////////////////
{
  call_pyrget(lun, move);
}

//////////////////////////////////////////////////////////////
void HepMCG4PythiaInterface::CallPyrset(G4int lun, G4int move)
//////////////////////////////////////////////////////////////
{
  call_pyrset(lun, move);
}

//////////////////////////////////////////////////////////////////////////
void HepMCG4PythiaInterface::PrintRandomStatus(std::ostream& ostr) const
//////////////////////////////////////////////////////////////////////////
{
  ostr << "# Pythia random numbers status" << G4endl;
  for (G4int j=0; j<6; j++) {
    ostr << "pydatr.mrpy[" << j << "]= " << pydatr.mrpy[j] << G4endl;
  }
  for (G4int k=0; k<100; k++) {
    ostr << "pydatr.rrpy[" << k << "]= " << pydatr.rrpy[k] << G4endl;
  }
}

////////////////////////////////////////////////
void HepMCG4PythiaInterface::SetUserParameters()
////////////////////////////////////////////////
{
  G4cout << "set user parameters of PYTHIA common." << G4endl
	 << "nothing to be done in default."
	 << G4endl;
}

/////////////////////////////////////////////////////////////
HepMC::GenEvent* HepMCG4PythiaInterface::GenerateHepMCEvent()
/////////////////////////////////////////////////////////////
{
  static G4int nevent= 0; // event counter

  call_pyevnt(); // generate one event with Pythia
  if(mpylist >=1 && mpylist<= 3) call_pylist(mpylist);
  
  call_pyhepc(1);

  HepMC::GenEvent* evt= hepevtio.read_next_event();
  evt-> set_event_number(nevent++);
  if(verbose>0) evt-> print();

  return evt;
}

//////////////////////////////////////////
void HepMCG4PythiaInterface::Print() const
//////////////////////////////////////////
{
  G4cout << "PythiaInterface::Print()..." << G4endl;
}

#endif
