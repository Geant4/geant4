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
//   G4HepMCPythiaInterface.cc
//   $Id: G4HepMCPythiaInterface.cc,v 1.1 2002-04-29 20:44:52 asaim Exp $
//
// ====================================================================
#include "G4HepMCPythiaInterface.hh"
#include "G4HepMCPythiaMessenger.hh"

#include "HepMC/GenEvent.h" 
#include "HepMC/PythiaWrapper6_152.h" 

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
G4HepMCPythiaInterface::G4HepMCPythiaInterface()
  : verbose(0), mpylist(0)
////////////////////////////////////////////////
{
#ifdef NEED_INITPYDATA  
  initpydata();
  // Some platforms may require the initialization of pythia PYDATA block 
  // data as external - if you get pythia initialization errors try 
  // commenting in/out the below call to initpydata().
#endif

  messenger= new G4HepMCPythiaMessenger(this);
}

/////////////////////////////////////////////////
G4HepMCPythiaInterface::~G4HepMCPythiaInterface()
/////////////////////////////////////////////////
{
  delete messenger;
}

/////////////////////////////////////////////////////
void G4HepMCPythiaInterface::CallPygive(G4String par)
/////////////////////////////////////////////////////
{
  call_pygive(par);
}

//////////////////////////////////////////////////////////////////////
void G4HepMCPythiaInterface::CallPyinit(G4String frame, G4String beam, 
					G4String target, G4double win)
//////////////////////////////////////////////////////////////////////
{
  call_pyinit(frame.c_str(), beam.c_str(), target.c_str(), win);
}

////////////////////////////////////////////////////
void G4HepMCPythiaInterface::CallPystat(G4int istat)
////////////////////////////////////////////////////
{
  call_pystat(istat);
}

///////////////////////////////////////////////////////
void G4HepMCPythiaInterface::SetRandomSeed(G4int iseed)
///////////////////////////////////////////////////////
{
  pydatr.mrpy[1-1]= iseed;
}

//////////////////////////////////////////////////////////////
void G4HepMCPythiaInterface::CallPyrget(G4int lun, G4int move)
//////////////////////////////////////////////////////////////
{
  call_pyrget(lun, move);
}

//////////////////////////////////////////////////////////////
void G4HepMCPythiaInterface::CallPyrset(G4int lun, G4int move)
//////////////////////////////////////////////////////////////
{
  call_pyrset(lun, move);
}

//////////////////////////////////////////////////////////////////////////
void G4HepMCPythiaInterface::PrintRandomStatus(G4std::ostream& ostr) const
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
void G4HepMCPythiaInterface::SetUserParameters()
////////////////////////////////////////////////
{
  G4cout << "set user parameters of PYTHIA common." << G4endl
	 << "nothing to be done in default."
	 << G4endl;
}

/////////////////////////////////////////////////////////////
HepMC::GenEvent* G4HepMCPythiaInterface::GenerateHepMCEvent()
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
void G4HepMCPythiaInterface::Print() const
//////////////////////////////////////////
{
  G4cout << "PythiaInterface::Print()..." << G4endl;
}

