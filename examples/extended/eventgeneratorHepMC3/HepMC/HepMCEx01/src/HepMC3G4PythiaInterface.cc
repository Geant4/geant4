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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMC3G4PythiaInterface.cc
/// \brief Implementation of the HepMC3G4PythiaInterface class
//
//

#ifdef G4LIB_USE_PYTHIA

#include "HepMC3G4PythiaInterface.hh"
#include "HepMC3G4PythiaMessenger.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "PythiaWrapper6_4.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// additional pythia calls
#define pygive pygive_
#define pyrget pyrget_
#define pyrset pyrset_

extern "C" {
  void pygive(const char*, int);
  void pyrget(int*, int*);
  void pyrset(int*, int*);
  struct HEPEVT hepevt_;
}

void call_pygive(G4String s) { pygive(s.c_str(), s.length()); }
void call_pyrget(int a, int b) { pyrget(&a, &b); }
void call_pyrset(int a, int b) { pyrset(&a, &b); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3G4PythiaInterface::HepMC3G4PythiaInterface()
  : verbose(0), mpylist(0)
{
#ifdef NEED_INITPYDATA
  initpydata();
  // Some platforms may require the initialization of pythia PYDATA block
  // data as external - if you get pythia initialization errors try
  // commenting in/out the below call to initpydata().
#endif
   fGenRunInfo=std::make_shared<HepMC3::GenRunInfo>();
   std::vector<std::string> weight_names={"Default"};
   fGenRunInfo->set_weight_names(weight_names);
   HepMC3::HEPEVT_Wrapper::set_hepevt_address((char*)(&hepevt_));
   
  messenger= new HepMC3G4PythiaMessenger(this);
printf("OK+++++++++");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3G4PythiaInterface::~HepMC3G4PythiaInterface()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::CallPygive(G4String par)
{
  call_pygive(par);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::CallPyinit(G4String frame, G4String beam,
                                        G4String target, G4double win)
{
  call_pyinit(frame.c_str(), beam.c_str(), target.c_str(), win);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::CallPystat(G4int istat)
{
  call_pystat(istat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::SetRandomSeed(G4int iseed)
{
  pydatr.mrpy[1-1]= iseed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::CallPyrget(G4int lun, G4int move)
{
  call_pyrget(lun, move);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::CallPyrset(G4int lun, G4int move)
{
  call_pyrset(lun, move);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::PrintRandomStatus(std::ostream& ostr) const
{
  ostr << "# Pythia random numbers status" << G4endl;
  for (G4int j=0; j<6; j++) {
    ostr << "pydatr.mrpy[" << j << "]= " << pydatr.mrpy[j] << G4endl;
  }
  for (G4int k=0; k<100; k++) {
    ostr << "pydatr.rrpy[" << k << "]= " << pydatr.rrpy[k] << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::SetUserParameters()
{
  G4cout << "set user parameters of PYTHIA common." << G4endl
         << "nothing to be done in default."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3::GenEvent* HepMC3G4PythiaInterface::GenerateHepMCEvent()
{
  static G4int nevent= 0; // event counter

  call_pyevnt(); // generate one event with Pythia
  if(mpylist >=1 && mpylist<= 3) call_pylist(mpylist);

  call_pyhepc(1);
  
  

  
         HepMC3::GenEvent* evt=new HepMC3::GenEvent(HepMC3::Units::GEV,HepMC3::Units::MM);
        for( int i=1; i<=HepMC3::HEPEVT_Wrapper::number_entries(); i++ )
            if (HepMC3::hepevtptr->jmohep[i-1][1]<HepMC3::hepevtptr->jmohep[i-1][0])  HepMC3::hepevtptr->jmohep[i-1][1]=HepMC3::hepevtptr->jmohep[i-1][0];
        HepMC3::HEPEVT_Wrapper::HEPEVT_to_GenEvent(evt);
        evt->set_run_info(fGenRunInfo);
        evt->weights()=std::vector<double>(fGenRunInfo->weight_names().size(),1.0);
  evt-> set_event_number(nevent++);
  if(verbose>0) HepMC3::Print::content(*evt);

  return evt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4PythiaInterface::Print() const
{
  G4cout << "PythiaInterface::Print()..." << G4endl;
}

#endif
