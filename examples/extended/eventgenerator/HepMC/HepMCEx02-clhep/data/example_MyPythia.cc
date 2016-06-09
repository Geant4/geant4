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
// $Id: example_MyPythia.cc,v 1.5 2006/06/29 17:12:32 gunter Exp $
// ====================================================================
//   example_MyPythia.cc
//
//   sample event generation using Pythia and HepMC(CLHEP)
// ====================================================================

#include "CLHEP/HepMC/CBhepevt.h"
#include "CLHEP/HepMC/WriteHepMC.h"
#include "CLHEP/HepMC/ReadHepMC.h"
#include "CLHEP/HepMC/include/PythiaWrapper6_2.h"

#include <string>
#include <fstream>

// global parameters
hepev2_t hepev2_;
hepev3_t hepev3_;
hepev4_t hepev4_;
hepev5_t hepev5_;

int main() { 
  // initialize Pythia...
  // select W+gamma process (process number 20) 
  pysubs.msel= 0;
  pysubs.msub[20-1]= 1;
  // set random number seed (mandatory!)
  pydatr.mrpy[0]= 55122 ;
  // tell Pythia not to write multiple copies of particles in event record.
  pypars.mstp[128-1]= 2;
  // example of setting a Pythia parameter: set the top mass 
  pydat2.pmas[1-1][6-1]= 175;  
  // call pythia initialization
  call_pyinit( "CMS", "p", "p", 14000. );
  
  // prepare HepMC-HEPEVT interface...
  HepMC::CBhepevt* hepevtio= new HepMC::CBhepevt;
  
  std::ofstream to;
  std::string fname= "example_MyPythia.dat";
  to.open(fname.c_str(), std::ios::out);

  HepMC::writeComment(to, "generate sample pythia events...");
  for ( int i = 1; i <= 10; i++ ) {
    if ( i%50==1 ) std::cout << "Processing Event Number " 
			     << i << std::endl;
    call_pyevnt();  // generate one event with Pythia
    // pythia pyhepc routine converts common PYJETS in common HEPEVT
    call_pyhepc( 1 );
    
    HepMC::GenEvent* evt= new HepMC::GenEvent();
    hepevtio-> toGenEvent(evt, false);
    hepevtio-> clean(); // clear HEPEVT
    
    // add some information to the event
    evt-> set_event_number(i);
    evt-> set_signal_process_id(20);
    
    // write the event out to the ascii file
    to << evt;
    
    // we also need to delete the created event from memory
    delete evt;
  }

  // write out some information from Pythia to the screen
  call_pystat(1);
  to.close();
  delete hepevtio;
  
  return 0;
}

