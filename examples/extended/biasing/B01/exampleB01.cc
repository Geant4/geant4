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
//
//
// $Id: exampleB01.cc,v 1.13 2002-09-17 13:59:14 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB01
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "CLHEP/Random/Random.h"
#include "B01SimulationFactory.hh"
#include "B01VSimulation.hh"

void usage(const G4String &simnames){
  G4cout << "exmpleB01: usage: " 
	 << "\n"
	 << "           exmpleB01 <simulationname> <numberofevents>"
	 << "\n"
	 << "           simulationname: \n" << simnames 
	 << "\n" << G4endl;
}

int main(int argc, char **argv)
{  
  G4long myseed = 345354;
  HepRandom::setTheSeed(myseed);



  // construct the simulation factory
  B01SimulationFactory simfac;

  G4int nevents;
  // check the argument
  if (argc != 3) {
    usage(simfac.GetSimulationNames());
    return 0;
  } 
  else {
    if (!simfac.SimulationExists(argv[1])) {
       usage(simfac.GetSimulationNames());
      return 0;
    }
    nevents = atoi(argv[2]);
  }

  // construct 
  B01VSimulation *sim = simfac.Create(argv[1]);
  if (sim) {
    sim->Construct();
    sim->Run(nevents);
    sim->PostRun(&G4cout);
    delete sim;
  }
  else {
    G4cout << "No simulation constructed"  << G4endl;
  }
  return 0;
}
