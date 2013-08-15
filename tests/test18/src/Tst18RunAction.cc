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
// $Id$
//

#include "Tst18RunActionMessenger.hh"
#include "Tst18RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include <fstream>
#include <iomanip>
#include <vector>

Tst18RunAction::Tst18RunAction()
 : fileName("rdmex2.log")
{
  runMessenger = new Tst18RunActionMessenger(this);
}


Tst18RunAction::~Tst18RunAction()
{
  delete runMessenger;
}


void Tst18RunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  G4int RunN = aRun->GetRunID();
  if ( RunN % 1000 == 0 ) 
    G4cout << "### Run : " << RunN << G4endl;

  Particles.clear();
  Energies.clear();
  Weights.clear();
  Times.clear();
  EventNumber = 0;
}


void Tst18RunAction::EndOfRunAction(const G4Run* )
{
  
  for (size_t i=0; i<Particles.size();i++) {

    G4cout
      << std::setiosflags(std::ios::fixed)
      << std::setprecision(3)
      << std::setiosflags(std::ios::right)
      << std::setw(12)
      << Energies[i]
      << std::setw(12)<<std::setprecision(4) 
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::right)
      << Weights[i]
      << std::setw(12)<<std::setprecision(4)
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::right)
      << Times[i] << "     "
      << Particles[i]
      << G4endl ;    
  }
  G4cout<< G4endl;

  if(Particles.size()/EventNumber == 1) {
    G4cout << " Abort " << G4endl;
    G4Exception("Tst18RunAction::EndOfRunAction", "Tst18_00",
		FatalException, "No Decay Products Found");
  }

}

void Tst18RunAction::FillParticleName(G4String name)
{
  Particles.push_back(name);
}

void Tst18RunAction::FillEnergy(G4double energy)
{
  Energies.push_back(energy);
}

void Tst18RunAction::FillWeight(G4double weight)
{
  Weights.push_back(weight);
}

void Tst18RunAction::FillTime(G4double time)
{
  Times.push_back(time);
}

void Tst18RunAction::FillEventNumber()
{
  EventNumber++;
}

