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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
/// \file RunInitObserver.cc
/// \brief Implementation of the RunInitObserver class

#include "RunInitObserver.hh"

G4ThreadLocal RunInitManager* RunInitManager::fgInstance(0);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager::RunInitManager()
{
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager::~RunInitManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitObserver::RunInitObserver()
{
  RunInitManager::Instance()->Insert(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitObserver::~RunInitObserver()
{
  // TODO Auto-generated destructor stub
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager* RunInitManager::Instance()
{
  if(fgInstance == 0) new RunInitManager();
  return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
RunInitManager::Initialize()
{
  std::vector<RunInitObserver*>::iterator it = fObservers.begin();
  std::vector<RunInitObserver*>::iterator end = fObservers.end();
  for(; it != end ; it++) (*it)->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
RunInitManager::Insert(RunInitObserver* observer)
{
  fObservers.push_back(observer);
}
