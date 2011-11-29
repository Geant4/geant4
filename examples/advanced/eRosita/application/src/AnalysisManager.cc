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
#include "AnalysisManager.hh"
#include <sstream>

AnalysisManager* AnalysisManager::instance = 0;


AnalysisManager* AnalysisManager::Instance() 
{
 
  // A new instance of AnalysisManager is created if it does not exist:
  if (instance == 0) 
    {
      instance = new AnalysisManager();
    }
  
  // The instance of AnalysisManager is returned:
  return instance;
}


void AnalysisManager::Destroy() 
{

  // The AnalysisManager instance is deleted if it exists:
  if (instance != 0) 
    {
      delete instance;
      instance = 0;
    }
}


AnalysisManager::AnalysisManager() 
{
  outFile.open("TrackerPhotonEnergy.out");
//   outFileT.open("eTot.out"); 
}


AnalysisManager::~AnalysisManager() 
{
  outFile.close();
  outFileT.close();
}


void AnalysisManager::Score(G4double eDep) 
{
//   outFile << id << " "
// 	  << eDep << " "
// 	  << x << " "
// 	  << y << " "
// 	  << z << " "
// 	  << std::endl;
  outFile << eDep << std::endl;
}

// void AnalysisManager::ScoreTot(G4double eTot) 
// {
//   outFileT << eTot << std::endl;
// }
