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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include "G4HumanPhantomEnergyDeposit.hh"
#ifdef G4ANALYSIS_USE
#include "G4HumanPhantomAnalysisManager.hh"
#endif
#include "globals.hh"
#include "G4UnitsTable.hh"

G4HumanPhantomEnergyDeposit::G4HumanPhantomEnergyDeposit()
{  
}

G4HumanPhantomEnergyDeposit::~G4HumanPhantomEnergyDeposit()
{
}

void G4HumanPhantomEnergyDeposit::Fill(G4String bodypartName, 
			G4double energyDeposit)
{
  energyTotal[bodypartName] += energyDeposit;
}

void G4HumanPhantomEnergyDeposit::TotalEnergyDeposit()
{
  
  G4double totalBody = energyTotal["HeadVolume"]+energyTotal["TrunkVolume"]+energyTotal["LegsVolume"];
  G4cout << "Energy in total body = " << G4BestUnit(totalBody,"Energy") << G4endl;

  std::map<std::string,G4double>::iterator i = energyTotal.begin();
  std::map<std::string,G4double>::iterator end = energyTotal.end();

  G4double k = 0.;;

  while(i!=end)
    {

      G4String bodypart = i->first;
      G4double energyDep = i->second;

      G4cout << "Energy Total in " <<bodypart << " = "  
	     << G4BestUnit(energyDep,"Energy") 
	     << G4endl;

#ifdef G4ANALYSIS_USE
  G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
  if (energyDep != 0.)analysis -> bodypartEnergyDep(k,energyDep/MeV);
#endif
      i++;
      k++;
    }
  
}
