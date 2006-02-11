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


#include "G4HumanPhantomEnergyDeposit.hh"
#include "G4HumanPhantomAnalysisManager.hh"

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
  totalBody = energyTotal["HeadVolume"]+energyTotal["TrunkVolume"]+energyTotal["LegsVolume"];
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
  if (energyDep != 0.)analysis -> bodypartEnergyDep(k,energyDep);
#endif

      i++;
      k++;
    }

#ifdef G4ANALYSIS_USE
  G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
  analysis->finish();
#endif

}
