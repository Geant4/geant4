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
//
// Created by ngoc hoang tran on 03/08/2023.
//

#include "G4ChemEquilibrium.hh"
#include "G4DNAMolecularReactionTable.hh"

G4ChemEquilibrium::G4ChemEquilibrium(const G4int& type, const G4double& time)
  : fEquilibriumDuration(time), fRectionType(type)
{}

void G4ChemEquilibrium::Initialize()
{
  MolType H2O =
    G4MoleculeTable::Instance()->GetConfiguration("H2O");
  MolType H3OpB =
    G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  MolType OHmB =
    G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");
  const auto& reactionList = G4DNAMolecularReactionTable::Instance()->
GetVectorOfReactionData();
  for(const auto& it : reactionList)
  {
    if(it->GetReactionType()==fRectionType)
    {
      if(it->GetReactant1() != H2O
          && it->GetReactant1() != H3OpB
          && it->GetReactant1() != OHmB)
      {
        fReactant1 = it->GetReactant1();
        fReactantB1 = it->GetReactant2();
      }else
      {
        fReactant1 = it->GetReactant2();
        fReactantB1 = it->GetReactant1();
      }
      for(const auto& itt : *(it->GetProducts()))
      {
        if(itt != H3OpB
            && itt != OHmB)
        {
          fReactant2 = itt;
        }else
        {
          fReactantB2 = itt;
        }
      }
      if(fVerbose > 1) {
        G4cout << "Equilibrium processes(ID) " << fRectionType << " : " << fReactant1->GetName()
               << " <=> " << fReactant2->GetName()
               << " Time to Equilibrium : " << fEquilibriumDuration / CLHEP::us
               << " Initial status : " << fAddEquilibrium << G4endl;
      }
      break ;
    }
  }
}

void G4ChemEquilibrium::PrintInfo() const
{
  G4cout<<"Equilibrium reactions : "<<fReactant1->GetName()
         <<" + "<<fReactantB1->GetName()
         <<" <=> "<<fReactant2->GetName()
         <<" + "<<fReactantB2->GetName()
         <<"  Status : "<<fAddEquilibrium
         <<" from "<<G4BestUnit(fEquilibriumTime,"Time")<<" to "
         <<G4BestUnit(fEquilibriumTime + fEquilibriumDuration,"Time")<<G4endl;
}


void G4ChemEquilibrium::SetEquilibrium(Reaction pReaction)
{
  if(pReaction->GetReactionType() != fRectionType)
  {
    std::vector<MolType> molVector;
    molVector.push_back(pReaction->GetReactant1());
    molVector.push_back(pReaction->GetReactant2());
    const G4int nbProducts = pReaction->GetNbProducts();
    if (nbProducts) {
      for (G4int j = 0; j < nbProducts; ++j) {
        auto product = pReaction->GetProduct(j);
        molVector.push_back(product);
      }
    }
    for(const auto& it : molVector)
    {
      if(it == fReactant1 || it == fReactant2 )
      {
        fAddEquilibrium = true;
        fEquilibriumTime = fGlobalTime;
        if(fVerbose >1) {
          G4cout << "Reaction type : " << pReaction->GetReactionType() << " : "
                 << pReaction->GetReactant1()->GetName() << " + "
                 << pReaction->GetReactant2()->GetName() << G4endl;
          G4cout << "SetEquilibrium : on " << fRectionType << "  fEquilibriumTime : "
                 << G4BestUnit(fEquilibriumTime, "Time")<<G4endl;
        }
        break;
      }
    }
  }
}