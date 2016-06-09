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
// File name:     RadmonDataAnalysisDepositedEnergy.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisDepositedEnergy.cc,v 1.3 2006/06/29 16:07:43 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//

#include "RadmonDataAnalysisDepositedEnergy.hh"
#include "RadmonHitsManager.hh"
#include "G4Step.hh"
#include "AIDA/ITuple.h"

                                                RadmonDataAnalysisDepositedEnergy :: RadmonDataAnalysisDepositedEnergy()
:
 RadmonVDataAnalysisWithLabel("DepositedEnergy"),
 indexHitEnergyDeposit(RadmonHitsManager::Instance()->ReserveIdByLabel("EnergyDeposit"))
{
}


 
 
 
RadmonVDataAnalysisWithLabel *                  RadmonDataAnalysisDepositedEnergy :: New(void) const
{
 return new RadmonDataAnalysisDepositedEnergy;
}





G4String                                        RadmonDataAnalysisDepositedEnergy :: ObtainColumnsDeclaration(const G4String & prefix)
{
 G4String result("double ");
 result+=prefix;
 result+="_EnergyDeposit";
 
 return result;
}





void                                            RadmonDataAnalysisDepositedEnergy :: InitializeFromTuple(const G4String & prefix, const AIDA::ITuple * tuple)
{
 indexTupleEnergyDeposit=tuple->findColumn(prefix+"_EnergyDeposit");
}





void                                            RadmonDataAnalysisDepositedEnergy :: StoreIntoTuple(RadmonHitsCollection * hitsCollection, AIDA::ITuple * tuple)
{
 G4int n(hitsCollection->entries());
 G4double totalEnergy(0);
 
 while (n>0)
 {
  n--;
  totalEnergy+=(*hitsCollection)[n]->RetrieveById(indexHitEnergyDeposit);
 }
 
 tuple->fill(indexTupleEnergyDeposit, totalEnergy/MeV);
}





void                                            RadmonDataAnalysisDepositedEnergy :: StoreIntoHit(G4Step * step, RadmonHit * hit)
{
 hit->StoreById(indexHitEnergyDeposit, step->GetTotalEnergyDeposit());
}

