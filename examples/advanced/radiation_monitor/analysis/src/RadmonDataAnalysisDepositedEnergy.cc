//
// File name:     RadmonDataAnalysisDepositedEnergy.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisDepositedEnergy.cc,v 1.1 2005-11-25 01:56:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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

