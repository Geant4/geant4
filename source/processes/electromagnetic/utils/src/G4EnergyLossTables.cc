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
//
// -------------------------------------------------------------------
// first version created by P.Urban , 06/04/1998
// modifications + "precise" functions added by L.Urban , 27/05/98
// modifications , TOF functions , 26/10/98, L.Urban
// cache mechanism in order to gain time, 11/02/99, L.Urban
// bug fixed , 12/04/99 , L.Urban
// 10.11.99: moved from RWT hash dictionary to STL map, G.Barrand, M.Maire
// 27.09.01 L.Urban , bug fixed (negative energy deposit)
// 26.10.01 all static functions moved from .icc files (mma)
// 15.01.03 Add interfaces required for "cut per region" (V.Ivanchenko)
// 12.03.03 Add warnings to obsolete interfaces (V.Ivanchenko)
// 10.04.03 Add call to G4LossTableManager is particle is not registered (V.Ivanchenko)
//
// -------------------------------------------------------------------

#include "G4EnergyLossTables.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4RegionStore.hh"
#include "G4LossTableManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EnergyLossTablesHelper *G4EnergyLossTables::t = nullptr;
G4EnergyLossTablesHelper *G4EnergyLossTables::null_loss = nullptr;
G4ParticleDefinition* G4EnergyLossTables::lastParticle = nullptr;
G4double G4EnergyLossTables::QQPositron = 1.0; // e_squared
G4double G4EnergyLossTables::Chargesquare ;
G4int    G4EnergyLossTables::oldIndex = -1 ;
G4double G4EnergyLossTables::rmin = 0. ;
G4double G4EnergyLossTables::rmax = 0. ;
G4double G4EnergyLossTables::Thigh = 0. ;
G4int    G4EnergyLossTables::let_counter = 0;
G4int    G4EnergyLossTables::let_max_num_warnings = 100;
G4bool   G4EnergyLossTables::first_loss = true;

G4EnergyLossTables::helper_map *G4EnergyLossTables::dict = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EnergyLossTablesHelper::G4EnergyLossTablesHelper(
  const G4PhysicsTable* aDEDXTable,
  const G4PhysicsTable* aRangeTable,
  const G4PhysicsTable* anInverseRangeTable,
  const G4PhysicsTable* aLabTimeTable,
  const G4PhysicsTable* aProperTimeTable,
  G4double aLowestKineticEnergy,
  G4double aHighestKineticEnergy,
  G4double aMassRatio,
  G4int aNumberOfBins)
  :
  theDEDXTable(aDEDXTable), theRangeTable(aRangeTable),
  theInverseRangeTable(anInverseRangeTable),
  theLabTimeTable(aLabTimeTable),
  theProperTimeTable(aProperTimeTable),
  theLowestKineticEnergy(aLowestKineticEnergy),
  theHighestKineticEnergy(aHighestKineticEnergy),
  theMassRatio(aMassRatio),
  theNumberOfBins(aNumberOfBins)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EnergyLossTablesHelper::G4EnergyLossTablesHelper()
{ 
  theLowestKineticEnergy = 0.0;
  theHighestKineticEnergy= 0.0;
  theMassRatio = 0.0;
  theNumberOfBins = 0;
  theDEDXTable = theRangeTable = theInverseRangeTable = theLabTimeTable 
    = theProperTimeTable = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EnergyLossTables::Register(
  const G4ParticleDefinition* p,
  const G4PhysicsTable* tDEDX,
  const G4PhysicsTable* tRange,
  const G4PhysicsTable* tInverseRange,
  const G4PhysicsTable* tLabTime,
  const G4PhysicsTable* tProperTime,
  G4double lowestKineticEnergy,
  G4double highestKineticEnergy,
  G4double massRatio,
  G4int NumberOfBins)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  if (!null_loss) null_loss = new G4EnergyLossTablesHelper;
  if (!t) t = new G4EnergyLossTablesHelper;

  (*dict)[p]= G4EnergyLossTablesHelper(tDEDX, tRange,tInverseRange,
                    tLabTime,tProperTime,lowestKineticEnergy,
		    highestKineticEnergy, massRatio,NumberOfBins);

  *t = GetTables(p) ;    // important for cache !!!!!
  lastParticle = (G4ParticleDefinition*) p ;
  Chargesquare = (p->GetPDGCharge())*(p->GetPDGCharge())/
                  QQPositron ;
  if (first_loss ) {
    *null_loss = G4EnergyLossTablesHelper(
                 nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, 0.0, 0.0, 0);
    first_loss = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4PhysicsTable* G4EnergyLossTables::GetDEDXTable(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  helper_map::iterator it;
  if((it=dict->find(p))==dict->end()) return nullptr;
  return (*it).second.theDEDXTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4PhysicsTable* G4EnergyLossTables::GetRangeTable(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  helper_map::iterator it;
  if((it=dict->find(p))==dict->end()) return nullptr;
  return (*it).second.theRangeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4PhysicsTable* G4EnergyLossTables::GetInverseRangeTable(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  helper_map::iterator it;
  if((it=dict->find(p))==dict->end()) return nullptr;
  return (*it).second.theInverseRangeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4PhysicsTable* G4EnergyLossTables::GetLabTimeTable(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  helper_map::iterator it;
  if((it=dict->find(p))==dict->end()) return nullptr;
  return (*it).second.theLabTimeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4PhysicsTable* G4EnergyLossTables::GetProperTimeTable(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  helper_map::iterator it;
  if((it=dict->find(p))==dict->end()) return nullptr;
  return (*it).second.theProperTimeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EnergyLossTablesHelper G4EnergyLossTables::GetTables(
  const G4ParticleDefinition* p)
{
  if (!dict) dict = new G4EnergyLossTables::helper_map;
  if (!null_loss) null_loss = new G4EnergyLossTablesHelper;

  helper_map::iterator it;
  if ((it=dict->find(p))==dict->end()) {
    return *null_loss;
  }
  return (*it).second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                   QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if (!dEdxTable) {
    ParticleHaveNoLoss(aParticle,"dEdx");
    return 0.0;
  }

  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     dEdx =(*dEdxTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut)
           *std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut);

  } else {

    dEdx = (*dEdxTable)(materialIndex)->GetValue(
	       scaledKineticEnergy,isOut);

  }

  return dEdx*Chargesquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* labtimeTable= t->theLabTimeTable;
  if (!labtimeTable) {
    ParticleHaveNoLoss(aParticle,"LabTime");
    return 0.0;
  }

  const G4double parlowen=0.4 , ppar=0.5-parlowen ;
  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double time;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     time = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
            (*labtimeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     time = (*labtimeTable)(materialIndex)->GetValue(
              t->theHighestKineticEnergy,isOut);

  } else {

    time = (*labtimeTable)(materialIndex)->GetValue(
               scaledKineticEnergy,isOut);

  }

  return time/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetDeltaLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* labtimeTable= t->theLabTimeTable;
  if (!labtimeTable) {
    ParticleHaveNoLoss(aParticle,"LabTime");
    return 0.0;
  }

  const G4double parlowen=0.4 , ppar=0.5-parlowen ;
  const G4double dToverT = 0.05 , facT = 1. -dToverT ;
  G4double timestart,timeend,deltatime,dTT;
  G4bool isOut;

  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergyStart*t->theMassRatio;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     timestart = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
                (*labtimeTable)(materialIndex)->GetValue(
                t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     timestart = (*labtimeTable)(materialIndex)->GetValue(
                t->theHighestKineticEnergy,isOut);

  } else {

    timestart = (*labtimeTable)(materialIndex)->GetValue(
                scaledKineticEnergy,isOut);

  }

  dTT = (KineticEnergyStart - KineticEnergyEnd)/KineticEnergyStart ;

  if( dTT < dToverT )
    scaledKineticEnergy = facT*KineticEnergyStart*t->theMassRatio;
  else
    scaledKineticEnergy = KineticEnergyEnd*t->theMassRatio;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     timeend = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
                (*labtimeTable)(materialIndex)->GetValue(
                t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     timeend = (*labtimeTable)(materialIndex)->GetValue(
                t->theHighestKineticEnergy,isOut);

  } else {

    timeend = (*labtimeTable)(materialIndex)->GetValue(
                scaledKineticEnergy,isOut);

  }

  deltatime = timestart - timeend ;

  if( dTT < dToverT )
    deltatime *= dTT/dToverT;

  return deltatime/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* propertimeTable= t->theProperTimeTable;
  if (!propertimeTable) {
    ParticleHaveNoLoss(aParticle,"ProperTime");
    return 0.0;
  }

  const G4double parlowen=0.4 , ppar=0.5-parlowen ;
  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double time;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     time = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
            (*propertimeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     time = (*propertimeTable)(materialIndex)->GetValue(
              t->theHighestKineticEnergy,isOut);

  } else {

    time = (*propertimeTable)(materialIndex)->GetValue(
               scaledKineticEnergy,isOut);

  }

  return time/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetDeltaProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* propertimeTable= t->theProperTimeTable;
  if (!propertimeTable) {
    ParticleHaveNoLoss(aParticle,"ProperTime");
    return 0.0;
  }

  const G4double parlowen=0.4 , ppar=0.5-parlowen ;
  const G4double dToverT = 0.05 , facT = 1. -dToverT ;
  G4double timestart,timeend,deltatime,dTT;
  G4bool isOut;

  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergyStart*t->theMassRatio;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     timestart = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
                (*propertimeTable)(materialIndex)->GetValue(
                t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     timestart = (*propertimeTable)(materialIndex)->GetValue(
                t->theHighestKineticEnergy,isOut);

  } else {

    timestart = (*propertimeTable)(materialIndex)->GetValue(
                scaledKineticEnergy,isOut);

  }

  dTT = (KineticEnergyStart - KineticEnergyEnd)/KineticEnergyStart ;

  if( dTT < dToverT )
    scaledKineticEnergy = facT*KineticEnergyStart*t->theMassRatio;
  else
    scaledKineticEnergy = KineticEnergyEnd*t->theMassRatio;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     timeend = std::exp(ppar*std::log(scaledKineticEnergy/t->theLowestKineticEnergy))*
                (*propertimeTable)(materialIndex)->GetValue(
                t->theLowestKineticEnergy,isOut);


  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     timeend = (*propertimeTable)(materialIndex)->GetValue(
                t->theHighestKineticEnergy,isOut);

  } else {

    timeend = (*propertimeTable)(materialIndex)->GetValue(
                scaledKineticEnergy,isOut);

  }

  deltatime = timestart - timeend ;

  if( dTT < dToverT )
    deltatime *= dTT/dToverT ;

  return deltatime/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetRange(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* rangeTable= t->theRangeTable;
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if (!rangeTable) {
    ParticleHaveNoLoss(aParticle,"Range");
    return 0.0;
  }

  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

    Range = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)*
            (*rangeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut)+
            (scaledKineticEnergy-t->theHighestKineticEnergy)/
            (*dEdxTable)(materialIndex)->GetValue(
              t->theHighestKineticEnergy,isOut);

  } else {

    Range = (*rangeTable)(materialIndex)->GetValue(
	       scaledKineticEnergy,isOut);

  }

  return Range/(Chargesquare*t->theMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetPreciseEnergyFromRange(
                                     const G4ParticleDefinition *aParticle,
                                           G4double range,
                                     const G4Material *aMaterial)
// it returns the value of the kinetic energy for a given range
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  const G4PhysicsTable*  inverseRangeTable= t->theInverseRangeTable;
  if (!inverseRangeTable) {
    ParticleHaveNoLoss(aParticle,"InverseRange");
    return 0.0;
  }

  G4double scaledrange,scaledKineticEnergy ;
  G4bool isOut ;

  G4int materialIndex = (G4int)aMaterial->GetIndex() ;

  if(materialIndex != oldIndex)
  {
    oldIndex = materialIndex ;
    rmin = (*inverseRangeTable)(materialIndex)->
                              GetLowEdgeEnergy(0) ;
    rmax = (*inverseRangeTable)(materialIndex)->
                   GetLowEdgeEnergy(t->theNumberOfBins-2) ;
    Thigh = (*inverseRangeTable)(materialIndex)->
                              GetValue(rmax,isOut) ;
  }

  scaledrange = range*Chargesquare*t->theMassRatio ;

  if(scaledrange < rmin)
  {
    scaledKineticEnergy = t->theLowestKineticEnergy*
                   scaledrange*scaledrange/(rmin*rmin) ;
  }
  else
  {
    if(scaledrange < rmax)
    {
      scaledKineticEnergy = (*inverseRangeTable)(materialIndex)->
                              GetValue( scaledrange,isOut) ;
    }
    else
    {
      scaledKineticEnergy = Thigh +
                      (scaledrange-rmax)*
                      (*dEdxTable)(materialIndex)->
                                 GetValue(Thigh,isOut) ;
    }
  }

  return scaledKineticEnergy/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4double G4EnergyLossTables::GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if (!dEdxTable) {
    ParticleHaveNoLoss(aParticle,"dEdx");
    return 0.0;
  }

  G4int materialIndex = (G4int)aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     dEdx = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)
            *(*dEdxTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut);

  } else {

      dEdx = (*dEdxTable)(materialIndex)->GetValue(
                          scaledKineticEnergy,isOut) ;

  }

  return dEdx*Chargesquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4double G4EnergyLossTables::GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  CPRWarning();
  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* rangeTable= t->theRangeTable;
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if (!rangeTable) {
    ParticleHaveNoLoss(aParticle,"Range");
    return 0.0;
  }
  G4int materialIndex = (G4int)aMaterial->GetIndex();

  G4double Thighr = t->theHighestKineticEnergy*t->theLowestKineticEnergy/
                   (*rangeTable)(materialIndex)->
                   GetLowEdgeEnergy(1) ;

  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

    Range = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)*
            (*rangeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>Thighr) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      Thighr,isOut)+
            (scaledKineticEnergy-Thighr)/
            (*dEdxTable)(materialIndex)->GetValue(
              Thighr,isOut);

  } else {

     Range = (*rangeTable)(materialIndex)->GetValue(
                       scaledKineticEnergy,isOut) ;

  }

  return Range/(Chargesquare*t->theMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple,
    G4bool check)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                   QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  
  if (!dEdxTable ) {
    if (check) return G4LossTableManager::Instance()->GetDEDX(aParticle,KineticEnergy,couple);
    else       ParticleHaveNoLoss(aParticle, "dEdx");
    return 0.0;
  }

  G4int materialIndex = couple->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     dEdx =(*dEdxTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut)
           *std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut);

  } else {

    dEdx = (*dEdxTable)(materialIndex)->GetValue(
	       scaledKineticEnergy,isOut);

  }

  return dEdx*Chargesquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetRange(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple,
    G4bool check)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  if(aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle ;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* rangeTable= t->theRangeTable;
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if (!rangeTable) {
    if(check) return G4LossTableManager::Instance()->GetRange(aParticle,KineticEnergy,couple);
    else      return DBL_MAX;      
      //ParticleHaveNoLoss(aParticle,"Range");
  }

  G4int materialIndex = couple->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

    Range = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)*
            (*rangeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut)+
            (scaledKineticEnergy-t->theHighestKineticEnergy)/
            (*dEdxTable)(materialIndex)->GetValue(
              t->theHighestKineticEnergy,isOut);

  } else {

    Range = (*rangeTable)(materialIndex)->GetValue(
	       scaledKineticEnergy,isOut);

  }

  return Range/(Chargesquare*t->theMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetPreciseEnergyFromRange(
                                     const G4ParticleDefinition *aParticle,
                                           G4double range,
                                     const G4MaterialCutsCouple *couple,
				           G4bool check)
// it returns the value of the kinetic energy for a given range
{
  if (!t) t = new G4EnergyLossTablesHelper;

  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  const G4PhysicsTable*  inverseRangeTable= t->theInverseRangeTable;
  
  if (!inverseRangeTable) {
    if(check) return G4LossTableManager::Instance()->GetEnergy(aParticle,range,couple);
    else      return DBL_MAX;      
    //    else      ParticleHaveNoLoss(aParticle,"InverseRange");
  }

  G4double scaledrange,scaledKineticEnergy ;
  G4bool isOut ;

  G4int materialIndex = couple->GetIndex() ;

  if(materialIndex != oldIndex)
  {
    oldIndex = materialIndex ;
    rmin = (*inverseRangeTable)(materialIndex)->
                              GetLowEdgeEnergy(0) ;
    rmax = (*inverseRangeTable)(materialIndex)->
                   GetLowEdgeEnergy(t->theNumberOfBins-2) ;
    Thigh = (*inverseRangeTable)(materialIndex)->
                              GetValue(rmax,isOut) ;
  }

  scaledrange = range*Chargesquare*t->theMassRatio ;

  if(scaledrange < rmin)
  {
    scaledKineticEnergy = t->theLowestKineticEnergy*
                   scaledrange*scaledrange/(rmin*rmin) ;
  }
  else
  {
    if(scaledrange < rmax)
    {
      scaledKineticEnergy = (*inverseRangeTable)(materialIndex)->
                              GetValue( scaledrange,isOut) ;
    }
    else
    {
      scaledKineticEnergy = Thigh +
                      (scaledrange-rmax)*
                      (*dEdxTable)(materialIndex)->
                                 GetValue(Thigh,isOut) ;
    }
  }

  return scaledKineticEnergy/t->theMassRatio ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if ( !dEdxTable )
    return G4LossTableManager::Instance()->GetDEDX(aParticle,KineticEnergy,couple);

  G4int materialIndex = couple->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

     dEdx = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)
            *(*dEdxTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t->theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t->theHighestKineticEnergy,isOut);

  } else {

      dEdx = (*dEdxTable)(materialIndex)->GetValue(
                          scaledKineticEnergy,isOut) ;

  }

  return dEdx*Chargesquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4EnergyLossTables::GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple)
{
  if (!t) t = new G4EnergyLossTablesHelper;

  if( aParticle != (const G4ParticleDefinition*) lastParticle)
  {
    *t= GetTables(aParticle);
    lastParticle = (G4ParticleDefinition*) aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
    oldIndex = -1 ;
  }
  const G4PhysicsTable* rangeTable= t->theRangeTable;
  const G4PhysicsTable*  dEdxTable= t->theDEDXTable;
  if ( !dEdxTable || !rangeTable)
    return G4LossTableManager::Instance()->GetDEDX(aParticle,KineticEnergy,couple);

  G4int materialIndex = couple->GetIndex();

  G4double Thighr = t->theHighestKineticEnergy*t->theLowestKineticEnergy/
                   (*rangeTable)(materialIndex)->
                   GetLowEdgeEnergy(1) ;

  G4double scaledKineticEnergy = KineticEnergy*t->theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t->theLowestKineticEnergy) {

    Range = std::sqrt(scaledKineticEnergy/t->theLowestKineticEnergy)*
            (*rangeTable)(materialIndex)->GetValue(
              t->theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>Thighr) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      Thighr,isOut)+
            (scaledKineticEnergy-Thighr)/
            (*dEdxTable)(materialIndex)->GetValue(
              Thighr,isOut);

  } else {

     Range = (*rangeTable)(materialIndex)->GetValue(
                       scaledKineticEnergy,isOut) ;

  }

  return Range/(Chargesquare*t->theMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EnergyLossTables::CPRWarning()
{
  if (let_counter <  let_max_num_warnings) {

    G4cout << G4endl;
    G4cout << "##### G4EnergyLossTable WARNING: The obsolete interface is used!" << G4endl;
    G4cout << "##### RESULTS ARE NOT GARANTEED!" << G4endl;
    G4cout << "##### Please, substitute G4Material by G4MaterialCutsCouple" << G4endl;
    G4cout << "##### Obsolete interface will be removed soon" << G4endl;
    G4cout << G4endl;
    let_counter++;

  } else if (let_counter == let_max_num_warnings) {

    G4cout << "##### G4EnergyLossTable WARNING closed" << G4endl;
    let_counter++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EnergyLossTables::ParticleHaveNoLoss(const G4ParticleDefinition*, 
				       const G4String& /*q*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
