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
// $Id: G4PSDoseDeposit_RegNav.cc,v 1.1.2.1 2009/03/03 13:43:58 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// G4PSDoseDeposit_RegNav
#include "G4PSDoseDeposit_RegNav.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"
#include "G4PhantomParameterisation.hh"
#include "G4RegularNavigationHelper.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4EmCalculator.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVParameterised.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only energy deposit.
// 
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSDoseDeposit_RegNav::G4PSDoseDeposit_RegNav(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1)
{
   theElossExt = new G4EnergyLossForExtrapolator(0);
   thePhantomParam = 0;
}

G4PSDoseDeposit_RegNav::~G4PSDoseDeposit_RegNav()
{;}

G4bool G4PSDoseDeposit_RegNav::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( edep == 0. ) return FALSE;
  G4double volume  = aStep->GetPreStepPoint()->GetPhysicalVolume()
    ->GetLogicalVolume()->GetSolid()->GetCubicVolume();
  G4double density = aStep->GetTrack()->GetMaterial()->GetDensity();

  G4bool verbose = 0;

  if( aStep != 0 ) { // it is 0 when called by GmScoringMgr after last event
    if( G4RegularNavigationHelper::theStepLengths.size() <= 1 || aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0)  { // we are only counting dose deposit
      G4double dose    = edep / ( density * volume );
      dose *= aStep->GetPreStepPoint()->GetWeight(); 
      G4int  index = GetIndex(aStep);
      EvtMap->add(index,dose); 
    } else {
      std::vector< std::pair<G4int,G4double> > rnsl = G4RegularNavigationHelper::theStepLengths; 
      if( !thePhantomParam ) thePhantomParam = GetPhantomParam(true);
      const G4ParticleDefinition* part = aStep->GetTrack()->GetDefinition();
      G4double kinEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();

      G4double stepLength = aStep->GetStepLength();
      G4double slSum = 0.;
      size_t ii;
      for( ii = 0; ii < rnsl.size(); ii++ ){
	G4double sl = rnsl[ii].second;
	slSum += sl;
      }

      //--- 1st: distribute true step length in each voxel: geom SL in each voxel is multiplied by a constant so that the sum gives the total true step length
      std::vector<G4double> stepLengths;
      G4double slRatio = stepLength/slSum;
      for( ii = 0; ii < rnsl.size(); ii++ ){
	G4double sl = rnsl[ii].second;
	stepLengths.push_back( sl * slRatio );
	if(verbose) G4cout << ii << " RN: 1st step length geom " << sl << " true " << stepLengths[ii] << G4endl;
      }

      //--- energy at each interaction
      slSum = 0.;
      G4EmCalculator emcalc;
      G4double totalELost = 0.;
      std::vector<G4double> kinEnergies;
      kinEnergies.push_back( kinEnergy );
      for( ii = 0; ii < rnsl.size(); ii++ ){
	const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	G4double elost = stepLengths[ii] * emcalc.GetDEDX(kinEnergy, part, mate);
	kinEnergy -= elost;
	kinEnergies.push_back( kinEnergy );
	totalELost += elost;
	if(verbose) G4cout << ii << " RN: 1st energy lost " << elost << " energy at interaction " << kinEnergy << " dEdx " << emcalc.GetDEDX(kinEnergy, part, mate) << G4endl;
      }
      //correct energies so that they reproduce the real step energy lost
      G4double enerRatio = -(aStep->GetDeltaEnergy()/totalELost);
      if(verbose) G4cout << ii << " RN: 1st energy ratio " << enerRatio << G4endl;
      for( ii = kinEnergies.size()-1; ii >= 1; ii-- ){
	G4double elost = ( kinEnergies[ii-1]-kinEnergies[ii] )*enerRatio;
	//if(verbose) G4cout << " kin old " << kinEnergies[ii-1] << " " << kinEnergies[ii] << " " << elost/enerRatio << G4endl;
	kinEnergies[ii] = kinEnergies[ii-1] - elost;
	if(verbose) G4cout << ii << " RN: 1st corrected energy lost " << elost << " energy at interaction " << kinEnergies[ii] << " orig elost " << elost/enerRatio << G4endl;
      }

      //----- Get step lengths corrected by changing geom2true correction
      //-- Get ratios for each energy 
      slSum = 0.;
      for( ii = 0; ii < rnsl.size(); ii++ ){
	const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );

	//correct step length
	stepLengths[ii] = theElossExt->TrueStepLength( kinEnergies[ii] , rnsl[ii].second , mate, part );

	if(verbose) G4cout << ii << " RN: 2nd step length " << stepLengths[ii] << " g2t " << rnsl[ii].second / stepLengths[ii]  << G4endl;
	slSum += stepLengths[ii];
      }
      //Correct so that they sum the total step lengths
      G4double slratio = aStep->GetStepLength()/slSum;
      if(verbose) G4cout << ii << " RN: 2nd step ratio " << slRatio << G4endl;
      for( ii = 0; ii < rnsl.size(); ii++ ){
	stepLengths[ii] *= slratio;
	if(verbose) G4cout << ii << " RN: 2nd corrected step length " << stepLengths[ii] << " g2t " << G4endl;
      }

      //---- Recalculate energy lost with this new step lengths
      kinEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
      totalELost = 0.;
      for( ii = 0; ii < rnsl.size(); ii++ ){
	const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	G4double elost = stepLengths[ii] * emcalc.GetDEDX(kinEnergy, part, mate);
	kinEnergy -= elost;
	kinEnergies[ii+1] = kinEnergy;
	totalELost += elost;
	if(verbose) G4cout << ii << " RN: 2nd energy lost " << elost << " energy at interaction " << kinEnergy << G4endl;
      }
      //correct energies so that they reproduce the real step energy lost
      enerRatio = -(aStep->GetDeltaEnergy()/totalELost);
      if(verbose) G4cout << ii << " RN: 2nd energy ratio " << enerRatio << G4endl;
      G4double wei = aStep->GetPreStepPoint()->GetWeight(); 

      for( ii = kinEnergies.size()-1; ii >= 1; ii-- ){
	G4double elost = ( kinEnergies[ii-1]-kinEnergies[ii] )*enerRatio;
	kinEnergies[ii] = kinEnergies[ii-1] - elost;
	G4double dose  = elost / ( density * volume );
	G4double valwei = dose*wei;
	G4int index = rnsl[ii-1].first;

	EvtMap->add(index,valwei);
	if(verbose) G4cout << ii << " RN: 2nd corrected energy lost " << elost << " energy at interaction " << kinEnergies[ii] << G4endl;
	if(verbose) G4cout << ii << " RN: deposited energy " << elost*wei << " dose " << valwei << " in " << index << G4endl;
	
      }
      
    }
  }
 
  return TRUE;
}

void G4PSDoseDeposit_RegNav::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSDoseDeposit_RegNav::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSDoseDeposit_RegNav::clear()
{
  EvtMap->clear();
}

void G4PSDoseDeposit_RegNav::DrawAll()
{;}

void G4PSDoseDeposit_RegNav::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  dose deposit: " << G4BestUnit(*(itr->second),"Dose")
	   << G4endl;
  }
}



//-----------------------------------------------------------------------
G4PhantomParameterisation* G4PSDoseDeposit_RegNav::GetPhantomParam(G4bool mustExist)
{
  G4PhantomParameterisation* paramreg = 0;

  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*>::iterator cite;
  for( cite = pvs->begin(); cite != pvs->end(); cite++ ) {
    //    G4cout << " PV " << (*cite)->GetName() << " " << (*cite)->GetTranslation() << G4endl;
    if( IsPhantomVolume( *cite ) ) {
      const G4PVParameterised* pvparam = static_cast<const G4PVParameterised*>(*cite);
      G4VPVParameterisation* param = pvparam->GetParameterisation();
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //      G4cout << "G4PhantomParameterisation volume found  " << (*cite)->GetName() << G4endl;
      paramreg = static_cast<G4PhantomParameterisation*>(param);
    }
  }
  
  if( !paramreg && mustExist ) G4Exception("GmRegularParamUtils::GetPhantomParam:  No G4PhantomParameterisation found ");
  
  return paramreg;
  
}


//-----------------------------------------------------------------------
G4bool G4PSDoseDeposit_RegNav::IsPhantomVolume( G4VPhysicalVolume* pv )
{
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;
  pv->GetReplicationData(axis,nReplicas,width,offset,consuming);
  EVolume type = (consuming) ? kReplica : kParameterised;
  if( type == kParameterised && pv->GetRegularStructureId() == 1 ) {  
    return TRUE;
  } else {
    return FALSE;
  }

} 

