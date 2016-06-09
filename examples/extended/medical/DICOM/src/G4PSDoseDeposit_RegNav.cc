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
//#define VERBOSE_DOSEDEP
//
// $Id: G4PSDoseDeposit_RegNav.cc,v 1.4 2009/12/16 17:54:21 gunter Exp $
// GEANT4 tag $Name: geant4-09-03 $
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
   theNIterations = 2;
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

#ifdef VERBOSE_DOSEDEP
  G4bool verbose = 1;
#endif

 if( aStep == 0 ) return FALSE; // it is 0 when called by GmScoringMgr after last event

#ifdef VERBOSE_DOSEDEP
 if( verbose ) G4cout << "GmG4PSDoseDeposit::FillScorer totalEdepo " << aStep->GetTotalEnergyDeposit() 
		      << " Nsteps " << G4RegularNavigationHelper::theStepLengths.size() << G4endl;
#endif
  //----- Do not distribute dose in voxels 
 G4double dose  = edep / ( density * volume );
 G4double wei = aStep->GetPreStepPoint()->GetWeight(); 

  if( G4RegularNavigationHelper::theStepLengths.size() <= 1 || 
      aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0)  { // we are only counting dose deposit
    dose *= wei;
    G4int  index = GetIndex(aStep);
    EvtMap->add(index,dose); 
#ifdef VERBOSE_DOSEDEP
    if( verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer  RN: energy lost " << dose << " index " << index << G4endl;
#endif
  } else {
    //----- Distribute dose in voxels 
    std::vector< std::pair<G4int,G4double> > rnsl = G4RegularNavigationHelper::theStepLengths; 
    //      G4double geomSL = rnsl[0].second;
    //      G4double SL = aStep->GetStepLength(); 
    if( !thePhantomParam ) thePhantomParam = GetPhantomParam(true);
    //      const G4Material* mate = thePhantomParam->GetMaterial( rnsl[0].first );
    const G4ParticleDefinition* part = aStep->GetTrack()->GetDefinition();
    G4double kinEnergyPreOrig = aStep->GetPreStepPoint()->GetKineticEnergy();
    G4double kinEnergyPre = kinEnergyPreOrig;
    
    G4double stepLength = aStep->GetStepLength();
    G4double slSum = 0.;
    unsigned int ii;
    for( ii = 0; ii < rnsl.size(); ii++ ){
      G4double sl = rnsl[ii].second;
      slSum += sl;
#ifdef VERBOSE_DOSEDEP
      if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: it\
er1 step length geom " << sl << G4endl;
#endif
    }
    
#ifdef VERBOSE_DOSEDEP
    if( verbose )
      G4cout << "GmG4PSDoseDeposit RN:  step length geom TOTAL " << slSum 
	     << " true TOTAL " << stepLength 
	     << " ratio " << stepLength/slSum 
	     << " Energy " << aStep->GetPreStepPoint()->GetKineticEnergy() 
	     << " Material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() 
	     << " Number of geom steps " << rnsl.size() << G4endl;
#endif
    //----- No iterations to correct elost and msc, distribute dose according to geometrical step length in each voxel
    if( theNIterations == 0 ) { 
      for( unsigned int ii = 0; ii < rnsl.size(); ii++ ){
	G4int index = G4RegularNavigationHelper::theStepLengths[ii].first;
	G4double sl = G4RegularNavigationHelper::theStepLengths[ii].second;
	G4double doseStep = dose * sl/slSum; //divide dose along steps, proportional to step lengthr
	G4double dosewei = doseStep*wei;
#ifdef VERBOSE_DOSEDEP
	if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii 
			    << " dose " << dosewei 
			    << " in " << index << G4endl;
#endif
	
	EvtMap->add(index, dosewei );
	
      }
    } else { //  1 or more iterations demanded
      
#ifdef VERBOSE_DOSEDEP
      // print corrected energy at iteration 0 
      if(verbose)  {
	G4double slSum = 0.;
	for( ii = 0; ii < rnsl.size(); ii++ ){
	  G4double sl = rnsl[ii].second;
	  slSum += sl;
	}
	
	for( ii = 0; ii < rnsl.size(); ii++ ){
	  G4cout  << "GmG4PSDoseDeposit::FillScorer "<< ii
		  << " RN: iter0 corrected energy lost " << aStep->GetTotalEnergyDeposit()*rnsl[ii].second/slSum  
		  << G4endl;
	}
      }
#endif
      G4double slRatio = stepLength/slSum;
#ifdef VERBOSE_DOSEDEP
      if(verbose) G4cout << "GmG4PSDoseDeposit::FillScorer  RN: iter" << iiter << " step ratio " << slRatio << G4endl;
#endif
      
      //--- energy at each interaction
      G4EmCalculator emcalc;
      G4double totalELost = 0.;
      std::vector<G4double> kinELost;
      std::vector<G4double> stepLengths;
      for( int iiter = 1; iiter <= theNIterations; iiter++ ) {
	//--- iter1: distribute true step length in each voxel: geom SL in each voxel is multiplied by a constant so that the sum gives the total true step length
	if( iiter == 1 ) {
	  for( ii = 0; ii < rnsl.size(); ii++ ){
	    G4double sl = rnsl[ii].second;
	    stepLengths.push_back( sl * slRatio );
#ifdef VERBOSE_DOSEDEP
	    if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: iter" << iiter << " corrected step length " << sl*slRatio << G4endl;
#endif
	  }
	  
	  for( ii = 0; ii < rnsl.size(); ii++ ){
	    const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	    G4double dEdx = 0.;
	    if( kinEnergyPre > 0. ) {  //t check this 
	      dEdx = emcalc.GetDEDX(kinEnergyPre, part, mate);
	    }
	    G4double elost = stepLengths[ii] * dEdx;
	    
#ifdef VERBOSE_DOSEDEP
	    if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: iter1 energy lost "  << elost 
				<< " energy at interaction " << kinEnergyPre 
				<< " = stepLength " << stepLengths[ii] 
				<< " * dEdx " << dEdx << G4endl;
#endif
	    kinEnergyPre -= elost;
	    kinELost.push_back( elost );
	    totalELost += elost;
	  }
	  
	} else{
	  //------ 2nd and other iterations
	  //----- Get step lengths corrected by changing geom2true correction
	  //-- Get ratios for each energy 
	  slSum = 0.;
	  kinEnergyPre = kinEnergyPreOrig;
	  for( ii = 0; ii < rnsl.size(); ii++ ){
	    const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	    stepLengths[ii] = theElossExt->TrueStepLength( kinEnergyPre, rnsl[ii].second , mate, part );
	    kinEnergyPre -= kinELost[ii];
	    
#ifdef VERBOSE_DOSEDEP
	    if(verbose) G4cout << "GmG4PSDoseDeposit::FillScorer" << ii 
			       << " RN: iter" << iiter << " step length geom " << stepLengths[ii] 
			       << " geom2true " << rnsl[ii].second / stepLengths[ii]  << G4endl;
#endif
	    
	    slSum += stepLengths[ii];
	  }
	  
	  //Correct step lengths so that they sum the total step length
	  G4double slratio = aStep->GetStepLength()/slSum;
#ifdef VERBOSE_DOSEDEP
	  if(verbose) G4cout << "GmG4PSDoseDeposit::FillScorer" << ii << " RN: iter" << iiter << " step ratio " << slRatio << G4endl;
#endif
	  for( ii = 0; ii < rnsl.size(); ii++ ){
	    stepLengths[ii] *= slratio;
#ifdef VERBOSE_DOSEDEP
	    if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: iter" << iiter << " corrected step length " << stepLengths[ii] << G4endl;
#endif
	    }
	  
	  //---- Recalculate energy lost with this new step lengths
	  G4double kinEnergyPre = aStep->GetPreStepPoint()->GetKineticEnergy();
	  totalELost = 0.;
	  for( ii = 0; ii < rnsl.size(); ii++ ){
	    const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	    G4double dEdx = 0.;
	    if( kinEnergyPre > 0. ) {
	      dEdx = emcalc.GetDEDX(kinEnergyPre, part, mate);
	    }
	    G4double elost = stepLengths[ii] * dEdx;
#ifdef VERBOSE_DOSEDEP
	    if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: iter" << iiter << " energy lost " << elost 
				<< " energy at interaction " << kinEnergyPre 
				<< " = stepLength " << stepLengths[ii] 
				<< " * dEdx " << dEdx << G4endl;
#endif
	    kinEnergyPre -= elost;
	    kinELost[ii] = elost;
	    totalELost += elost;
	  }
	  
	}
	
	//correct energies so that they reproduce the real step energy lost
	G4double enerRatio = (aStep->GetTotalEnergyDeposit()/totalELost);
	
#ifdef VERBOSE_DOSEDEP
	if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii << " RN: iter" << iiter << " energy ratio " << enerRatio << G4endl;
#endif
	
#ifdef VERBOSE_DOSEDEP
	G4double elostTot = 0.; 
#endif
	for( ii = 0; ii < kinELost.size(); ii++ ){
	  kinELost[ii] *= enerRatio;
#ifdef VERBOSE_DOSEDEP
	  elostTot += kinELost[ii];
	  if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer "<< ii << " RN: iter" << iiter << " corrected energy lost " << kinELost[ii] 
			      << " orig elost " << kinELost[ii]/enerRatio 
			      << " energy before interaction " << kinEnergyPreOrig-elostTot+kinELost[ii]
			      << " energy after interaction " << kinEnergyPreOrig-elostTot
			      << G4endl;
#endif
	}
      }
      
      //---- Compute the dose (for N iterations)
      G4double dosesum = 0.;
      G4double volume  = aStep->GetPreStepPoint()->GetPhysicalVolume()
	->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      G4double density = aStep->GetTrack()->GetMaterial()->GetDensity();
      
      for( ii = 0; ii < kinELost.size(); ii++ ){
	G4double dose  = kinELost[ii] / ( density * volume );
	G4double dosewei = dose*wei;
	G4int index = rnsl[ii].first;
#ifdef VERBOSE_DOSEDEP
	if(verbose) G4cout  << "GmG4PSDoseDeposit::FillScorer"<< ii 
			    << " dose " << dosewei
			    << " in " << index
			    << " RN: deposited energy " << kinELost[ii]*wei 
			    << G4endl;
#endif
	
	EvtMap->add(index, dosewei );
	dosesum += dosewei;
	
      }
      
#ifdef VERBOSE_DOSEDEP
      if(verbose) {
	//      G4cout << "DOSE-DOSESUM " << (dose-dosesum)/dose << " DOSE " << dose << " DOSESUM " << dosesum << G4endl;
	if( (dose-dosesum)/dose > 1.E-9 ) G4cout << ScoringVerb(debugVerb) << "GmG4PSDoseDeposit::FillScorer"<< "ERRORDOSE-DOSESUM " << (dose-dosesum)/dose << " DOSE " << dose << " DOSESUM " << dosesum << G4endl;
      }
#endif
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

