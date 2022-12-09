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
#include "G4EnergySplitter.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"
#include "G4RegularNavigationHelper.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4EmCalculator.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Step.hh"
#include "G4PVParameterised.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//
// Created: 
// 
///////////////////////////////////////////////////////////////////////////////

G4EnergySplitter::G4EnergySplitter()
{
   theElossExt = new G4EnergyLossForExtrapolator(0);
   thePhantomParam = 0;
   theNIterations = 2;
}

G4EnergySplitter::~G4EnergySplitter()
{
	delete theElossExt;
}

G4int G4EnergySplitter::SplitEnergyInVolumes(const G4Step* aStep )
{
  theEnergies.clear();

  G4double edep = aStep->GetTotalEnergyDeposit();
 
#ifdef VERBOSE_ENERSPLIT
  G4bool verbose = 1;
  if( verbose ) G4cout << "G4EnergySplitter::SplitEnergyInVolumes totalEdepo " << aStep->GetTotalEnergyDeposit() 
		       << " Nsteps " << G4RegularNavigationHelper::Instance()->GetStepLengths().size() << G4endl;
#endif    
  if( G4RegularNavigationHelper::Instance()->GetStepLengths().size() == 0 ||
      aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0)  { // we are only counting dose deposit
    return (G4int)theEnergies.size();
  }
  if( G4RegularNavigationHelper::Instance()->GetStepLengths().size() == 1 ) {
    theEnergies.push_back(edep);
    return (G4int)theEnergies.size();
  }

  if( !thePhantomParam ) GetPhantomParam(TRUE);
  
  if( aStep == 0 ) return FALSE; // it is 0 when called by GmScoringMgr after last event
  
  //----- Distribute energy deposited in voxels 
  std::vector< std::pair<G4int,G4double> > rnsl = G4RegularNavigationHelper::Instance()->GetStepLengths(); 

  const G4ParticleDefinition* part = aStep->GetTrack()->GetDefinition();
  G4double kinEnergyPreOrig = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4double kinEnergyPre = kinEnergyPreOrig;
  
  G4double stepLength = aStep->GetStepLength();
  G4double slSum = 0.;
  unsigned int ii;
  for( ii = 0; ii < rnsl.size(); ++ii ){
    G4double sl = rnsl[ii].second;
    slSum += sl;
#ifdef VERBOSE_ENERSPLIT
    if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter1 step length geom " << sl << G4endl;
#endif
  }
  
#ifdef VERBOSE_ENERSPLIT
  if( verbose )
    G4cout << "G4EnergySplitter RN:  step length geom TOTAL " << slSum 
	   << " true TOTAL " << stepLength 
	   << " ratio " << stepLength/slSum 
	   << " Energy " << aStep->GetPreStepPoint()->GetKineticEnergy() 
	   << " Material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() 
	   << " Number of geom steps " << rnsl.size() << G4endl;
#endif
  //----- No iterations to correct elost and msc => distribute energy deposited according to geometrical step length in each voxel
  if( theNIterations == 0 ) { 
    for( ii = 0; ii < rnsl.size(); ++ii ){
      G4double sl = rnsl[ii].second;
      G4double edepStep = edep * sl/slSum; //divide edep along steps, proportional to step length
#ifdef VERBOSE_ENERSPLIT
      if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii 
			  << " edep " << edepStep << G4endl;
#endif
	
      theEnergies.push_back(edepStep);
	
    }
  } else { //  1 or more iterations demanded
      
#ifdef VERBOSE_ENERSPLIT
    // print corrected energy at iteration 0 
    if(verbose)  {
      G4double slSum = 0.;
      for( ii = 0; ii < rnsl.size(); ++ii ){
	G4double sl = rnsl[ii].second;
	slSum += sl;
      }
      for( ii = 0; ii < rnsl.size(); ii++ ){
	G4cout  << "G4EnergySplitter::SplitEnergyInVolumes "<< ii
		<< " RN: iter0 corrected energy lost " << edep*rnsl[ii].second/slSum  
		<< G4endl;
      }
    }
#endif

    G4double slRatio = stepLength/slSum;
#ifdef VERBOSE_ENERSPLIT
    if(verbose) G4cout << "G4EnergySplitter::SplitEnergyInVolumes  RN: iter 0, step ratio " << slRatio << G4endl;
#endif
      
    //--- energy at each interaction
    G4EmCalculator emcalc;
    G4double totalELost = 0.;
    std::vector<G4double> stepLengths;
    for( G4int iiter = 1; iiter <= theNIterations; ++iiter ) {
      //--- iter1: distribute true step length in each voxel: geom SL in each voxel is multiplied by a constant so that the sum gives the total true step length
      if( iiter == 1 ) {
	for( ii = 0; ii < rnsl.size(); ++ii ){
	  G4double sl = rnsl[ii].second;
	  stepLengths.push_back( sl * slRatio );
#ifdef VERBOSE_ENERSPLIT
	  if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter" << iiter << " corrected step length " << sl*slRatio << G4endl;
#endif
	}
	
	for( ii = 0; ii < rnsl.size(); ++ii ){
	  const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	  G4double dEdx = 0.;
	  if( kinEnergyPre > 0. ) {  //t check this 
	    dEdx = emcalc.GetDEDX(kinEnergyPre, part, mate);
	  }
	  G4double elost = stepLengths[ii] * dEdx;
	  
#ifdef VERBOSE_ENERSPLIT
	  if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter1 energy lost "  << elost 
			      << " energy at interaction " << kinEnergyPre 
			      << " = stepLength " << stepLengths[ii] 
			      << " * dEdx " << dEdx << G4endl;
#endif
	  kinEnergyPre -= elost;
	  theEnergies.push_back( elost );
	  totalELost += elost;
	}
	
      } else{
	//------ 2nd and other iterations
	//----- Get step lengths corrected by changing geom2true correction
	//-- Get ratios for each energy 
	slSum = 0.;
	kinEnergyPre = kinEnergyPreOrig;
	for( ii = 0; ii < rnsl.size(); ++ii ){
	  const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	  stepLengths[ii] = theElossExt->TrueStepLength( kinEnergyPre, rnsl[ii].second , mate, part );
	  kinEnergyPre -= theEnergies[ii];
	  
#ifdef VERBOSE_ENERSPLIT
	  if(verbose) G4cout << "G4EnergySplitter::SplitEnergyInVolumes" << ii 
			     << " RN: iter" << iiter << " step length geom " << stepLengths[ii] 
			     << " geom2true " << rnsl[ii].second / stepLengths[ii]  << G4endl;
#endif
	  
	  slSum += stepLengths[ii];
	}
	
	//Correct step lengths so that they sum the total step length
	G4double slratio = aStep->GetStepLength()/slSum;
#ifdef VERBOSE_ENERSPLIT
	if(verbose) G4cout << "G4EnergySplitter::SplitEnergyInVolumes" << ii << " RN: iter" << iiter << " step ratio " << slRatio << G4endl;
#endif
	for( ii = 0; ii < rnsl.size(); ++ii ){
	  stepLengths[ii] *= slratio;
#ifdef VERBOSE_ENERSPLIT
	  if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter" << iiter << " corrected step length " << stepLengths[ii] << G4endl;
#endif
	}
	
	//---- Recalculate energy lost with this new step lengths
        kinEnergyPre = aStep->GetPreStepPoint()->GetKineticEnergy();
	totalELost = 0.;
	for( ii = 0; ii < rnsl.size(); ++ii ){
	  const G4Material* mate = thePhantomParam->GetMaterial( rnsl[ii].first );
	  G4double dEdx = 0.;
	  if( kinEnergyPre > 0. ) {
	    dEdx = emcalc.GetDEDX(kinEnergyPre, part, mate);
	  }
	  G4double elost = stepLengths[ii] * dEdx;
#ifdef VERBOSE_ENERSPLIT
	  if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter" << iiter << " energy lost " << elost 
			      << " energy at interaction " << kinEnergyPre 
			      << " = stepLength " << stepLengths[ii] 
			      << " * dEdx " << dEdx << G4endl;
#endif
	  kinEnergyPre -= elost;
	  theEnergies[ii] = elost;
	  totalELost += elost;
	}
	
      }
      
      //correct energies so that they reproduce the real step energy lost
      G4double enerRatio = (edep/totalELost);
      
#ifdef VERBOSE_ENERSPLIT
      if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes"<< ii << " RN: iter" << iiter << " energy ratio " << enerRatio << G4endl;
#endif
	
#ifdef VERBOSE_ENERSPLIT
      G4double elostTot = 0.; 
#endif
      for( ii = 0; ii < theEnergies.size(); ++ii ){
	theEnergies[ii] *= enerRatio;
#ifdef VERBOSE_ENERSPLIT
	elostTot += theEnergies[ii];
	if(verbose) G4cout  << "G4EnergySplitter::SplitEnergyInVolumes "<< ii << " RN: iter" << iiter << " corrected energy lost " << theEnergies[ii] 
			    << " orig elost " << theEnergies[ii]/enerRatio 
			    << " energy before interaction " << kinEnergyPreOrig-elostTot+theEnergies[ii]
			    << " energy after interaction " << kinEnergyPreOrig-elostTot
			    << G4endl;
#endif
      }
    }
    
  }
  
  return (G4int)theEnergies.size();
}


//-----------------------------------------------------------------------
void G4EnergySplitter::GetPhantomParam(G4bool mustExist)
{
  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  for( auto cite = pvs->cbegin(); cite != pvs->cend(); ++cite ) {
    //    G4cout << " PV " << (*cite)->GetName() << " " << (*cite)->GetTranslation() << G4endl;
    if( IsPhantomVolume( *cite ) ) {
      const G4PVParameterised* pvparam = static_cast<const G4PVParameterised*>(*cite);
      G4VPVParameterisation* param = pvparam->GetParameterisation();
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //      G4cout << "G4PhantomParameterisation volume found  " << (*cite)->GetName() << G4endl;
      thePhantomParam = static_cast<G4PhantomParameterisation*>(param);
    }
  }
  
  if( !thePhantomParam && mustExist )
    G4Exception("G4EnergySplitter::GetPhantomParam",
                "PhantomParamError", FatalException,
                "No G4PhantomParameterisation found !");
}


//-----------------------------------------------------------------------
G4bool G4EnergySplitter::IsPhantomVolume( G4VPhysicalVolume* pv )
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

//-----------------------------------------------------------------------
void G4EnergySplitter::GetLastVoxelID( G4int& voxelID)
{	
  voxelID = (*(G4RegularNavigationHelper::Instance()->GetStepLengths().cbegin())).first;
}

//-----------------------------------------------------------------------
void G4EnergySplitter::GetFirstVoxelID( G4int& voxelID)
{
  voxelID =  (*(G4RegularNavigationHelper::Instance()->GetStepLengths().crbegin())).first;
}

//-----------------------------------------------------------------------
void G4EnergySplitter::GetVoxelID( G4int stepNo, G4int& voxelID )
{
  if( stepNo < 0 || stepNo >= G4int(G4RegularNavigationHelper::Instance()->GetStepLengths().size()) ) {
  G4Exception("G4EnergySplitter::GetVoxelID",
	      "Invalid stepNo, smaller than 0 or bigger or equal to number of voxels traversed",
	      FatalErrorInArgument,
	      G4String("stepNo = " + G4UIcommand::ConvertToString(stepNo) + ", number of voxels = " + G4UIcommand::ConvertToString(G4int(G4RegularNavigationHelper::Instance()->GetStepLengths().size())) ).c_str());
  }
  std::vector< std::pair<G4int,G4double> >::const_iterator ite = G4RegularNavigationHelper::Instance()->GetStepLengths().cbegin();
  advance( ite, stepNo );
  voxelID = (*ite).first;

}


//-----------------------------------------------------------------------
void G4EnergySplitter::GetStepLength( G4int stepNo, G4double& stepLength )
{
  std::vector< std::pair<G4int,G4double> >::const_iterator ite = G4RegularNavigationHelper::Instance()->GetStepLengths().cbegin();
  advance( ite, stepNo );
  stepLength = (*ite).second;
}
