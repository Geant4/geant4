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
// MODULE:              G4RadioactiveDecay.cc
//
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// Documentation avaialable at http://www.space.dera.gov.uk/space_env/rdm.html
//   These include:
//       User Requirement Document (URD)
//       Software Specification Documents (SSD)
//       Software User Manual (SUM)
//       Technical Note (TN) on the physics and algorithms
//
//    The test and example programs are not included in the public release of 
//    G4 but they can be downloaded from
//      http://www.space.qinetiq.com/space_env/rdm.html
// 
// CHANGE HISTORY
// --------------
//
// 13 Oct  2015, L.G. Sarmiento Neutron emission added
//
// 06 Aug  2014, L.G. Sarmiento Proton decay mode added mimicking the alpha decay
//
// 03 Oct  2012, V. Ivanchenko removed internal table for mean free path 
//                             similar to what is done for as G4Decay
// 10 July 2012, L. Desorgher
//			-In LoadDecayTable:
//                          Add LoadedNuclei.push_back(theParentNucleus.GetParticleName());
//			also for the case where user data files are used. Correction for bug
//			1324. Changes proposed by Joa L.
//
//
// 01 May 2012, L. Desorgher
//			-Force the reading of user file to theIsotopeTable
// 			-Merge the development by Fan Lei for activation computation
//
// 17 Oct 2011, L. Desorgher
//			-Add possibility for the user to load its own decay file.
//			-Set halflifethreshold negative by default to allow the tracking of all
//			   excited nuclei resulting from a radioactive decay
//
// 01 June 2011, M. Kelsey -- Add directional biasing interface to allow for
//		"collimation" of decay daughters.
// 16 February 2006, V.Ivanchenko fix problem in IsApplicable connected with
//            8.0 particle design
// 18 October 2002, F. Lei
//            in the case of beta decay, added a check of the end-energy 
//            to ensure it is > 0.
//            ENSDF occationally have beta decay entries with zero energies
//
// 27 Sepetember 2001, F. Lei
//            verboselevel(0) used in constructor
//
// 01 November 2000, F.Lei
//            added " ee = e0 +1. ;" as line 763
//            tagged as "radiative_decay-V02-00-02"              
// 28 October 2000, F Lei 
//            added fast beta decay mode. Many files have been changed.
//            tagged as "radiative_decay-V02-00-01"
//
// 25 October 2000, F Lei, DERA UK
//            1) line 1185 added 'const' to work with tag "Track-V02-00-00"
//            tagged as "radiative_decay-V02-00-00"
// 14 April 2000, F Lei, DERA UK
// 0.b.4 release. Changes are:
//            1) Use PhotonEvaporation instead of DiscreteGammaDeexcitation
//            2) VR: Significant efficiency improvement
// 
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecaymessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4ITDecay.hh"
#include "G4BetaDecayType.hh"
#include "G4BetaMinusDecay.hh"
#include "G4BetaPlusDecay.hh"
#include "G4ECDecay.hh"
#include "G4AlphaDecay.hh"
#include "G4ProtonDecay.hh"
#include "G4NeutronDecay.hh"
#include "G4VDecayChannel.hh"
#include "G4NuclearDecay.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4Fragment.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4BetaDecayType.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4LevelManager.hh"
#include "G4ThreeVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"

#include "G4HadronicProcessType.hh"
#include "G4HadronicException.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhotonEvaporation.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace CLHEP;

const G4double G4RadioactiveDecay::levelTolerance = 10.0*eV;
const G4ThreeVector G4RadioactiveDecay::origin(0.,0.,0.);

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
G4Mutex G4RadioactiveDecay::radioactiveDecayMutex = G4MUTEX_INITIALIZER;
DecayTableMap* G4RadioactiveDecay::master_dkmap = 0;
#endif

G4RadioactiveDecay::G4RadioactiveDecay(const G4String& processName)
 : G4VRestDiscreteProcess(processName, fDecay), isInitialised(false),
   forceDecayDirection(0.,0.,0.), forceDecayHalfAngle(0.*deg), verboseLevel(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4RadioactiveDecay constructor: processName = " << processName
           << G4endl;
  }
#endif

  SetProcessSubType(fRadioactiveDecay);

  theRadioactiveDecaymessenger = new G4RadioactiveDecaymessenger(this);
  pParticleChange = &fParticleChangeForRadDecay;

  // Set up photon evaporation for use in G4ITDecay
  photonEvaporation = new G4PhotonEvaporation();
  // photonEvaporation->SetVerboseLevel(2);
  photonEvaporation->RDMForced(true);
  photonEvaporation->SetICM(true);

  // Reset the list of user defined data files
  theUserRadioactiveDataFiles.clear();

  // Instantiate the map of decay tables
#ifdef G4MULTITHREADED
  G4AutoLock lk(&G4RadioactiveDecay::radioactiveDecayMutex);
  if(!master_dkmap) master_dkmap = new DecayTableMap;
#endif
  dkmap = new DecayTableMap;

  // Apply default values.
  NSourceBin  = 1;
  SBin[0]     = 0.* s;
  SBin[1]     = 1.* s;
  SProfile[0] = 1.;
  SProfile[1] = 0.;
  NDecayBin   = 1;
  DBin[0]     = 0. * s ;
  DBin[1]     = 1. * s;
  DProfile[0] = 1.;
  DProfile[1] = 0.;
  decayWindows[0] = 0;
  G4RadioactivityTable* rTable = new G4RadioactivityTable() ;
  theRadioactivityTables.push_back(rTable);
  NSplit      = 1;
  AnalogueMC  = true ;
  FBeta       = false ;
  BRBias      = true ;
  applyICM    = true ;
  applyARM    = true ;
  halflifethreshold = nanosecond;

  // RDM applies to all logical volumes by default
  isAllVolumesMode = true;
  SelectAllVolumes();
}


G4RadioactiveDecay::~G4RadioactiveDecay()
{
  delete theRadioactiveDecaymessenger;
  delete photonEvaporation;
  for (DecayTableMap::iterator i = dkmap->begin(); i != dkmap->end(); i++) {
    delete i->second;
  }
  dkmap->clear();
  delete dkmap;
}


G4bool G4RadioactiveDecay::IsApplicable(const G4ParticleDefinition& aParticle)
{
  // All particles other than G4Ions, are rejected by default
  if (((const G4Ions*)(&aParticle))->GetExcitationEnergy() > 0.) {return true;}
  if (aParticle.GetParticleName() == "GenericIon") {
    return true;
  } else if (!(aParticle.GetParticleType() == "nucleus")
             || aParticle.GetPDGLifeTime() < 0. ) {
    return false;
  }

  // Determine whether the nuclide falls into the correct A and Z range
  G4int A = ((const G4Ions*) (&aParticle))->GetAtomicMass();
  G4int Z = ((const G4Ions*) (&aParticle))->GetAtomicNumber();

  if (A > theNucleusLimits.GetAMax() || A < theNucleusLimits.GetAMin())
    {return false;}
  else if (Z > theNucleusLimits.GetZMax() || Z < theNucleusLimits.GetZMin())
    {return false;}
  return true;
}

G4DecayTable* G4RadioactiveDecay::GetDecayTable(const G4ParticleDefinition* aNucleus)
{
  G4String key = aNucleus->GetParticleName();
  DecayTableMap::iterator table_ptr = dkmap->find(key);

  G4DecayTable* theDecayTable = 0;
  if (table_ptr == dkmap->end() ) {                   // If table not there,     
    theDecayTable = LoadDecayTable(*aNucleus);        // load from file and
    if(theDecayTable) (*dkmap)[key] = theDecayTable;  // store in library 
  } else {
    theDecayTable = table_ptr->second;
  }

  return theDecayTable;
}


void G4RadioactiveDecay::SelectAVolume(const G4String aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++) {
    volume=(*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      ValidVolumes.push_back(aVolume);
      std::sort(ValidVolumes.begin(), ValidVolumes.end());
      // sort need for performing binary_search
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	G4cout << " RDM Applies to : " << aVolume << G4endl; 
#endif
    } else if(i == theLogicalVolumes->size()) {
      G4cerr << "SelectAVolume: "<< aVolume
             << " is not a valid logical volume name" << G4endl; 
    }
  }
}


void G4RadioactiveDecay::DeselectAVolume(const G4String aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume=(*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      std::vector<G4String>::iterator location;
      location = std::find(ValidVolumes.begin(),ValidVolumes.end(),aVolume);
      if (location != ValidVolumes.end()) {
        ValidVolumes.erase(location);
        std::sort(ValidVolumes.begin(), ValidVolumes.end());
        isAllVolumesMode =false;
      } else {
        G4cerr << " DeselectVolume:" << aVolume << " is not in the list "
               << G4endl; 
      }	  
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 0)
        G4cout << " DeselectVolume: " << aVolume << " is removed from list "
               << G4endl; 
#endif
    } else if (i ==  theLogicalVolumes->size()) {
      G4cerr << " DeselectVolume:" << aVolume
             << "is not a valid logical volume name" << G4endl; 
    }
  }
}


void G4RadioactiveDecay::SelectAllVolumes() 
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  ValidVolumes.clear();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout << " RDM Applies to all Volumes"  << G4endl;
#endif
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume = (*theLogicalVolumes)[i];
    ValidVolumes.push_back(volume->GetName());    
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout << "       RDM Applies to Volume " << volume->GetName() << G4endl;
#endif
  }
  std::sort(ValidVolumes.begin(), ValidVolumes.end());
  // sort needed in order to allow binary_search
  isAllVolumesMode=true;
}


void G4RadioactiveDecay::DeselectAllVolumes() 
{
  ValidVolumes.clear();
  isAllVolumesMode=false;
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "RDM removed from all volumes" << G4endl; 
#endif
}


G4bool
G4RadioactiveDecay::IsRateTableReady(const G4ParticleDefinition& aParticle)
{
  // Check whether the radioactive decay rates table for the ion has already
  // been calculated.
  G4String aParticleName = aParticle.GetParticleName();
  for (size_t i = 0; i < theDecayRateTableVector.size(); i++) {
    if (theDecayRateTableVector[i].GetIonName() == aParticleName) return true;
  }
  return false;
}

// GetDecayRateTable
// retrieve the decayratetable for the specified aParticle

void
G4RadioactiveDecay::GetDecayRateTable(const G4ParticleDefinition& aParticle)
{
  G4String aParticleName = aParticle.GetParticleName();

  for (size_t i = 0; i < theDecayRateTableVector.size(); i++) {
    if (theDecayRateTableVector[i].GetIonName() == aParticleName) {
      theDecayRateVector = theDecayRateTableVector[i].GetItsRates();
    }
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) {
    G4cout << "The DecayRate Table for " << aParticleName << " is selected."
           <<  G4endl;
  }
#endif
}

// ConvolveSourceTimeProfile performs the convolution of the source time profile
// function with a single exponential characterized by a decay constant in the 
// decay chain.  The time profile is treated as a step function so that the 
// convolution integral can be done bin-by-bin.
// Input time and mean life (tau) are in ns.   

G4double
G4RadioactiveDecay::ConvolveSourceTimeProfile(const G4double t, const G4double tau)
{
  long double convolvedTime = 0.L;
  G4int nbin;
  if ( t > SBin[NSourceBin]) {
    nbin  = NSourceBin;
  } else {
    nbin = 0;

    G4int loop = 0;
    G4ExceptionDescription ed;
    ed << " While count exceeded " << G4endl;
    while (t > SBin[nbin]) {  /* Loop checking, 01.09.2015, D.Wright */
      loop++;
      if (loop > 1000) {
        G4Exception("G4RadioactiveDecay::ConvolveSourceTimeProfile()",
                    "HAD_RDM_100", JustWarning, ed);
        break;
      }
 
      nbin++;
    }
    nbin--;
  }
  long double lt = t ;
  long double ltau = tau;
  // G4cout << " Convolve: tau = " << tau << G4endl;
  if (nbin > 0) {
    for (G4int i = 0; i < nbin; i++) {
      convolvedTime += (long double)SProfile[i] *
       (std::exp(-(lt-(long double)SBin[i+1])/ltau)-std::exp(-(lt-(long double)SBin[i])/ltau));
    }
  }
  convolvedTime +=  (long double)SProfile[nbin] * (1.L-std::exp(-(lt-(long double)SBin[nbin])/ltau));
  // Is the above line necessary?  If so, the 1.L looks incorrect - should be an exp
  // Also, it looks like the final integral should be multiplied by ltau

  if (convolvedTime < 0.)  {
    G4cout << " Convolved time =: " << convolvedTime << " reset to zero! " << G4endl;
    G4cout << " t = " << t << " tau = " << tau << G4endl;
    G4cout << SBin[nbin] << " " << SBin[0] << G4endl;
    convolvedTime = 0.;
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << " Convolved time: " << convolvedTime << G4endl;
#endif
  return (G4double)convolvedTime ;
}

/*
// Other implementation tests to avoid use of  long double
G4double G4RadioactiveDecay::GetTaoTime(const G4double t, const G4double tao)
{
  long double taotime =0.L;
  G4int nbin;
  if ( t > SBin[NSourceBin]) {
    nbin  = NSourceBin;}
  else {
    nbin = 0;
    while (t > SBin[nbin]) nbin++;
    nbin--;}
  long double lt = t ;
  long double ltao = tao;
  long double factor,factor1,dt1,dt;
  if (nbin > 0) {
	    for (G4int i = 0; i < nbin; i++)
	      { long double s1=SBin[i];
	        long double s2=SBin[i+1];
	    	dt1=(s2-s1)/ltao;
	      	if (dt1 <50.) {
	      		factor1=std::exp(dt1)-1.;
	      		if (factor1<dt1) factor1 =dt1;
	      		dt=(lt-s1)/ltao;
	      		factor=std::exp(-dt);
	      	}
	      	else {
	      		factor1=1.-std::exp(-dt1);
	      		dt=(lt-s2)/ltao;
	      		factor=std::exp(-dt);
	      	}
	      	G4cout<<(long double) SProfile[i] *factor*factor1<<'\t'<<std::endl;
	      	long double test  = (long double)SProfile[i] * (std::exp(-(lt-(long double)SBin[i+1])/ltao)-std::exp(-(lt-(long double)SBin[i])/ltao));
	      	G4cout<<test<<std::endl;
	      	taotime += (long double) SProfile[i] *factor*factor1;
	   }
   }
  long double s=SBin[nbin];
  dt1=(lt-s)/ltao;
  factor=1.-std::exp(-dt1);
  taotime += (long double)  SProfile[nbin] *factor;
 if (taotime < 0.)  {
    G4cout <<" Tao time =: " <<taotime << " reset to zero!"<<G4endl;
    G4cout <<" t = " << t <<" tao = " <<tao <<G4endl;
    G4cout << SBin[nbin] << " " <<SBin[0] << G4endl;
    taotime = 0.;
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Tao time: " <<taotime <<G4endl;}
#endif
  return (G4double)taotime ;
}



G4double G4RadioactiveDecay::GetTaoTime(const G4double t, const G4double tao)
{
  G4double taotime =0.;
  G4int nbin;
  if ( t > SBin[NSourceBin]) {
    nbin  = NSourceBin;}
  else {
    nbin = 0;
    while (t > SBin[nbin]) nbin++;
    nbin--;}
  G4double lt = t ;
  G4double ltao = tao;
  G4double factor,factor1,dt1,dt;

  if (nbin > 0) {
    for (G4int i = 0; i < nbin; i++)
      { dt1=(SBin[i+1]-SBin[i])/ltao;
      	if (dt1 <50.) {
      		factor1=std::exp(dt1)-1.;
      		if (factor1<dt1) factor1 =dt1;
      		dt=(lt-SBin[i])/ltao;
      		factor=std::exp(-(lt-SBin[i])/ltao);
      		G4cout<<factor<<'\t'<<factor1<<std::endl;
      	}
      	else {
      		factor1=1.-std::exp(-dt1);
      		factor=std::exp(-(lt-SBin[i+1])/ltao);
      	}
      	G4cout<<factor<<'\t'<<factor1<<std::endl;
      	taotime += SProfile[i] *factor*factor1;
      	G4cout<<taotime<<std::endl;
      }
  }
  dt1=(lt-SBin[nbin])/ltao;
  factor=1.-std::exp(-dt1);
  if (factor<(dt1-0.5*dt1*dt1)) factor =dt1-0.5*dt1*dt1;


  taotime +=  SProfile[nbin] *factor;
  G4cout<<factor<<'\t'<<taotime<<std::endl;
  if (taotime < 0.)  {
    G4cout <<" Tao time =: " <<taotime << " reset to zero!"<<G4endl;
    G4cout <<" t = " << t <<" tao = " <<tao <<G4endl;
    G4cout << SBin[nbin] << " " <<SBin[0] << G4endl;
    taotime = 0.;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Tao time: " <<taotime <<G4endl;}
#endif
  return (G4double)taotime ;
}
*/
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetDecayTime                                                              //
//    Randomly select a decay time for the decay process, following the       //
//    supplied decay time bias scheme.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecay::GetDecayTime()
{
  G4double decaytime = 0.;
  G4double rand = G4UniformRand();
  G4int i = 0;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
  while ( DProfile[i] < rand) {  /* Loop checking, 01.09.2015, D.Wright */
    i++;
    loop++;
    if (loop > 100000) {
      G4Exception("G4RadioactiveDecay::GetDecayTime()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
  }

  rand = G4UniformRand();
  decaytime = DBin[i] + rand*(DBin[i+1]-DBin[i]);
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout <<" Decay time: " <<decaytime/s <<"[s]" <<G4endl;
#endif
  return  decaytime;	    
}


G4int G4RadioactiveDecay::GetDecayTimeBin(const G4double aDecayTime)
{
  G4int i = 0;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
  while ( aDecayTime > DBin[i] ) {   /* Loop checking, 01.09.2015, D.Wright */
    i++;
    loop++;
    if (loop > 100000) {
      G4Exception("G4RadioactiveDecay::GetDecayTimeBin()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
  }

  return  i;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanLifeTime (required by the base class)                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecay::GetMeanLifeTime(const G4Track& theTrack,
                                             G4ForceCondition*)
{
  // For variance reduction the time is set to 0 so as to force the particle
  // to decay immediately.
  // In analogueMC mode it returns the particle's mean-life.

  G4double meanlife = 0.;
  if (AnalogueMC) {
    const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
    const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();
    G4double theLife = theParticleDef->GetPDGLifeTime();
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 2) {
       G4cout << "G4RadioactiveDecay::GetMeanLifeTime() " << G4endl;
       G4cout << "KineticEnergy: " << theParticle->GetKineticEnergy()/GeV
              << " GeV, Mass: " << theParticle->GetMass()/GeV
              << " GeV, Life time: " << theLife/ns << " ns " << G4endl;
    }
#endif
    if (theParticleDef->GetPDGStable()) {meanlife = DBL_MAX;}
    else if (theLife < 0.0) {meanlife = DBL_MAX;}
    else {meanlife = theLife;}
    // Set meanlife to zero for excited istopes which are not in the
    // RDM database
    if (((const G4Ions*)(theParticleDef))->GetExcitationEnergy() > 0. &&
                                          meanlife == DBL_MAX) {meanlife = 0.;}
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << " mean life time: " << meanlife/s << " s " << G4endl;
#endif

  return  meanlife;
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  GetMeanFreePath for decay in flight                               //
//                                                                    //
////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecay::GetMeanFreePath (const G4Track& aTrack, G4double,
                                              G4ForceCondition*)
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
  G4double tau = aParticleDef->GetPDGLifeTime();
  G4double aMass = aParticle->GetMass();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "G4RadioactiveDecay::GetMeanFreePath() " << G4endl;
    G4cout << "  KineticEnergy: " << aParticle->GetKineticEnergy()/GeV
           << " GeV, Mass: " << aMass/GeV << " GeV, tau: " << tau << " ns "
           << G4endl;
  }
#endif
  G4double pathlength = DBL_MAX;
  if (tau != -1) {
    // Ion can decay

    if (tau < -1000.0) {
      pathlength = DBL_MIN;  // nuclide had very short lifetime or wasn't in table

    } else if (tau < 0.0) {
      G4cout << aParticleDef->GetParticleName() << " has lifetime " << tau << G4endl;
      G4ExceptionDescription ed;
      ed << "Ion has negative lifetime " << tau
         << " but is not stable.  Setting mean free path to DBL_MAX" << G4endl; 
      G4Exception("G4RadioactiveDecay::GetMeanFreePath()", "HAD_RDM_011",
                   JustWarning, ed);
      pathlength = DBL_MAX;

    } else {
      // Calculate mean free path
      G4double betaGamma = aParticle->GetTotalMomentum()/aMass;
      pathlength = c_light*tau*betaGamma;

      if (pathlength < DBL_MIN) {
        pathlength = DBL_MIN;
#ifdef G4VERBOSE
        if (GetVerboseLevel() > 2) {
          G4cout << "G4Decay::GetMeanFreePath: "
                 << aParticleDef->GetParticleName()
                 << " stops, kinetic energy = "
                 << aParticle->GetKineticEnergy()/keV <<" keV " << G4endl;
        }
#endif
      }
    }
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "mean free path: "<< pathlength/m << " m" << G4endl;
  }
#endif
  return  pathlength;
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  BuildPhysicsTable - initialization of atomic de-excitation        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

void G4RadioactiveDecay::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (!isInitialised) {
    isInitialised = true;
    G4LossTableManager* theManager = G4LossTableManager::Instance();
    G4VAtomDeexcitation* p = theManager->AtomDeexcitation();
    if (!p) {
      G4ExceptionDescription ed;
      ed << " Atomic deexcitation is not defined."; 
      G4Exception("G4RadioactiveDecay::BuildPhysicsTable", "HAD_RDM_001",
                  FatalException, ed);
      /*
      p = new G4UAtomicDeexcitation();
      p->SetFluo(true);
      p->SetAuger(true);
      p->InitialiseAtomicDeexcitation();
      theManager->SetAtomDeexcitation(p);
      */
    }

    G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
    param->SetUseFilesNEW(true);
    param->SetCorrelatedGamma(true);
  }
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  LoadDecayTable loads the decay scheme from the RadioactiveDecay database  // 
//  for the parent nucleus.                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4DecayTable*
G4RadioactiveDecay::LoadDecayTable(const G4ParticleDefinition& theParentNucleus)
{
  // Generate input data file name using Z and A of the parent nucleus
  // file containing radioactive decay data.
  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();

  G4double levelEnergy = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4Ions::G4FloatLevelBase floatingLevel =
    ((const G4Ions*)(&theParentNucleus))->GetFloatLevelBase();

#ifdef G4MULTITHREADED
  G4AutoLock lk(&G4RadioactiveDecay::radioactiveDecayMutex);

  G4String key = theParentNucleus.GetParticleName();
  DecayTableMap::iterator master_table_ptr = master_dkmap->find(key);

  if (master_table_ptr != master_dkmap->end() ) {   // If table is there              
    return master_table_ptr->second;
  }
#endif

  //Check if data have been provided by the user
  G4String file = theUserRadioactiveDataFiles[1000*A+Z];

  if (file == "") {
    if (!getenv("G4RADIOACTIVEDATA") ) {
      G4cout << "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files."
             << G4endl;
      throw G4HadronicException(__FILE__, __LINE__, " Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files.");
    }
    G4String dirName = getenv("G4RADIOACTIVEDATA");

    std::ostringstream os;
    os << dirName << "/z" << Z << ".a" << A << '\0';
    file = os.str();
  }

  G4DecayTable* theDecayTable = new G4DecayTable();
  G4bool found(false);     // True if energy level matches one in table

  std::ifstream DecaySchemeFile;
  DecaySchemeFile.open(file);

  if (DecaySchemeFile.good()) {
    // Initialize variables used for reading in radioactive decay data
    G4bool floatMatch(false);
    const G4int nMode = 9;
    G4double modeTotalBR[nMode] = {0.0};
    G4double modeSumBR[nMode];
    for (G4int i = 0; i < nMode; i++) {
      modeSumBR[i] = 0.0;
    }

    char inputChars[120]={' '};
    G4String inputLine;
    G4String recordType("");
    G4String floatingFlag("");
    G4String daughterFloatFlag("");
    G4Ions::G4FloatLevelBase daughterFloatLevel;
    G4RadioactiveDecayMode theDecayMode;
    G4double decayModeTotal(0.0);
    G4double parentExcitation(0.0);
    G4double a(0.0);
    G4double b(0.0);
    G4double c(0.0);
    G4double dummy(0.0);
    G4BetaDecayType betaType(allowed);

    // Loop through each data file record until you identify the decay
    // data relating to the nuclide of concern.

    G4bool complete(false);  // bool insures only one set of values read for any
                             // given parent energy level
    G4int loop = 0;
    G4ExceptionDescription ed;
    ed << " While count exceeded " << G4endl;
 
    while (!complete && !DecaySchemeFile.getline(inputChars, 120).eof()) {  /* Loop checking, 01.09.2015, D.Wright */
      loop++;
      if (loop > 100000) {
        G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_100", JustWarning, ed);
        break;
      }
 
      inputLine = inputChars;
      inputLine = inputLine.strip(1);
      if (inputChars[0] != '#' && inputLine.length() != 0) {
        std::istringstream tmpStream(inputLine);

        if (inputChars[0] == 'P') {
          // Nucleus is a parent type.  Check excitation level to see if it
          // matches that of theParentNucleus
          tmpStream >> recordType >> parentExcitation >> floatingFlag >> dummy;
          // "dummy" takes the place of half-life
          //  Now read in from ENSDFSTATE in particle category

          if (found) {
            complete = true;
          } else {
            // Take first level which matches excitation energy regardless of floating level
            found = (std::abs(parentExcitation*keV - levelEnergy) < levelTolerance);
            if (floatingLevel != noFloat) {
              // If floating level specificed, require match of both energy and floating level
              floatMatch = (floatingLevel == G4Ions::FloatLevelBase(floatingFlag.back()) );
              if (!floatMatch) found = false;
            }
          }

        } else if (found) {
          // The right part of the radioactive decay data file has been found.  Search
          // through it to determine the mode of decay of the subsequent records.

          // Store for later the total decay probability for each decay mode 
          if (inputLine.length() < 72) {
            tmpStream >> theDecayMode >> dummy >> decayModeTotal;
            switch (theDecayMode) {
              case IT:
                {
                G4ITDecay* anITChannel = new G4ITDecay(&theParentNucleus, decayModeTotal,
                                                       0.0, 0.0, photonEvaporation);
                anITChannel->SetHLThreshold(halflifethreshold);
                anITChannel->SetARM(applyARM);
                theDecayTable->Insert(anITChannel);
//                anITChannel->DumpNuclearInfo();
                }
                break;
              case BetaMinus:
                modeTotalBR[1] = decayModeTotal; break;
              case BetaPlus:
                modeTotalBR[2] = decayModeTotal; break;
              case KshellEC:
                modeTotalBR[3] = decayModeTotal; break;
              case LshellEC:
                modeTotalBR[4] = decayModeTotal; break;
              case MshellEC:
                modeTotalBR[5] = decayModeTotal; break;
              case Alpha:
                modeTotalBR[6] = decayModeTotal; break;
              case Proton:
                modeTotalBR[7] = decayModeTotal; break;
              case Neutron:
                modeTotalBR[8] = decayModeTotal; break;
              case BDProton:
                break;
              case BDNeutron:
                break;
              case Beta2Minus:
                break;
              case Beta2Plus:
                break;
              case Proton2:
                break;
              case Neutron2:
                break;
              case SpFission:
                break;
              case RDM_ERROR:

              default:
                G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_000",
                            FatalException, "Selected decay mode does not exist");
            }  // switch

          } else {
            if (inputLine.length() < 84) {
              tmpStream >> theDecayMode >> a >> daughterFloatFlag >> b >> c;
              betaType = allowed;
            } else {
              tmpStream >> theDecayMode >> a >> daughterFloatFlag >> b >> c >> betaType;
            }

            // Allowed transitions are the default. Forbidden transitions are
            // indicated in the last column.
            a /= 1000.;
            c /= 1000.;
            daughterFloatLevel = G4Ions::FloatLevelBase(daughterFloatFlag.back());

            switch (theDecayMode) {
              case BetaMinus:
              {
                G4BetaMinusDecay* aBetaMinusChannel =
                  new G4BetaMinusDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                       daughterFloatLevel, betaType);
//              aBetaMinusChannel->DumpNuclearInfo();
                aBetaMinusChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aBetaMinusChannel);
                modeSumBR[1] += b;
              }
              break;

              case BetaPlus:
              {
                G4BetaPlusDecay* aBetaPlusChannel =
                  new G4BetaPlusDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                      daughterFloatLevel, betaType);
//              aBetaPlusChannel->DumpNuclearInfo();
                aBetaPlusChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aBetaPlusChannel);
                modeSumBR[2] += b;
              }
              break;

              case KshellEC:  // K-shell electron capture
              {
                G4ECDecay* aKECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, KshellEC);
//              aKECChannel->DumpNuclearInfo();
                aKECChannel->SetHLThreshold(halflifethreshold);
                aKECChannel->SetARM(applyARM);
                theDecayTable->Insert(aKECChannel);
                modeSumBR[3] += b;
              }
              break;

              case LshellEC:  // L-shell electron capture
              {
                G4ECDecay* aLECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, LshellEC);
//              aLECChannel->DumpNuclearInfo();
                aLECChannel->SetHLThreshold(halflifethreshold);
                aLECChannel->SetARM(applyARM);
                theDecayTable->Insert(aLECChannel);
                modeSumBR[4] += b;
              }
              break;

              case MshellEC:  // M-shell electron capture
              {
                G4ECDecay* aMECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, MshellEC);
//              aMECChannel->DumpNuclearInfo();
                aMECChannel->SetHLThreshold(halflifethreshold);
                aMECChannel->SetARM(applyARM);
                theDecayTable->Insert(aMECChannel);
                modeSumBR[5] += b;
              }
              break;

              case Alpha:
              {
                G4AlphaDecay* anAlphaChannel =
                  new G4AlphaDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                   daughterFloatLevel);
//              anAlphaChannel->DumpNuclearInfo();
                anAlphaChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(anAlphaChannel);
                modeSumBR[6] += b;
              }
              break;

	      case Proton:
              {
                G4ProtonDecay* aProtonChannel =
                  new G4ProtonDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                    daughterFloatLevel);
//              aProtonChannel->DumpNuclearInfo();
                aProtonChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aProtonChannel);
                modeSumBR[7] += b;
              }
              break;

              case Neutron:
              {
                G4NeutronDecay* aNeutronChannel =
                  new G4NeutronDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                     daughterFloatLevel);
//              aNeutronChannel->DumpNuclearInfo();
                aNeutronChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aNeutronChannel);
                modeSumBR[8] += b;
              }
              break;

              case BDProton:
                  // Not yet implemented
                  // G4cout << " beta-delayed proton decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case BDNeutron:
                  // Not yet implemented
                  // G4cout << " beta-delayed neutron decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Beta2Minus:
                  // Not yet implemented
                  // G4cout << " Double beta- decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Beta2Plus:
                  // Not yet implemented
                  // G4cout << " Double beta+ decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Proton2:
                  // Not yet implemented
                  // G4cout << " Double proton decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Neutron2:
                  // Not yet implemented
                  // G4cout << " Double beta- decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case SpFission:
            	  // Not yet implemented
            	  //G4cout<<"Sp fission channel"<<a<<'\t'<<b<<'\t'<<c<<std::endl;
            	  break;
              case RDM_ERROR:

              default:
                G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_000",
                            FatalException, "Selected decay mode does not exist");
            }  // switch
          }  // line < 72
        }  // if char == P
      }  // if char != #
    }  // While

    // Go through the decay table and make sure that the branching ratios are
    // correctly normalised.

    G4VDecayChannel* theChannel = 0;
    G4NuclearDecay* theNuclearDecayChannel = 0;
    G4String mode = "";

    G4double theBR = 0.0;
    for (G4int i = 0; i < theDecayTable->entries(); i++) {
      theChannel = theDecayTable->GetDecayChannel(i);
      theNuclearDecayChannel = static_cast<G4NuclearDecay*>(theChannel);
      theDecayMode = theNuclearDecayChannel->GetDecayMode();

      if (theDecayMode != IT) {
	theBR = theChannel->GetBR();
	theChannel->SetBR(theBR*modeTotalBR[theDecayMode]/modeSumBR[theDecayMode]);
      }
    }
  }  // decay file exists

  DecaySchemeFile.close();

  if (!found && levelEnergy > 0) {
    // Case where IT cascade for excited isotopes has no entries in RDM database
    // Decay mode is isomeric transition.
    G4ITDecay* anITChannel = new G4ITDecay(&theParentNucleus, 1.0, 0.0, 0.0,
                                           photonEvaporation);
    anITChannel->SetHLThreshold(halflifethreshold);
    anITChannel->SetARM(applyARM);
    theDecayTable->Insert(anITChannel);
  }

  if (theDecayTable && GetVerboseLevel() > 1) {
    theDecayTable->DumpInfo();
  }

#ifdef G4MULTITHREADED
  //(*master_dkmap)[key] = theDecayTable;                  // store in master library 
#endif
  return theDecayTable;
}

void
G4RadioactiveDecay::AddUserDecayDataFile(G4int Z, G4int A, G4String filename)
{
  if (Z < 1 || A < 2) G4cout << "Z and A not valid!" << G4endl;

  std::ifstream DecaySchemeFile(filename);
  if (DecaySchemeFile) {
    G4int ID_ion = A*1000 + Z;
    theUserRadioactiveDataFiles[ID_ion] = filename;
  } else {
    G4cout << "The file " << filename << " does not exist!" << G4endl;
  }
}


void
G4RadioactiveDecay::SetDecayRate(G4int theZ, G4int theA, G4double theE, 
                                 G4int theG, std::vector<G4double> theCoefficients, 
                                 std::vector<G4double> theTaos)
//  Why not make this a method of G4RadioactiveDecayRate? (e.g. SetParameters)
{ 
  //fill the decay rate vector 
  theDecayRate.SetZ(theZ);
  theDecayRate.SetA(theA);
  theDecayRate.SetE(theE);
  theDecayRate.SetGeneration(theG);
  theDecayRate.SetDecayRateC(theCoefficients);
  theDecayRate.SetTaos(theTaos);
}


void
G4RadioactiveDecay::AddDecayRateTable(const G4ParticleDefinition& theParentNucleus)
{
  // Use extended Bateman equation to calculate the radioactivities of all
  // progeny of theParentNucleus.  The coefficients required to do this are 
  // calculated using the method of P. Truscott (Ph.D. thesis and
  // DERA Technical Note DERA/CIS/CIS2/7/36/4/10) 11 January 2000.
  // Coefficients are then added to the decay rate table vector

  // Create and initialise variables used in the method.
  theDecayRateVector.clear();

  G4int nGeneration = 0;

  std::vector<G4double> taos;

  // Dimensionless A coefficients of Eqs. 4.24 and 4.25 of the TN
  std::vector<G4double> Acoeffs;

  // According to Eq. 4.26 the first coefficient (A_1:1) is -1
  Acoeffs.push_back(-1.);

  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();
  G4double E = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4double tao = theParentNucleus.GetPDGLifeTime();
  if (tao < 0.) tao = 1e-100;
  taos.push_back(tao);
  G4int nEntry = 0;

  // Fill the decay rate container (G4RadioactiveDecayRate) with the parent 
  // isotope data
  SetDecayRate(Z,A,E,nGeneration,Acoeffs,taos);   // Fill TP with parent lifetime

  // store the decay rate in decay rate vector
  theDecayRateVector.push_back(theDecayRate);
  nEntry++;

  // Now start treating the secondary generations.
  G4bool stable = false;
  G4int i;
  G4int j;
  G4VDecayChannel* theChannel = 0;
  G4NuclearDecay* theNuclearDecayChannel = 0;

  G4ITDecay* theITChannel = 0;
  G4BetaMinusDecay* theBetaMinusChannel = 0;
  G4BetaPlusDecay* theBetaPlusChannel = 0;
  G4AlphaDecay* theAlphaChannel = 0;
  G4ProtonDecay* theProtonChannel = 0;
  G4NeutronDecay* theNeutronChannel = 0;
  G4RadioactiveDecayMode theDecayMode;
  G4double theBR = 0.0;
  G4int AP = 0;
  G4int ZP = 0;
  G4int AD = 0;
  G4int ZD = 0;
  G4double EP = 0.;
  std::vector<G4double> TP;
  std::vector<G4double> RP;   // A coefficients of the previous generation
  G4ParticleDefinition *theDaughterNucleus;
  G4double daughterExcitation;
  G4double nearestEnergy = 0.0;
  G4int nearestLevelIndex = 0;
  G4ParticleDefinition *aParentNucleus;
  G4IonTable* theIonTable;
  G4DecayTable* parentDecayTable;
  G4double theRate;
  G4double TaoPlus;
  G4int nS = 0;        // Running index of first decay in a given generation
  G4int nT = nEntry;   // Total number of decays accumulated over entire history
  const G4int nMode = 9;
  G4double brs[nMode];
  //
  theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
 
  while (!stable) {   /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 10000) {
      G4Exception("G4RadioactiveDecay::AddDecayRateTable()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
    nGeneration++;
    for (j = nS; j < nT; j++) {
      // First time through, get data for parent nuclide
      ZP = theDecayRateVector[j].GetZ();
      AP = theDecayRateVector[j].GetA();
      EP = theDecayRateVector[j].GetE();
      RP = theDecayRateVector[j].GetDecayRateC();
      TP = theDecayRateVector[j].GetTaos();
      if (GetVerboseLevel() > 0) {
        G4cout << "G4RadioactiveDecay::AddDecayRateTable : daughters of ("
               << ZP << ", " << AP << ", " << EP
               << ") are being calculated,  generation = " << nGeneration
               << G4endl;
      }
//      G4cout << " Taus = " << G4endl;
//      for (G4int ii = 0; ii < TP.size(); ii++) G4cout << TP[ii] << ", " ;
//      G4cout << G4endl;

      aParentNucleus = theIonTable->GetIon(ZP,AP,EP);
      parentDecayTable = GetDecayTable(aParentNucleus);

      G4DecayTable* summedDecayTable = new G4DecayTable();
      // This instance of G4DecayTable is for accumulating BRs and decay
      // channels.  It will contain one decay channel per type of decay
      // (alpha, beta, etc.); its branching ratio will be the sum of all
      // branching ratios for that type of decay of the parent.  If the
      // halflife of a particular channel is longer than some threshold,
      // that channel will be inserted specifically and its branching
      // ratio will not be included in the above sums.
      // This instance is not used to perform actual decays.

      for (G4int k = 0; k < nMode; k++) brs[k] = 0.0;

      // Go through the decay table and sum all channels having the same decay mode
      for (i = 0; i < parentDecayTable->entries(); i++) {
        theChannel = parentDecayTable->GetDecayChannel(i);
        theNuclearDecayChannel = static_cast<G4NuclearDecay*>(theChannel);
        theDecayMode = theNuclearDecayChannel->GetDecayMode();
        daughterExcitation = theNuclearDecayChannel->GetDaughterExcitation();
        theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus() ;

        AD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
        ZD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();  
        const G4LevelManager* levelManager =
          G4NuclearLevelData::GetInstance()->GetLevelManager(ZD,AD);

        if (levelManager->NumberOfTransitions() ) {
          nearestEnergy = levelManager->NearestLevelEnergy(daughterExcitation);
          if (std::abs(daughterExcitation - nearestEnergy) < levelTolerance) {
            // Level half-life is in ns and the threshold is set to 1 micros
            // by default, user can set it via the UI command
            nearestLevelIndex = levelManager->NearestLevelIndex(daughterExcitation);
            if (levelManager->LifeTime(nearestLevelIndex)*ns >= halflifethreshold){
              // save the metastable nucleus 
              summedDecayTable->Insert(theChannel);
            } else {
              brs[theDecayMode] += theChannel->GetBR();
            }
          } else {
            brs[theDecayMode] += theChannel->GetBR();
          }
        } else {
          brs[theDecayMode] += theChannel->GetBR();
        }
      } // Combine decay channels (loop i)
	    
      brs[2] = brs[2]+brs[3]+brs[4]+brs[5];  // Combine beta+ and EC 
      brs[3] = brs[4] =brs[5] =  0.0;
      for (i= 0; i<nMode; i++){            // loop over decay modes
        if (brs[i] > 0.) {
          switch ( i ) {
          case 0:
            // Decay mode is isomeric transition
            theITChannel = new G4ITDecay(aParentNucleus, brs[0], 0.0, 0.0,
                                         photonEvaporation);

            summedDecayTable->Insert(theITChannel);
            break;

          case 1:
            // Decay mode is beta-
            theBetaMinusChannel = new G4BetaMinusDecay(aParentNucleus, brs[1],
                                                       0.*MeV, 0.*MeV,
                                                       noFloat, allowed);
            summedDecayTable->Insert(theBetaMinusChannel);
            break;

          case 2:
            // Decay mode is beta+ + EC.
            theBetaPlusChannel = new G4BetaPlusDecay(aParentNucleus, brs[2],    // DHW: April 2015
                                                     0.*MeV, 0.*MeV,
                                                     noFloat, allowed);
            summedDecayTable->Insert(theBetaPlusChannel);
            break;

          case 6:
            // Decay mode is alpha.
            theAlphaChannel = new G4AlphaDecay(aParentNucleus, brs[6], 0.*MeV,
                                               0.*MeV, noFloat);
            summedDecayTable->Insert(theAlphaChannel);
            break;

	  case 7:
            // Decay mode is proton.
            theProtonChannel = new G4ProtonDecay(aParentNucleus, brs[7], 0.*MeV,
                                                 0.*MeV, noFloat);
            summedDecayTable->Insert(theProtonChannel);
            break;
	  case 8:
            // Decay mode is neutron.
            theNeutronChannel = new G4NeutronDecay(aParentNucleus, brs[8], 0.*MeV,
                                                 0.*MeV, noFloat);
            summedDecayTable->Insert(theNeutronChannel);
            break;

          default:
            break;
          }
        }
      }
      // loop over all branches in summedDecayTable
      //
      for (i = 0; i < summedDecayTable->entries(); i++){
        theChannel = summedDecayTable->GetDecayChannel(i);
        theNuclearDecayChannel = static_cast<G4NuclearDecay*>(theChannel);
        theBR = theChannel->GetBR();
        theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus();

        // First check if the decay of the original nucleus is an IT channel,
        // if true create a new ground-state nucleus
        if (theNuclearDecayChannel->GetDecayMode() == IT && nGeneration == 1) {

          A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
          Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
          theDaughterNucleus=theIonTable->GetIon(Z,A,0.);
        }
        if (IsApplicable(*theDaughterNucleus) && theBR && 
            aParentNucleus != theDaughterNucleus) { 
          // need to make sure daughter has decay table
          parentDecayTable = GetDecayTable(theDaughterNucleus);

          if (parentDecayTable->entries() ) {
            A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
            Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
            E = ((const G4Ions*)(theDaughterNucleus))->GetExcitationEnergy();

            TaoPlus = theDaughterNucleus->GetPDGLifeTime();
            if (TaoPlus <= 0.)  TaoPlus = 1e-100;

            // first set the taos, one simply need to add to the parent ones
            taos.clear();
            taos = TP;   // load lifetimes of all previous generations 
            size_t k;
            //check that TaoPlus differs from other taos from at least 1.e5 relative difference
            //for (k = 0; k < TP.size(); k++){
            //if (std::abs((TaoPlus-TP[k])/TP[k])<1.e-5 ) TaoPlus=1.00001*TP[k];
            //}
            taos.push_back(TaoPlus);  // add daughter lifetime to list
            // now calculate the coefficiencies
            //
            // they are in two parts, first the less than n ones
            // Eq 4.24 of the TN
            Acoeffs.clear();
            long double ta1,ta2;
            ta2 = (long double)TaoPlus;
            for (k = 0; k < RP.size(); k++){
              ta1 = (long double)TP[k];    // loop over lifetimes of all previous generations
              if (ta1 == ta2) {
                theRate = 1.e100;
              } else {
                theRate = ta1/(ta1-ta2);
              }
              theRate = theRate * theBR * RP[k];
              Acoeffs.push_back(theRate);
            }

            // the second part: the n:n coefficiency
            // Eq 4.25 of the TN.  Note Yn+1 is zero apart from Y1 which is -1
            // as treated at line 1013 
            theRate = 0.;
            long double aRate, aRate1;
            aRate1 = 0.L;
            for (k = 0; k < RP.size(); k++){
              ta1 = (long double)TP[k];
              if (ta1 == ta2 ) {
                aRate = 1.e100;
              } else {
                aRate = ta2/(ta1-ta2);
              }
              aRate = aRate * (long double)(theBR * RP[k]);
              aRate1 += aRate;
            }
            theRate = -aRate1;
            Acoeffs.push_back(theRate); 	      
            SetDecayRate (Z,A,E,nGeneration,Acoeffs,taos);
            theDecayRateVector.push_back(theDecayRate);
            nEntry++;
          } // there are entries in the table
        } // nuclide is OK to decay
      } // end of loop (i) over decay table branches 
      //      delete summedDecayTable;

    } // Getting contents of decay rate vector (end loop on j)
    nS = nT;
    nT = nEntry;
    if (nS == nT) stable = true;
  } // while nuclide is not stable

  // end of while loop
  // the calculation completed here


  // fill the first part of the decay rate table
  // which is the name of the original particle (isotope)
  theDecayRateTable.SetIonName(theParentNucleus.GetParticleName()); 

  // now fill the decay table with the newly completed decay rate vector
  theDecayRateTable.SetItsRates(theDecayRateVector);

  // finally add the decayratetable to the tablevector
  theDecayRateTableVector.push_back(theDecayRateTable);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  SetSourceTimeProfile                                                      //
//     read in the source time profile function (histogram)                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void G4RadioactiveDecay::SetSourceTimeProfile(G4String filename)
{
  std::ifstream infile ( filename, std::ios::in );
  if (!infile) {
    G4ExceptionDescription ed;
    ed << " Could not open file " << filename << G4endl; 
    G4Exception("G4RadioactiveDecay::SetSourceTimeProfile()", "HAD_RDM_001",
                FatalException, ed);
  }

  G4double bin, flux;
  NSourceBin = -1;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
 
  while (infile >> bin >> flux ) {  /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 10000) {
      G4Exception("G4RadioactiveDecay::SetSourceTimeProfile()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
 
    NSourceBin++;
    if (NSourceBin > 99) {
      G4Exception("G4RadioactiveDecay::SetSourceTimeProfile()", "HAD_RDM_002",
                  FatalException, "Input source time file too big (>100 rows)");

    } else {
      SBin[NSourceBin] = bin * s;
      SProfile[NSourceBin] = flux;
    }
  }
  SetAnalogueMonteCarlo(0);
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout <<" Source Timeprofile Nbin = " << NSourceBin <<G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  SetDecayBiasProfile                                                       //
//    read in the decay bias scheme function (histogram)                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void G4RadioactiveDecay::SetDecayBias(G4String filename)
{
  
  std::ifstream infile(filename, std::ios::in);
  if (!infile) G4Exception("G4RadioactiveDecay::SetDecayBias()", "HAD_RDM_003",
                           FatalException, "Unable to open bias data file" );

  G4double bin, flux;
  G4int dWindows = 0;
  G4int i ;

  theRadioactivityTables.clear();

  NDecayBin = -1;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
 
  while (infile >> bin >> flux ) {  /* Loop checking, 01.09.2015, D.Wright */
    NDecayBin++;
    loop++;
    if (loop > 10000) {
      G4Exception("G4RadioactiveDecay::SetDecayBias()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
 
    if (NDecayBin > 99) {
      G4Exception("G4RadioactiveDecay::SetDecayBias()", "HAD_RDM_004",
                  FatalException, "Input bias file too big (>100 rows)" );
    } else {
      DBin[NDecayBin] = bin * s;
      DProfile[NDecayBin] = flux;
      if (flux > 0.) {
        decayWindows[NDecayBin] = dWindows;
        dWindows++;
        G4RadioactivityTable *rTable = new G4RadioactivityTable() ;
        theRadioactivityTables.push_back(rTable);
      }
    }
  }
  for ( i = 1; i<= NDecayBin; i++) DProfile[i] += DProfile[i-1];
  for ( i = 0; i<= NDecayBin; i++) DProfile[i] /= DProfile[NDecayBin];
  // converted to accumulated probabilities

  SetAnalogueMonteCarlo(0);
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout <<" Decay Bias Profile  Nbin = " << NDecayBin <<G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DecayIt                                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4VParticleChange*
G4RadioactiveDecay::DecayIt(const G4Track& theTrack, const G4Step&)
{
  // Initialize G4ParticleChange object, get particle details and decay table

  fParticleChangeForRadDecay.Initialize(theTrack);
  fParticleChangeForRadDecay.ProposeWeight(theTrack.GetWeight());
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();

  // First check whether RDM applies to the current logical volume
  if (!isAllVolumesMode) {
    if (!std::binary_search(ValidVolumes.begin(), ValidVolumes.end(),
                     theTrack.GetVolume()->GetLogicalVolume()->GetName())) {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
        G4cout <<"G4RadioactiveDecay::DecayIt : "
               << theTrack.GetVolume()->GetLogicalVolume()->GetName()
               << " is not selected for the RDM"<< G4endl;
        G4cout << " There are " << ValidVolumes.size() << " volumes" << G4endl;
        G4cout << " The Valid volumes are " << G4endl;
        for (size_t i = 0; i< ValidVolumes.size(); i++)
                                  G4cout << ValidVolumes[i] << G4endl;
      }
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

      // Kill the parent particle.
      fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }
  }

  // Now check if particle is valid for RDM
  if (!(IsApplicable(*theParticleDef) ) ) { 
    // Particle is not an ion or is outside the nucleuslimits for decay

#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4RadioactiveDecay::DecayIt : "
             << theParticleDef->GetParticleName() 
             << " is not a valid nucleus for the RDM"<< G4endl;
    }
#endif
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;
  }
  G4DecayTable* theDecayTable = GetDecayTable(theParticleDef);

  if (theDecayTable == 0 || theDecayTable->entries() == 0) {
    // No data in the decay table.  Set particle change parameters
    // to indicate this.
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 0) {
      G4cerr <<"G4RadioactiveDecay::DecayIt : decay table not defined  for ";
      G4cerr <<theParticleDef->GetParticleName() <<G4endl;
    }
#endif
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle.
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;

  } else { 
    // Data found.  Try to decay nucleus
    G4double energyDeposit = 0.0;
    G4double finalGlobalTime = theTrack.GetGlobalTime();
    G4double finalLocalTime = theTrack.GetLocalTime();
    G4int index;
    G4ThreeVector currentPosition;
    currentPosition = theTrack.GetPosition();

    // Check whether use Analogue or VR implementation
    if (AnalogueMC) {
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 0)
        G4cout <<"DecayIt:  Analogue MC version " << G4endl;
# endif

      G4DecayProducts* products = DoDecay(*theParticleDef);

      // Check if the product is the same as input and kill the track if
      // necessary to prevent infinite loop (11/05/10, F.Lei)
      if ( products->entries() == 1) {
        fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
        fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill);
        fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
        ClearNumberOfInteractionLengthLeft();
        return &fParticleChangeForRadDecay;
      }

      // Get parent particle information and boost the decay products to the
      // laboratory frame based on this information.

      //The Parent Energy used for the boost should be the total energy of
      // the nucleus of the parent ion without the energy of the shell electrons
      // (correction for bug 1359 by L. Desorgher)
      G4double ParentEnergy = theParticle->GetKineticEnergy()
                            + theParticle->GetParticleDefinition()->GetPDGMass();
      G4ThreeVector ParentDirection(theParticle->GetMomentumDirection());

      if (theTrack.GetTrackStatus() == fStopButAlive) {
        //this condition seems to be always True, further investigation is needed (L.Desorgher)

        // The particle is decayed at rest.
        // since the time is still for rest particle in G4 we need to add the
        // additional time lapsed between the particle come to rest and the
        // actual decay.  This time is simply sampled with the mean-life of
        // the particle.  But we need to protect the case PDGTime < 0.
        // (F.Lei 11/05/10)
        G4double temptime = -std::log( G4UniformRand())
                            *theParticleDef->GetPDGLifeTime();
        if (temptime < 0.) temptime = 0.; 
        finalGlobalTime += temptime;
        finalLocalTime += temptime;
        energyDeposit += theParticle->GetKineticEnergy();
      }
      products->Boost(ParentEnergy, ParentDirection);

      // Add products in theParticleChangeForRadDecay.
      G4int numberOfSecondaries = products->entries();
      fParticleChangeForRadDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout <<"G4RadioactiveDecay::DecayIt : Decay vertex :";
        G4cout <<" Time: " <<finalGlobalTime/ns <<"[ns]";
        G4cout <<" X:" <<(theTrack.GetPosition()).x() /cm <<"[cm]";
        G4cout <<" Y:" <<(theTrack.GetPosition()).y() /cm <<"[cm]";
        G4cout <<" Z:" <<(theTrack.GetPosition()).z() /cm <<"[cm]";
        G4cout << G4endl;
        G4cout <<"G4Decay::DecayIt  : decay products in Lab. Frame" <<G4endl;
        products->DumpInfo();
        products->IsChecked();
      }
#endif
      for (index=0; index < numberOfSecondaries; index++) {
        G4Track* secondary = new G4Track(products->PopProducts(),
                                         finalGlobalTime, currentPosition);
        secondary->SetGoodForTrackingFlag();
        secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
        fParticleChangeForRadDecay.AddSecondary(secondary);
      }
      delete products;
      // end of analogue MC algorithm

    } else {
      // Variance Reduction Method
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
        G4cout << "DecayIt: Variance Reduction version " << G4endl;
#endif
      if (!IsRateTableReady(*theParticleDef)) {
        // if the decayrates are not ready, calculate them and 
        // add to the rate table vector 
        AddDecayRateTable(*theParticleDef);
      }
      //retrieve the rates 
      GetDecayRateTable(*theParticleDef);

      // declare some of the variables required in the implementation
      G4ParticleDefinition* parentNucleus;
      G4IonTable* theIonTable;
      G4int PZ;
      G4int PA;
      G4double PE;
      G4String keyName;
      std::vector<G4double> PT;
      std::vector<G4double> PR;
      G4double taotime;
      long double decayRate;

      size_t i;
      size_t j;
      G4int numberOfSecondaries;
      G4int totalNumberOfSecondaries = 0;
      G4double currentTime = 0.;
      G4int ndecaych;
      G4DynamicParticle* asecondaryparticle;
      std::vector<G4DynamicParticle*> secondaryparticles;
      std::vector<G4double> pw;
      std::vector<G4double> ptime;
      pw.clear();
      ptime.clear();

      //now apply the nucleus splitting
      for (G4int n = 0; n < NSplit; n++) {
        // Get the decay time following the decay probability function 
        // suppllied by user  
        G4double theDecayTime = GetDecayTime();
        G4int nbin = GetDecayTimeBin(theDecayTime);

        // calculate the first part of the weight function  
        G4double weight1 = 1.; 
        if (nbin == 1) {
          weight1 = 1./DProfile[nbin-1] 
                    *(DBin[nbin]-DBin[nbin-1])/NSplit;
        } else if (nbin > 1) {
          weight1 = 1./(DProfile[nbin]-DProfile[nbin-2])
                    *(DBin[nbin]-DBin[nbin-1])/NSplit;
        }

        // it should be calculated in seconds
        weight1 /= s ;
	    
        // loop over all the possible secondaries of the nucleus
        // the first one is itself.
        for (i = 0; i < theDecayRateVector.size(); i++) {
          PZ = theDecayRateVector[i].GetZ();
          PA = theDecayRateVector[i].GetA();
          PE = theDecayRateVector[i].GetE();
          PT = theDecayRateVector[i].GetTaos();
          PR = theDecayRateVector[i].GetDecayRateC();

          // Calculate the decay rate of the isotope
          // decayRate is the radioactivity of isotope (PZ,PA,PE) at the 
          // time 'theDecayTime'
          // it will be used to calculate the statistical weight of the 
          // decay products of this isotope

//          G4cout <<"PA= "<< PA << " PZ= " << PZ << " PE= "<< PE  <<G4endl;
          decayRate = 0.L;
          for (j = 0; j < PT.size(); j++) {
//            G4cout <<  " RDM::DecayIt: tau input to Convolve: " <<  PT[j] << G4endl; 
            taotime = ConvolveSourceTimeProfile(theDecayTime,PT[j]);
//            taotime = GetTaoTime(theDecayTime,PT[j]);
            decayRate -= PR[j] * (long double)taotime;
            // Eq.4.23 of of the TN
            // note the negative here is required as the rate in the
            // equation is defined to be negative, 
            // i.e. decay away, but we need positive value here.

            // G4cout << j << "\t"<< PT[j]/s <<"\t"<<PR[j]<< "\t"
            //        << decayRate << G4endl;		
          }

          // add the isotope to the radioactivity tables
          //  G4cout <<theDecayTime/s <<"\t"<<nbin<<G4endl;
          //  G4cout << theTrack.GetWeight() <<"\t"<<weight1<<"\t"<<decayRate<< G4endl;
          theRadioactivityTables[decayWindows[nbin-1]]->AddIsotope(PZ,PA,PE,weight1*decayRate,theTrack.GetWeight());

          // Now calculate the statistical weight
          // One needs to fold the source bias function with the decaytime
          // also need to include the track weight! (F.Lei, 28/10/10)
          G4double weight = weight1*decayRate*theTrack.GetWeight();

          // decay the isotope 
          theIonTable = (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
          parentNucleus = theIonTable->GetIon(PZ,PA,PE);

          // Create a temprary products buffer.
          // Its contents to be transfered to the products at the end of the loop
          G4DecayProducts* tempprods = 0;

          // Decide whether to apply branching ratio bias or not	     
          if (BRBias) {
            G4DecayTable* decayTable = GetDecayTable(parentNucleus);

            ndecaych = G4int(decayTable->entries()*G4UniformRand());
            G4VDecayChannel* theDecayChannel = decayTable->GetDecayChannel(ndecaych);
            if (theDecayChannel == 0) {
              // Decay channel not found.
#ifdef G4VERBOSE
              if (GetVerboseLevel()>0) {
                G4cerr << " G4RadioactiveDecay::DoIt : cannot determine decay channel ";
                G4cerr << " for this nucleus; decay as if no biasing active ";
                G4cerr << G4endl;
                decayTable ->DumpInfo();
              }
#endif
              tempprods = DoDecay(*parentNucleus);  // DHW 6 Dec 2010 - do decay as if no biasing
                                                    //           to avoid deref of temppprods = 0
            } else {
              // A decay channel has been identified, so execute the DecayIt.
              G4double tempmass = parentNucleus->GetPDGMass();
              tempprods = theDecayChannel->DecayIt(tempmass);
              weight *= (theDecayChannel->GetBR())*(decayTable->entries());
            }
          } else {
            tempprods = DoDecay(*parentNucleus);
          }


          // save the secondaries for buffers
          numberOfSecondaries = tempprods->entries();
          currentTime = finalGlobalTime + theDecayTime;
          for (index = 0; index < numberOfSecondaries; index++) {
            asecondaryparticle = tempprods->PopProducts();
            if (asecondaryparticle->GetDefinition()->GetBaryonNumber() < 5) {
              pw.push_back(weight);
              ptime.push_back(currentTime);
              secondaryparticles.push_back(asecondaryparticle);
            }
            //Generate gammas and XRays from  excited nucleus, added by L.Desorgher
            else if (((const G4Ions*)(asecondaryparticle->GetDefinition()))->GetExcitationEnergy()>0. && weight>0.){//Compute the gamma
              G4ParticleDefinition* apartDef =asecondaryparticle->GetDefinition();
              AddDeexcitationSpectrumForBiasMode(apartDef,weight,currentTime,pw,ptime,secondaryparticles);
            }
          }
          delete tempprods;

        } // end of i loop
      } // end of n loop 

      // now deal with the secondaries in the two stl containers
      // and submmit them back to the tracking manager
      totalNumberOfSecondaries = pw.size();
      fParticleChangeForRadDecay.SetNumberOfSecondaries(totalNumberOfSecondaries);
      for (index=0; index < totalNumberOfSecondaries; index++) { 
        G4Track* secondary = new G4Track(secondaryparticles[index],
                                         ptime[index], currentPosition);
        secondary->SetGoodForTrackingFlag(); 	   
        secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
        secondary->SetWeight(pw[index]); 	   
        fParticleChangeForRadDecay.AddSecondary(secondary); 
      }
      // make sure the original track is set to stop and its kinematic energy collected
      // 
      //theTrack.SetTrackStatus(fStopButAlive);
      //energyDeposit += theParticle->GetKineticEnergy();

    } // End of Variance Reduction 

    // Kill the parent particle
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(energyDeposit); 
    fParticleChangeForRadDecay.ProposeLocalTime(finalLocalTime);
    // Reset NumberOfInteractionLengthLeft.
    ClearNumberOfInteractionLengthLeft();

    return &fParticleChangeForRadDecay ;
  }
} 


G4DecayProducts*
G4RadioactiveDecay::DoDecay(const G4ParticleDefinition& theParticleDef)
{
  G4DecayProducts* products = 0;
  G4DecayTable* theDecayTable = GetDecayTable(&theParticleDef);

  // Choose a decay channel.
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "Select a channel..." << G4endl;
#endif

  // G4DecayTable::SelectADecayChannel checks to see if sum of daughter masses
  // exceeds parent mass. Pass it the parent mass + maximum Q value to account
  // for difference in mass defect.
  G4double parentPlusQ = theParticleDef.GetPDGMass() + 30.*MeV;
  G4VDecayChannel* theDecayChannel = theDecayTable->SelectADecayChannel(parentPlusQ);

  if (theDecayChannel == 0) {
    // Decay channel not found.
    G4ExceptionDescription ed;
    ed << " Cannot determine decay channel for " << theParticleDef.GetParticleName() << G4endl;
    G4Exception("G4RadioactiveDecay::DoDecay", "HAD_RDM_013",
                FatalException, ed);
  } else {
    // A decay channel has been identified, so execute the DecayIt.
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cerr << "G4RadioactiveDecay::DoIt : selected decay channel  addr:";
      G4cerr << theDecayChannel << G4endl;
    }
#endif
    products = theDecayChannel->DecayIt(theParticleDef.GetPDGMass() );

    // Apply directional bias if requested by user
    CollimateDecay(products);
  }

  return products;
}


// Apply directional bias for "visible" daughters (e+-, gamma, n, p, alpha)

void G4RadioactiveDecay::CollimateDecay(G4DecayProducts* products) {
  if (origin == forceDecayDirection) return;	// No collimation requested
  if (180.*deg == forceDecayHalfAngle) return;
  if (0 == products || 0 == products->entries()) return;

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "Begin of CollimateDecay..." << G4endl;
#endif

  // Particles suitable for directional biasing (for if-blocks below)
  static const G4ParticleDefinition* electron = G4Electron::Definition();
  static const G4ParticleDefinition* positron = G4Positron::Definition();
  static const G4ParticleDefinition* neutron  = G4Neutron::Definition();
  static const G4ParticleDefinition* gamma    = G4Gamma::Definition();
  static const G4ParticleDefinition* alpha    = G4Alpha::Definition();
  static const G4ParticleDefinition* proton   = G4Proton::Definition();

  G4ThreeVector newDirection;		// Re-use to avoid memory churn
  for (G4int i=0; i<products->entries(); i++) {
    G4DynamicParticle* daughter = (*products)[i];
    const G4ParticleDefinition* daughterType =
                                  daughter->GetParticleDefinition();
    if (daughterType == electron || daughterType == positron ||
	daughterType == neutron || daughterType == gamma ||
	daughterType == alpha || daughterType == proton) CollimateDecayProduct(daughter);
  }
}

void G4RadioactiveDecay::CollimateDecayProduct(G4DynamicParticle* daughter) {
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "CollimateDecayProduct for daughter "
	   << daughter->GetParticleDefinition()->GetParticleName() << G4endl;
  }
#endif

  G4ThreeVector collimate = ChooseCollimationDirection();
  if (origin != collimate) daughter->SetMomentumDirection(collimate);
}


// Choose random direction within collimation cone

G4ThreeVector G4RadioactiveDecay::ChooseCollimationDirection() const {
  if (origin == forceDecayDirection) return origin;	// Don't do collimation
  if (forceDecayHalfAngle == 180.*deg) return origin;

  G4ThreeVector dir = forceDecayDirection;

  // Return direction offset by random throw
  if (forceDecayHalfAngle > 0.) {
    // Generate uniform direction around central axis
    G4double phi = 2.*pi*G4UniformRand();
    G4double cosMin = std::cos(forceDecayHalfAngle);
    G4double cosTheta = (1.-cosMin)*G4UniformRand() + cosMin;	// [cosMin,1.)
    
    dir.setPhi(dir.phi()+phi);
    dir.setTheta(dir.theta()+std::acos(cosTheta));
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    G4cout << " ChooseCollimationDirection returns " << dir << G4endl;
#endif

  return dir;
}

//Add gamma,Xray,conversion,and auger electrons for bias mode
void G4RadioactiveDecay::AddDeexcitationSpectrumForBiasMode(G4ParticleDefinition* apartDef,
                                            G4double weight,G4double currentTime,
                                            std::vector<double>& weights_v,
                                            std::vector<double>& times_v,
                                            std::vector<G4DynamicParticle*>& secondaries_v)
{
  G4double elevel=((const G4Ions*)(apartDef))->GetExcitationEnergy();
  G4double life_time=apartDef->GetPDGLifeTime();
  while (life_time <halflifethreshold && elevel>0.) {
    G4ITDecay* anITChannel = new G4ITDecay(apartDef, 100., elevel,elevel,
                                           photonEvaporation);
    G4DecayProducts* pevap_products = anITChannel->DecayIt(0.);
    G4int nb_pevapSecondaries = pevap_products->entries();
    for (G4int ind = 0; ind < nb_pevapSecondaries; ind++) {
		G4DynamicParticle* a_pevap_secondary= pevap_products->PopProducts();
		//Gammas,electrons, alphas coming from excited state
		if (a_pevap_secondary->GetDefinition()->GetBaryonNumber() < 5) {
		  weights_v.push_back(weight);
		  times_v.push_back(currentTime);
		  secondaries_v.push_back(a_pevap_secondary);
		}
		//New excited or ground state
		else {
		  apartDef =a_pevap_secondary->GetDefinition();
		  elevel=((const G4Ions*)(apartDef))->GetExcitationEnergy();
		  life_time=apartDef->GetPDGLifeTime();
		}
    }
    delete anITChannel;
  }
}


