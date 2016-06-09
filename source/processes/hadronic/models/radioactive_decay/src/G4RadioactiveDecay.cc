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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecaymessenger.hh"

#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4ITDecayChannel.hh"
#include "G4BetaMinusDecayChannel.hh"
#include "G4BetaPlusDecayChannel.hh"
#include "G4KshellECDecayChannel.hh"
#include "G4LshellECDecayChannel.hh"
#include "G4MshellECDecayChannel.hh"
#include "G4AlphaDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4RIsotopeTable.hh"
#include "G4BetaFermiFunction.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"

#include "G4HadTmpUtil.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace CLHEP;

const G4double   G4RadioactiveDecay::levelTolerance =2.0*keV;

///////////////////////////////////////////////////////////////////////////////
//
//
// Constructor
//
G4RadioactiveDecay::G4RadioactiveDecay (const G4String& processName)
 :G4VRestDiscreteProcess(processName, fDecay), HighestBinValue(10.0),
  LowestBinValue(1.0e-3), TotBin(200), verboseLevel(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout <<"G4RadioactiveDecay constructor    Name: ";
    G4cout <<processName << G4endl;   }
#endif

  theRadioactiveDecaymessenger = new G4RadioactiveDecaymessenger(this);
  theIsotopeTable              = new G4RIsotopeTable();
  aPhysicsTable                = 0;
  pParticleChange              = &fParticleChangeForRadDecay;

  //
  // Now register the Isotopetable with G4IonTable.
  //
  G4IonTable *theIonTable =
    (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4VIsotopeTable *aVirtualTable = theIsotopeTable;
  theIonTable->RegisterIsotopeTable(aVirtualTable);
  //
  //
  // Reset the contents of the list of nuclei for which decay scheme data
  // have been loaded.
  //
  LoadedNuclei.clear();
  //
  //
  // Apply default values.
  //
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
  halflifethreshold = 0.*second;
  //
  // RDM applies to xall logical volumes as default
  SelectAllVolumes();
}
////////////////////////////////////////////////////////////////////////////////
//
//
// Destructor
//
G4RadioactiveDecay::~G4RadioactiveDecay()
{
  if (aPhysicsTable != 0)
    {
      aPhysicsTable->clearAndDestroy();
      delete aPhysicsTable;
    }
  //    delete theIsotopeTable;
    delete theRadioactiveDecaymessenger;
}

///////////////////////////////////////////////////////////////////////////////
//
//
// IsApplicable
//
G4bool G4RadioactiveDecay::IsApplicable( const G4ParticleDefinition &
					 aParticle)
{
  //
  //
  // All particles, other than G4Ions, are rejected by default.
  //
  if (((const G4Ions*)(&aParticle))->GetExcitationEnergy() > 0.) {return true;}
  if (aParticle.GetParticleName() == "GenericIon")      {return true;}
  else if (!(aParticle.GetParticleType() == "nucleus") || aParticle.GetPDGLifeTime() < 0. ) {return false;}
  //
  //
  // Determine whether the nuclide falls into the correct A and Z range.
  //
  G4int A = ((const G4Ions*) (&aParticle))->GetAtomicMass();
  G4int Z = ((const G4Ions*) (&aParticle))->GetAtomicNumber();
  if( A> theNucleusLimits.GetAMax() || A< theNucleusLimits.GetAMin())
    {return false;}
  else if( Z> theNucleusLimits.GetZMax() || Z< theNucleusLimits.GetZMin())
    {return false;}
  return true;
}
///////////////////////////////////////////////////////////////////////////////
//
//
// IsLoaded
//
G4bool G4RadioactiveDecay::IsLoaded(const G4ParticleDefinition &
				    aParticle)
{
  //
  //
  // Check whether the radioactive decay data on the ion have already been
  // loaded.
  //
  return std::binary_search(LoadedNuclei.begin(),
			    LoadedNuclei.end(),
			    aParticle.GetParticleName());
}
////////////////////////////////////////////////////////////////////////////////
//
//
// SelectAVolume
//
void G4RadioactiveDecay::SelectAVolume(const G4String aVolume)
{

  G4LogicalVolumeStore *theLogicalVolumes;
  G4LogicalVolume *volume;
  theLogicalVolumes=G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume=(*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      ValidVolumes.push_back(aVolume);
      std::sort(ValidVolumes.begin(), ValidVolumes.end());
      // sort need for performing binary_search
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	G4cout << " RDM Applies to : " << aVolume << G4endl; 
#endif
    }else if(i ==  theLogicalVolumes->size())
      {
	G4cerr << "SelectAVolume: "<<aVolume << " is not a valid logical volume name"<< G4endl; 
      }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
// DeSelectAVolume
//
void G4RadioactiveDecay::DeselectAVolume(const G4String aVolume)
{
  G4LogicalVolumeStore *theLogicalVolumes;
  G4LogicalVolume *volume;
  theLogicalVolumes=G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume=(*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      std::vector<G4String>::iterator location;
      location = std::find(ValidVolumes.begin(),ValidVolumes.end(),aVolume);
      if (location != ValidVolumes.end()) {
	ValidVolumes.erase(location);
	std::sort(ValidVolumes.begin(), ValidVolumes.end());
      }else{
	G4cerr << " DeselectVolume:" << aVolume << " is not in the list"<< G4endl; 
      }	  
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	G4cout << " DeselectVolume: " << aVolume << " is removed from list"<<G4endl; 
#endif
    }else if(i ==  theLogicalVolumes->size())
      {
	G4cerr << " DeselectVolume:" << aVolume << "is not a valid logical volume name"<< G4endl; 
      }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
// SelectAllVolumes
//
void G4RadioactiveDecay::SelectAllVolumes() 
{
  G4LogicalVolumeStore *theLogicalVolumes;
  G4LogicalVolume *volume;
  theLogicalVolumes=G4LogicalVolumeStore::GetInstance();
  ValidVolumes.clear();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout << " RDM Applies to all Volumes"  << G4endl;
#endif
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume=(*theLogicalVolumes)[i];
    ValidVolumes.push_back(volume->GetName());    
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout << "         RDM Applies to Volume "  << volume->GetName() << G4endl;
#endif
  }
  std::sort(ValidVolumes.begin(), ValidVolumes.end());
  // sort needed in order to allow binary_search
}
////////////////////////////////////////////////////////////////////////////////
//
//
// DeSelectAllVolumes
//
void G4RadioactiveDecay::DeselectAllVolumes() 
{
  ValidVolumes.clear();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout << " RDM removed from all volumes" << G4endl; 
#endif

}

///////////////////////////////////////////////////////////////////////////////
//
//
// IsRateTableReady
//
G4bool G4RadioactiveDecay::IsRateTableReady(const G4ParticleDefinition &
					    aParticle)
{
  //
  //
  // Check whether the radioactive decay rates table on the ion has already
  // been calculated.
  //
  G4String aParticleName = aParticle.GetParticleName();
  for (size_t i = 0; i < theDecayRateTableVector.size(); i++)
    {
      if (theDecayRateTableVector[i].GetIonName() == aParticleName)
	return true ;
    }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
//
//
// GetDecayRateTable
//
// retrieve the decayratetable for the specified aParticle
//
void G4RadioactiveDecay::GetDecayRateTable(const G4ParticleDefinition &
					   aParticle)
{

  G4String aParticleName = aParticle.GetParticleName();

  for (size_t i = 0; i < theDecayRateTableVector.size(); i++)
    {
      if (theDecayRateTableVector[i].GetIonName() == aParticleName)
	{
	  theDecayRateVector = theDecayRateTableVector[i].GetItsRates() ;
	}
    }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    {
      G4cout <<"The DecayRate Table for "
	     << aParticleName << " is selected." <<  G4endl;
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
// GetTaoTime
//
// to perform the convolution of the source time profile function with the 
// decay constants in the decay chain. 
//
G4double G4RadioactiveDecay::GetTaoTime(G4double t, G4double tao)
{
  G4double taotime =0.;
  G4int nbin;
  if ( t > SBin[NSourceBin]) {
    nbin  = NSourceBin;}
  else {
    nbin = 0;
    while (t > SBin[nbin]) nbin++;
    nbin--;}
  if (nbin > 0) { 
    for (G4int i = 0; i < nbin; i++) 
      {
	taotime += SProfile[i] * (std::exp(-(t-SBin[i+1])/tao)-std::exp(-(t-SBin[i])/tao));
      }
  }
  taotime +=  SProfile[nbin] * (1-std::exp(-(t-SBin[nbin])/tao));
  if (taotime < 0.)  {
    G4cout <<" Tao time: " <<taotime << " reset to zero!"<<G4endl;
    taotime = 0.;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Tao time: " <<taotime <<G4endl;}
#endif
  return taotime ;
}

////////////////////////////////////////////////////////////////////////////////
//
//
// GetDecayTime
//
// Randomly select a decay time for the decay process, following the supplied
// decay time bias scheme.
//
G4double G4RadioactiveDecay::GetDecayTime()
{
  G4double decaytime = 0.;
  G4double rand = G4UniformRand();
  G4int i = 0;
  while ( DProfile[i] < rand) i++;
  rand = G4UniformRand();
  decaytime = DBin[i] + rand*(DBin[i+1]-DBin[i]);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Decay time: " <<decaytime/s <<"[s]" <<G4endl;}
#endif
  return  decaytime;	    
}

////////////////////////////////////////////////////////////////////////////////
//
//
// GetDecayTimeBin
//
//
//
G4int G4RadioactiveDecay::GetDecayTimeBin(const G4double aDecayTime)
{
  G4int i = 0;
  while ( aDecayTime > DBin[i] ) i++;
  return  i;
}
////////////////////////////////////////////////////////////////////////////////
//
//
// GetMeanLifeTime
//
// this is required by the base class 
//
G4double G4RadioactiveDecay::GetMeanLifeTime(const G4Track& theTrack,
					     G4ForceCondition* )
{
  // For varience reduction implementation the time is set to 0 so as to 
  // force the particle to decay immediately.
  // in analogueMC mode it return the particles meanlife.
  // 
  G4double meanlife = 0.;
  if (AnalogueMC) {
    const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
    G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();
    G4double theLife = theParticleDef->GetPDGLifeTime();

#ifdef G4VERBOSE
    if (GetVerboseLevel()>2)
      {
	G4cout <<"G4RadioactiveDecay::GetMeanLifeTime() " <<G4endl;
	G4cout <<"KineticEnergy:" <<theParticle->GetKineticEnergy()/GeV <<"[GeV]";
	G4cout <<"Mass:" <<theParticle->GetMass()/GeV <<"[GeV]"; 
	G4cout <<"Life time: " <<theLife/ns <<"[ns]" << G4endl;
      }
#endif
    if (theParticleDef->GetPDGStable()) {meanlife = DBL_MAX;}
    else if (theLife < 0.0) {meanlife = DBL_MAX;}
    else {meanlife = theLife;}
    // set the meanlife to zero for excited istopes which is not in the RDM database
    if (((const G4Ions*)(theParticleDef))->GetExcitationEnergy() > 0. && meanlife == DBL_MAX) {meanlife = 0.;}
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"mean life time: " <<meanlife/s <<"[s]" <<G4endl;}
#endif

  return  meanlife;
}
////////////////////////////////////////////////////////////////////////////////
//
//
// GetMeanFreePath
//
// it is of similar functions to GetMeanFreeTime 
//
G4double G4RadioactiveDecay::GetMeanFreePath (const G4Track& aTrack,
					      G4double, G4ForceCondition*)
{
  // constants
  G4bool isOutRange ;

  // get particle
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

  // returns the mean free path in GEANT4 internal units
  G4double pathlength;
  G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
  G4double aCtau = c_light * aParticleDef->GetPDGLifeTime();
  G4double aMass = aParticle->GetMass();

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
    G4cout << "G4RadioactiveDecay::GetMeanFreePath() "<< G4endl;
    G4cout << "KineticEnergy:" << aParticle->GetKineticEnergy()/GeV <<"[GeV]";
    G4cout << "Mass:" << aMass/GeV <<"[GeV]";
    G4cout << "c*Tau:" << aCtau/m <<"[m]" <<G4endl;
  }
#endif

  // check if the particle is stable?
  if (aParticleDef->GetPDGStable()) {
    pathlength = DBL_MAX;

  } else if (aCtau < 0.0) {
    pathlength =  DBL_MAX;

    //check if the particle has very short life time ?
  } else if (aCtau < DBL_MIN) {
    pathlength =  DBL_MIN;

    //check if zero mass
  } else if (aMass <  DBL_MIN)  {
    pathlength =  DBL_MAX;
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cerr << " Zero Mass particle " << G4endl;
    }
#endif
  } else {
    //calculate the mean free path
    // by using normalized kinetic energy (= Ekin/mass)
    G4double   rKineticEnergy = aParticle->GetKineticEnergy()/aMass;
    if ( rKineticEnergy > HighestBinValue) {
      // beta >> 1
      pathlength = ( rKineticEnergy + 1.0)* aCtau;
    } else if ( rKineticEnergy > LowestBinValue) {
      // check if aPhysicsTable exists
      if (aPhysicsTable == 0) BuildPhysicsTable(*aParticleDef);
      // beta is in the range valid for PhysicsTable
      pathlength = aCtau *
	((*aPhysicsTable)(0))-> GetValue(rKineticEnergy,isOutRange);
    } else if ( rKineticEnergy < DBL_MIN ) {
      // too slow particle
#ifdef G4VERBOSE
      if (GetVerboseLevel()>2) {
	G4cout << "G4Decay::GetMeanFreePath()   !!particle stops!!";
	G4cout << aParticleDef->GetParticleName() << G4endl;
	G4cout << "KineticEnergy:" << aParticle->GetKineticEnergy()/GeV <<"[GeV]";
      }
#endif
      pathlength = DBL_MIN;
    } else {
      // beta << 1
      pathlength = (aParticle->GetTotalMomentum())/aMass*aCtau ;
    }
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "mean free path: "<< pathlength/m << "[m]" << G4endl;
  }
#endif
  return  pathlength;
}
////////////////////////////////////////////////////////////////////////////////
//
//
//
//
void G4RadioactiveDecay::BuildPhysicsTable(const G4ParticleDefinition&)
{
  // if aPhysicsTableis has already been created, do nothing
  if (aPhysicsTable != 0) return;

  // create  aPhysicsTable
  if (GetVerboseLevel()>1) G4cerr <<" G4Decay::BuildPhysicsTable() "<< G4endl;
  aPhysicsTable = new G4PhysicsTable(1);

  //create physics vector
  G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
						       LowestBinValue,
						       HighestBinValue,
						       TotBin);

  G4double beta, gammainv;
  // fill physics Vector
  G4int i;
  for ( i = 0 ; i < TotBin ; i++ ) {
    gammainv = 1.0/(aVector->GetLowEdgeEnergy(i) + 1.0);
    beta  = std::sqrt((1.0 - gammainv)*(1.0 +gammainv));
    aVector->PutValue(i, beta/gammainv);
  }
  aPhysicsTable->insert(aVector);
}
///////////////////////////////////////////////////////////////////////////////
//
// LoadDecayTable
// 
//     To load the decay scheme from the Radioactivity database for 
//     theParentNucleus.
//
G4DecayTable *G4RadioactiveDecay::LoadDecayTable (G4ParticleDefinition &theParentNucleus)
{
  //
  //
  // Create and initialise variables used in the method.
  //
  G4DecayTable *theDecayTable = new G4DecayTable();
  //
  //
  // Determine the filename of the file containing radioactive decay data.  Open
  // it.
  //
  G4int A    = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z    = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();
  G4double E = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();

  if ( !getenv("G4RADIOACTIVEDATA") ) {
    G4cout << "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files." << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
			      "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files.");
  }
  G4String dirName = getenv("G4RADIOACTIVEDATA");
  LoadedNuclei.push_back(theParentNucleus.GetParticleName());
  std::sort( LoadedNuclei.begin(), LoadedNuclei.end() );
  // sort needed to allow binary_search

  std::ostringstream os;
  os <<dirName <<"/z" <<Z <<".a" <<A <<'\0';
  G4String file = os.str();


  std::ifstream DecaySchemeFile(file);
  G4bool found(false);
  if (DecaySchemeFile) { 
    //
    //
    // Initialise variables used for reading in radioactive decay data.
    //
    G4int nMode = 7;
    G4bool modeFirstRecord[7];
    G4double modeTotalBR[7] = {0.0};
    G4double modeSumBR[7];
    G4int i;
    for (i=0; i<nMode; i++) {
      modeFirstRecord[i] = true;
      modeSumBR[i] = 0.0;
    }

    G4bool complete(false);
    char inputChars[80]={' '};
    G4String inputLine;
    G4String recordType("");
    G4RadioactiveDecayMode theDecayMode;
    G4double a(0.0);
    G4double b(0.0);
    G4double c(0.0);
    G4double n(1.0);
    G4double e0;
    //
    //
    // Go through each record in the data file until you identify the decay
    // data relating to the nuclide of concern.
    //
    //  while (!complete && -DecaySchemeFile.getline(inputChars, 80).eof() != EOF)
    while (!complete && !DecaySchemeFile.getline(inputChars, 80).eof()) {
      inputLine = inputChars;
      //    G4String::stripType stripend(1);
      //    inputLine = inputLine.strip(stripend);
      inputLine = inputLine.strip(1);
      if (inputChars[0] != '#' && inputLine.length() != 0) {
	std::istringstream tmpStream(inputLine);
	if (inputChars[0] == 'P') {
	  //
	  //
	  // Nucleus is a parent type.  Check the excitation level to see if it matches
	  // that of theParentNucleus
	  //
	  tmpStream >>recordType >>a >>b;
	  if (found) {complete = true;}
	  else {found = (std::abs(a*keV - E)<levelTolerance);}
	} else if (found) {
	  //
	  //
	  // The right part of the radioactive decay data file has been found.  Search
	  // through it to determine the mode of decay of the subsequent records.
	  //
	  if (inputChars[0] == 'W') {
#ifdef G4VERBOSE
	    if (GetVerboseLevel()>0) {
	      // a comment line identified and print out the message
	      //
	      G4cout << " Warning in G4RadioactiveDecay::LoadDecayTable " << G4endl;
	      G4cout << "   In data file " << file << G4endl;
	      G4cout << "   " << inputLine << G4endl;
	    }
#endif
	  }	
	  else {
	    tmpStream >>theDecayMode >>a >>b >>c;
	    a/=1000.;
	    c/=1000.;

	    //	cout<< "The decay mode is [LoadTable] "<<theDecayMode<<G4endl;

	    switch (theDecayMode) {
	    case IT:
	      {
		//
		//
		// Decay mode is isomeric transition.
		//
		G4ITDecayChannel *anITChannel = new G4ITDecayChannel
		  (GetVerboseLevel(), (const G4Ions*) &theParentNucleus, b);
		anITChannel->SetICM(applyICM);
		anITChannel->SetARM(applyARM);
		anITChannel->SetHLThreshold(halflifethreshold);
		theDecayTable->Insert(anITChannel);
		break;
	      }
	    case BetaMinus:
	      {
		//
		//
		// Decay mode is beta-.
		//
		if (modeFirstRecord[1])
		  {modeFirstRecord[1] = false; modeTotalBR[1] = b;}
		else {
		  if (c > 0.) {
		    // to work out the Fermi function normalization factor first
		    G4BetaFermiFunction* aBetaFermiFunction = new G4BetaFermiFunction (A, (Z+1));
		    e0 = c*MeV/0.511;
		    n = aBetaFermiFunction->GetFFN(e0);

		    // now to work out the histogram and initialise the random generator
		    G4int npti = 100;				
		    G4double* pdf = new G4double[npti];
		    G4int ptn;
		    G4double g,e,ee,f;
		    ee = e0+1.;
		    for (ptn=0; ptn<npti; ptn++) {
		      // e =e0*(ptn+1.)/102.;
		      // bug fix (#662) (flei, 22/09/2004)
		      e =e0*(ptn+0.5)/100.;
		      g = e+1.;
		      f = std::sqrt(g*g-1)*(ee-g)*(ee-g)*g;
		      // Special treatment for K-40 (problem #1068) (flei,06/05/2010) 
		      if (theParentNucleus.GetParticleName() == "K40[0.0]") f *=
									      (std::pow((g*g-1),3)+std::pow((ee-g),6)+7*(g*g-1)*std::pow((ee-g),2)*(g*g-1+std::pow((ee-g),2)));
		      pdf[ptn] = f*aBetaFermiFunction->GetFF(e);
		    }		  
		    RandGeneral* aRandomEnergy = new RandGeneral( pdf, npti);  
		    G4BetaMinusDecayChannel *aBetaMinusChannel = new
		      G4BetaMinusDecayChannel (GetVerboseLevel(), &theParentNucleus,
					       b, c*MeV, a*MeV, n, FBeta, aRandomEnergy);
		    aBetaMinusChannel->SetICM(applyICM);
		    aBetaMinusChannel->SetARM(applyARM);
		    aBetaMinusChannel->SetHLThreshold(halflifethreshold);
		    theDecayTable->Insert(aBetaMinusChannel);
		    modeSumBR[1] += b;
		    delete[] pdf;
		    delete aBetaFermiFunction;
		  }
		}
		break;
	      case BetaPlus:
		//
		//
		// Decay mode is beta+.
		//
		if (modeFirstRecord[2])
		  {modeFirstRecord[2] = false; modeTotalBR[2] = b;}
		else {
		  //  e0 = c*MeV/0.511;
		  // bug fix (#662) (flei, 22/09/2004)
		  // need to test e0 as there are some data files (e.g. z67.a162) which have entries for beta+
		  // with Q < 2Me
		  //
		  e0 = c*MeV/0.511 -2.;
		  if (e0 > 0.) {
		    G4BetaFermiFunction* aBetaFermiFunction = new G4BetaFermiFunction (A, -(Z-1));

		    n = aBetaFermiFunction->GetFFN(e0);

		    // now to work out the histogram and initialise the random generator
		    G4int npti = 100;				
		    G4double* pdf = new G4double[npti];
		    G4int ptn;
		    G4double g,e,ee,f;
		    ee = e0+1.;
		    for (ptn=0; ptn<npti; ptn++) {
		      // e =e0*(ptn+1.)/102.;
		      // bug fix (#662) (flei, 22/09/2004)
		      e =e0*(ptn+0.5)/100.;
		      g = e+1.;
		      f = std::sqrt(g*g-1)*(ee-g)*(ee-g)*g;
		      pdf[ptn] = f*aBetaFermiFunction->GetFF(e);
		    }
		    RandGeneral* aRandomEnergy = new RandGeneral( pdf, npti);  
		    G4BetaPlusDecayChannel *aBetaPlusChannel = new 
		      G4BetaPlusDecayChannel (GetVerboseLevel(), &theParentNucleus,
					      b, c*MeV, a*MeV, n, FBeta, aRandomEnergy);
		    aBetaPlusChannel->SetICM(applyICM);
		    aBetaPlusChannel->SetARM(applyARM);
		    aBetaPlusChannel->SetHLThreshold(halflifethreshold);
		    theDecayTable->Insert(aBetaPlusChannel);
		    modeSumBR[2] += b;

		    delete[] pdf;
		    delete aBetaFermiFunction;	      
		  }
		}
		break;
	      case KshellEC:
		//
		//
		// Decay mode is K-electron capture.
		//
		if (modeFirstRecord[3])
		  {modeFirstRecord[3] = false; modeTotalBR[3] = b;}
		else {
		  G4KshellECDecayChannel *aKECChannel = new
		    G4KshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
		  aKECChannel->SetICM(applyICM);
		  aKECChannel->SetARM(applyARM);
		  aKECChannel->SetHLThreshold(halflifethreshold);
		  theDecayTable->Insert(aKECChannel);
		  //delete aKECChannel;
		  modeSumBR[3] += b;
		}
		break;
	      case LshellEC:
		//
		//
		// Decay mode is L-electron capture.
		//
		if (modeFirstRecord[4])
		  {modeFirstRecord[4] = false; modeTotalBR[4] = b;}
		else {
		  G4LshellECDecayChannel *aLECChannel = new
		    G4LshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
		  aLECChannel->SetICM(applyICM);
		  aLECChannel->SetARM(applyARM);
		  aLECChannel->SetHLThreshold(halflifethreshold);
		  theDecayTable->Insert(aLECChannel);
		  //delete aLECChannel;
		  modeSumBR[4] += b;
		}
		break;
	      case MshellEC:
		//
		//
		// Decay mode is M-electron capture. In this implementation it is added to L-shell case
		//
		if (modeFirstRecord[5])
		  {modeFirstRecord[5] = false; modeTotalBR[5] = b;}
		else {
		  G4MshellECDecayChannel *aMECChannel = new
		    G4MshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
		  aMECChannel->SetICM(applyICM);
		  aMECChannel->SetARM(applyARM);
		  aMECChannel->SetHLThreshold(halflifethreshold);
		  theDecayTable->Insert(aMECChannel);
		  modeSumBR[5] += b;
		}
		break;
	      case Alpha:
		//
		//
		// Decay mode is alpha.
		//
		if (modeFirstRecord[6])
		  {modeFirstRecord[6] = false; modeTotalBR[6] = b;}
		else {
		  G4AlphaDecayChannel *anAlphaChannel = new
		    G4AlphaDecayChannel (GetVerboseLevel(), &theParentNucleus,
					 b, c*MeV, a*MeV);
		  anAlphaChannel->SetICM(applyICM);
		  anAlphaChannel->SetARM(applyARM);
		  anAlphaChannel->SetHLThreshold(halflifethreshold);
		  theDecayTable->Insert(anAlphaChannel);
		  //	      delete anAlphaChannel;
		  modeSumBR[6] += b;
		}
		break;
	      case ERROR:
	      default:
		G4Exception("G4RadioactiveDecay::LoadDecayTable()", "601",
			    FatalException, "Error in decay mode selection");

	      }
	    }
	  }
	}
      }
    }

    // Go through the decay table and make sure that the branching ratios are
    // correctly normalised.

    G4VDecayChannel* theChannel = 0;
    G4NuclearDecayChannel* theNuclearDecayChannel = 0;
    G4String mode = "";

    // G4int j = 0;  DHW 19 May 2011: variable set but not used

    G4double theBR = 0.0;
    for (G4int i = 0; i < theDecayTable->entries(); i++) {
      theChannel = theDecayTable->GetDecayChannel(i);
      theNuclearDecayChannel = static_cast<G4NuclearDecayChannel*>(theChannel);
      theDecayMode = theNuclearDecayChannel->GetDecayMode();

      // j = 0;   DHW 19 May 2011: variable set but not used

      if (theDecayMode != IT) {
	theBR = theChannel->GetBR();
	theChannel->SetBR(theBR*modeTotalBR[theDecayMode]/modeSumBR[theDecayMode]);
      }
    }
 
  }   // if (DecaySchemeFile)	
  DecaySchemeFile.close();

		
  if (!found && E > 0.) {
    // cases where IT cascade for exited isotopes without entry in RDM database
    // Decay mode is isomeric transition.
    //
    G4ITDecayChannel *anITChannel = new G4ITDecayChannel
      (GetVerboseLevel(), (const G4Ions*) &theParentNucleus, 1);
    anITChannel->SetICM(applyICM);
    anITChannel->SetARM(applyARM);
    anITChannel->SetHLThreshold(halflifethreshold);
    theDecayTable->Insert(anITChannel);
  } 
  if (!theDecayTable) {
    //
    // There is no radioactive decay data for this nucleus.  Return a null
    // decay table.
    //
    G4cerr <<"G4RadoactiveDecay::LoadDecayTable() : cannot find ion radioactive decay file " <<G4endl;
    theDecayTable = 0;
    return theDecayTable;
  }	
  if (theDecayTable && GetVerboseLevel()>1)
    {
      G4cout <<"G4RadioactiveDecay::LoadDecayTable()" << G4endl;
      G4cout << "  No. of  entries: "<< theDecayTable->entries() <<G4endl;
      theDecayTable ->DumpInfo();
    }

  return theDecayTable;
}


void
G4RadioactiveDecay::SetDecayRate(G4int theZ, G4int theA, G4double theE, 
                                 G4int theG, std::vector<G4double> theRates, 
                                 std::vector<G4double> theTaos)
{ 
  //fill the decay rate vector 
  theDecayRate.SetZ(theZ);
  theDecayRate.SetA(theA);
  theDecayRate.SetE(theE);
  theDecayRate.SetGeneration(theG);
  theDecayRate.SetDecayRateC(theRates);
  theDecayRate.SetTaos(theTaos);
}

//////////////////////////////////////////////////////////////////////////
// 
void
G4RadioactiveDecay::AddDecayRateTable(const G4ParticleDefinition& theParentNucleus)
{
  // 1) To calculate all the coefficiecies required to derive the
  //    radioactivities for all progeny of theParentNucleus
  //
  // 2) Add the coefficiencies to the decay rate table vector 
  //

  //
  // Create and initialise variables used in the method.
  //

  theDecayRateVector.clear();

  G4int nGeneration = 0;
  std::vector<G4double> rates;
  std::vector<G4double> taos;

  // start rate is -1.
  // Eq.4.26 of the Technical Note
  rates.push_back(-1.);
  //
  //
  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();
  G4double E = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4double tao = theParentNucleus.GetPDGLifeTime();
  if (tao < 0.) tao = 1e-30;
  taos.push_back(tao);
  G4int nEntry = 0;

  //fill the decay rate with the intial isotope data
  SetDecayRate(Z,A,E,nGeneration,rates,taos);

  // store the decay rate in decay rate vector
  theDecayRateVector.push_back(theDecayRate);
  nEntry++;

  // now start treating the sencondary generations..

  G4bool stable = false;
  G4int i;
  G4int j;
  G4VDecayChannel* theChannel = 0;
  G4NuclearDecayChannel* theNuclearDecayChannel = 0;
  G4ITDecayChannel* theITChannel = 0;
  G4BetaMinusDecayChannel *theBetaMinusChannel = 0;
  G4BetaPlusDecayChannel *theBetaPlusChannel = 0;
  G4AlphaDecayChannel *theAlphaChannel = 0;
  G4RadioactiveDecayMode theDecayMode;
  G4double theBR = 0.0;
  G4int AP = 0;
  G4int ZP = 0;
  G4int AD = 0;
  G4int ZD = 0;
  G4double EP = 0.;
  std::vector<G4double> TP;
  std::vector<G4double> RP;
  G4ParticleDefinition *theDaughterNucleus;
  G4double daughterExcitation;
  G4ParticleDefinition *aParentNucleus;
  G4IonTable* theIonTable;
  G4DecayTable *aTempDecayTable;
  G4double theRate;
  G4double TaoPlus;
  G4int nS = 0;
  G4int nT = nEntry;
  G4double brs[7];
  //
  theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  while (!stable) {
    nGeneration++;
    for (j = nS; j< nT; j++) {
      ZP = theDecayRateVector[j].GetZ();
      AP = theDecayRateVector[j].GetA();
      EP = theDecayRateVector[j].GetE();
      RP = theDecayRateVector[j].GetDecayRateC();
      TP = theDecayRateVector[j].GetTaos();      
      if (GetVerboseLevel()>0){
	G4cout <<"G4RadioactiveDecay::AddDecayRateTable : "
	       << " daughters of ("<< ZP <<", "<<AP<<", "
	       << EP <<") "
	       << " are being calculated. "	  
	       <<" generation = "
	       << nGeneration << G4endl;
      }
      aParentNucleus = theIonTable->GetIon(ZP,AP,EP);
      if (!IsLoaded(*aParentNucleus)){
	aParentNucleus->SetDecayTable(LoadDecayTable(*aParentNucleus));
      }
  	
      G4DecayTable* theDecayTable = new G4DecayTable();
      aTempDecayTable = aParentNucleus->GetDecayTable();
      for (i=0; i< 7; i++) brs[i] = 0.0;

      //
      // Go through the decay table and to combine the same decay channels
      //
      for (i=0; i<aTempDecayTable->entries(); i++) {
	theChannel             = aTempDecayTable->GetDecayChannel(i);
	theNuclearDecayChannel = static_cast<G4NuclearDecayChannel *>(theChannel);
	theDecayMode           = theNuclearDecayChannel->GetDecayMode();
	daughterExcitation = theNuclearDecayChannel->GetDaughterExcitation ();
	theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus () ;
	AD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
	ZD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();  
	G4NuclearLevelManager* levelManager = 
           G4NuclearLevelStore::GetInstance()->GetManager(ZD, AD);
	if (levelManager->NumberOfLevels() ) {
	  const G4NuclearLevel* level =
              levelManager->NearestLevel (daughterExcitation);

	  if (std::abs(daughterExcitation - level->Energy()) < levelTolerance) {
	    // Level half-life is in ns and the threshold is set to 1 micros by default, user can set it via the UI command
	    if (level->HalfLife()*ns >= halflifethreshold ){    
	      // save the metastable nucleus 
	      theDecayTable->Insert(theChannel);
	    } 
	    else{
	      brs[theDecayMode] += theChannel->GetBR();
	    }
	  }
	  else {
	    brs[theDecayMode] += theChannel->GetBR();
	  }
	}
	else{
	  brs[theDecayMode] += theChannel->GetBR();
	}
      }	    
      brs[2] = brs[2]+brs[3]+brs[4]+brs[5];
      brs[3] = brs[4] =brs[5] =  0.0;
      for (i= 0; i<7; i++){
	if (brs[i] > 0.) {
	  switch ( i ) {
	  case 0:
	    //
	    //
	    // Decay mode is isomeric transition.
	    //

	    theITChannel =  new G4ITDecayChannel
	      (0, (const G4Ions*) aParentNucleus, brs[0]);

	    theDecayTable->Insert(theITChannel);
	    break;

	  case 1:
	    //
	    //
	    // Decay mode is beta-.
	    //
	    theBetaMinusChannel = new G4BetaMinusDecayChannel (0, aParentNucleus,
							       brs[1], 0.*MeV, 0.*MeV, 1, false, 0);
	    theDecayTable->Insert(theBetaMinusChannel);

	    break;

	  case 2:
	    //
	    //
	    // Decay mode is beta+ + EC.
	    //
	    theBetaPlusChannel = new G4BetaPlusDecayChannel (GetVerboseLevel(), aParentNucleus,
							     brs[2], 0.*MeV, 0.*MeV, 1, false, 0);
	    theDecayTable->Insert(theBetaPlusChannel);
	    break;		      

	  case 6:
	    //
	    //
	    // Decay mode is alpha.
	    //
	    theAlphaChannel = new G4AlphaDecayChannel(GetVerboseLevel(),
                                                      aParentNucleus,
						      brs[6], 0.*MeV, 0.*MeV);
	    theDecayTable->Insert(theAlphaChannel);
	    break;

	  default:
	    break;
	  }
	}
      }
      // 
      // loop over all branches in theDecayTable
      //
      for (i = 0; i < theDecayTable->entries(); i++){
	theChannel = theDecayTable->GetDecayChannel(i);
	theNuclearDecayChannel = static_cast<G4NuclearDecayChannel*>(theChannel);
	theBR = theChannel->GetBR();
	theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus();
	//  first check if the decay of the original nucleus is an IT channel, if true create a new groud-level nucleus
	if (theNuclearDecayChannel->GetDecayMode() == IT && nGeneration == 1) {
	  A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
	  Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
	  theDaughterNucleus=theIonTable->GetIon(Z,A,0.);
        }
        if (IsApplicable(*theDaughterNucleus) &&
            theBR && 
            aParentNucleus != theDaughterNucleus) { 
	  // need to make sure daugher has decaytable
	  if (!IsLoaded(*theDaughterNucleus)){
	    theDaughterNucleus->SetDecayTable(LoadDecayTable(*theDaughterNucleus));
	  }
	  if (theDaughterNucleus->GetDecayTable()->entries() ) {
	    //
	    A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
	    Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
	    E = ((const G4Ions*)(theDaughterNucleus))->GetExcitationEnergy();

	    TaoPlus = theDaughterNucleus->GetPDGLifeTime();
	    //		cout << TaoPlus <<G4endl;
	    if (TaoPlus <= 0.)  TaoPlus = 1e-30; 

	    // first set the taos, one simply need to add to the parent ones
	    taos.clear();
	    taos = TP;
	    taos.push_back(TaoPlus);
	    // now calculate the coefficiencies
	    //
	    // they are in two parts, first the less than n ones
	    // Eq 4.24 of the TN
	    rates.clear();
	    size_t k;
	    for (k = 0; k < RP.size(); k++){
	      if ((TP[k]-TaoPlus) == 0.) {
		theRate =  1e30;
	      } else {
		theRate = TP[k]/(TP[k]-TaoPlus) * theBR * RP[k];
	      }
	      rates.push_back(theRate);
	    }
	    //
	    // the sencond part: the n:n coefficiency
	    // Eq 4.25 of the TN.  Note Yn+1 is zero apart from Y1 which is -1 as treated at line 1013
	    // 
	    theRate = 0.;
	    G4double aRate;
	    for (k = 0; k < RP.size(); k++){
	      if ((TP[k]-TaoPlus) == 0.) {
		aRate =  1e30;
	      } else {
		aRate = TaoPlus/(TP[k]-TaoPlus) * theBR * RP[k];
	      }
	      theRate -= aRate;
	    }
	    rates.push_back(theRate); 	      
	    SetDecayRate (Z,A,E,nGeneration,rates,taos);
	    theDecayRateVector.push_back(theDecayRate);
	    nEntry++;
	  }
	} // end of testing daughter nucleus
      } // end of i loop( the branches) 
      delete theDecayTable;

    } //end of for j loop
    nS = nT;
    nT = nEntry;
    if (nS == nT) stable = true;
  }

  //end of while loop
  // the calculation completed here


  // fill the first part of the decay rate table
  // which is the name of the original particle (isotope) 
  //
  theDecayRateTable.SetIonName(theParentNucleus.GetParticleName()); 
  //
  //
  // now fill the decay table with the newly completed decay rate vector
  //

  theDecayRateTable.SetItsRates(theDecayRateVector);

  //
  // finally add the decayratetable to the tablevector
  //
  theDecayRateTableVector.push_back(theDecayRateTable);
}

///////////////////////////////////////////////////////////////////////////////
//
// SetSourceTimeProfile
//
//  read in the source time profile function (histogram)
//

void G4RadioactiveDecay::SetSourceTimeProfile(G4String filename)
{
  std::ifstream infile ( filename, std::ios::in );
  if ( !infile ) G4Exception(__FILE__, G4inttostring(__LINE__), FatalException,  "Unable to open source data file" );

  G4double bin, flux;
  NSourceBin = -1;
  while (infile >> bin >> flux ) {
    NSourceBin++;
    if (NSourceBin > 99) {
      G4Exception(__FILE__, G4inttostring(__LINE__), FatalException,
                  "input source data file too big (>100 rows)" );
    } else {
      SBin[NSourceBin] = bin * s;
      SProfile[NSourceBin] = flux;
    }
  }
  SetAnalogueMonteCarlo(0);
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Source Timeprofile Nbin = " << NSourceBin <<G4endl;}
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
// SetDecayBiasProfile
//
// read in the decay bias scheme function (histogram)
//
void G4RadioactiveDecay::SetDecayBias(G4String filename)
{
  std::ifstream infile ( filename, std::ios::in);
  if ( !infile ) G4Exception(__FILE__, G4inttostring(__LINE__), FatalException,  "Unable to open bias data file" );

  G4double bin, flux;
  G4int dWindows = 0;
  G4int i ;

  theRadioactivityTables.clear();
  //  for (i = 0; i<100; i++) decayWindows[i] = -1;

  NDecayBin = -1;
  while (infile >> bin >> flux ) {
    NDecayBin++;
    if (NDecayBin > 99) {
      G4Exception(__FILE__, G4inttostring(__LINE__), FatalException,
                  "input bias data file too big (>100 rows)" );
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
  //
  SetAnalogueMonteCarlo(0);
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Decay Bias Profile  Nbin = " << NDecayBin <<G4endl;}
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
// DecayIt
//
G4VParticleChange*
G4RadioactiveDecay::DecayIt(const G4Track& theTrack, const G4Step&)
{
  // Initialize the G4ParticleChange object. Get the particle details and the
  // decay table.

  fParticleChangeForRadDecay.Initialize(theTrack);
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  G4ParticleDefinition *theParticleDef = theParticle->GetDefinition();

  // First check whether RDM applies to the current logical volume

  if (!std::binary_search(ValidVolumes.begin(), ValidVolumes.end(), 
			  theTrack.GetVolume()->GetLogicalVolume()->GetName())) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout <<"G4RadioactiveDecay::DecayIt : "
             << theTrack.GetVolume()->GetLogicalVolume()->GetName()
             << " is not selected for the RDM"<< G4endl;
      G4cout << " There are " << ValidVolumes.size() << " volumes" << G4endl;
      G4cout << " The Valid volumes are " << G4endl;
      for (size_t i = 0; i< ValidVolumes.size(); i++) G4cout << ValidVolumes[i] << G4endl;
    }
#endif
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle.

    fParticleChangeForRadDecay.ProposeTrackStatus( fStopAndKill ) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;
  }

  // now check is the particle is valid for RDM

  if (!(IsApplicable(*theParticleDef))) { 
    //
    // The particle is not a Ion or outside the nucleuslimits for decay
    //
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr <<"G4RadioactiveDecay::DecayIt : "
	     <<theParticleDef->GetParticleName() 
	     << " is not a valid nucleus for the RDM"<< G4endl;
    }
#endif
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    //
    // Kill the parent particle.
    //
    fParticleChangeForRadDecay.ProposeTrackStatus( fStopAndKill ) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;
  }

  if (!IsLoaded(*theParticleDef))
    theParticleDef->SetDecayTable(LoadDecayTable(*theParticleDef));

  G4DecayTable* theDecayTable = theParticleDef->GetDecayTable();

  if  (theDecayTable == 0 || theDecayTable->entries() == 0 )
    {
      //
      //
      // There are no data in the decay table.  Set the particle change parameters
      // to indicate this.
      //
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	{
	  G4cerr <<"G4RadioactiveDecay::DecayIt : decay table not defined  for ";
	  G4cerr <<theParticleDef->GetParticleName() <<G4endl;
	}
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      //
      //
      // Kill the parent particle.
      //
      fParticleChangeForRadDecay.ProposeTrackStatus( fStopAndKill ) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }
  else 
    { 
      //
      // now try to  decay it
      //
      G4double energyDeposit = 0.0;
      G4double finalGlobalTime = theTrack.GetGlobalTime();
      G4int index;
      G4ThreeVector currentPosition;
      currentPosition = theTrack.GetPosition();

      // check whether use Analogue or VR implementation
      //
      if (AnalogueMC){
	//
	// Aanalogue MC 
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0)
	  {
	    G4cout <<"DecayIt:  Analogue MC version " << G4endl;
	  }
#endif
	//
	G4DecayProducts *products = DoDecay(*theParticleDef);
	//
	// check if the product is the same as input and kill the track if necessary to prevent infinite loop
	// (11/05/10, F.Lei)
	//
	if ( products->entries() == 1) {
	  fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
	  fParticleChangeForRadDecay.ProposeTrackStatus( fStopAndKill ) ;
	  fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
	  ClearNumberOfInteractionLengthLeft();
	  return &fParticleChangeForRadDecay;
	}	
	//
	// Get parent particle information and boost the decay products to the
	// laboratory frame based on this information.
	//
	G4double ParentEnergy = theParticle->GetTotalEnergy();
	G4ThreeVector ParentDirection(theParticle->GetMomentumDirection());

	if (theTrack.GetTrackStatus() == fStopButAlive )
	  {
	    //
	    //
	    // The particle is decayed at rest.
	    //
	    // since the time is still for rest particle in G4 we need to add the additional 
	    // time lapsed between the particle come to rest and the actual decay. This time 
	    // is simply sampled with the mean-life of the particle.
	    // but we need to protect the case PDGTime < 0. (F.Lei 11/05/10)
	    G4double temptime =  -std::log( G4UniformRand()) * theParticleDef->GetPDGLifeTime();
	    if (temptime <0.) temptime =0.; 
	    finalGlobalTime += temptime ;
	    energyDeposit += theParticle->GetKineticEnergy();
	  }
	else
	  {
	    //
	    //
	    // The particle is decayed in flight (PostStep case).
	    //
	    products->Boost( ParentEnergy, ParentDirection);
	  }
	//
	//
	// Add products in theParticleChangeForRadDecay.
	//
	G4int numberOfSecondaries = products->entries();
	fParticleChangeForRadDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
	if (GetVerboseLevel()>1) {
	  G4cout <<"G4RadioactiveDecay::DecayIt : Decay vertex :";
	  G4cout <<" Time: " <<finalGlobalTime/ns <<"[ns]";
	  G4cout <<" X:" <<(theTrack.GetPosition()).x() /cm <<"[cm]";
	  G4cout <<" Y:" <<(theTrack.GetPosition()).y() /cm <<"[cm]";
	  G4cout <<" Z:" <<(theTrack.GetPosition()).z() /cm <<"[cm]";
	  G4cout <<G4endl;
	  G4cout <<"G4Decay::DecayIt  : decay products in Lab. Frame" <<G4endl;
	  products->DumpInfo();
	}
#endif
	for (index=0; index < numberOfSecondaries; index++) 
	  {
	    G4Track* secondary = new G4Track
	      (products->PopProducts(), finalGlobalTime, currentPosition);
	    secondary->SetGoodForTrackingFlag();
	    secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
	    fParticleChangeForRadDecay.AddSecondary(secondary);
	  }
	delete products;
	//
	// end of analogue MC algarithm
	//
      }
      else {
	//
	// Varaice Reduction Method
	//
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0)
	  {
	    G4cout <<"DecayIt:  Variance Reduction version " << G4endl;
	  }
#endif
	if (!IsRateTableReady(*theParticleDef)) {
	  // if the decayrates are not ready, calculate them and 
	  // add to the rate table vector 
	  AddDecayRateTable(*theParticleDef);
	}
	//retrieve the rates 
	GetDecayRateTable(*theParticleDef);
	//
	// declare some of the variables required in the implementation
	//
	G4ParticleDefinition* parentNucleus;
	G4IonTable *theIonTable;
	G4int PZ;
	G4int PA;
	G4double PE;
	std::vector<G4double> PT;
	std::vector<G4double> PR;
	G4double taotime;
	G4double decayRate;

	size_t i;
	size_t j;
	G4int numberOfSecondaries;
	G4int totalNumberOfSecondaries = 0;
	G4double currentTime = 0.;
	G4int ndecaych;
	G4DynamicParticle* asecondaryparticle;
	//	G4DecayProducts* products = 0;
	std::vector<G4DynamicParticle*> secondaryparticles;
	std::vector<G4double> pw;
	std::vector<G4double> ptime;
	pw.clear();
	ptime.clear();
	//now apply the nucleus splitting
	//
	//
	for (G4int n = 0; n < NSplit; n++)
	  {
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

            for (i = 0; i<theDecayRateVector.size(); i++){
	      PZ = theDecayRateVector[i].GetZ();
	      PA = theDecayRateVector[i].GetA();
	      PE = theDecayRateVector[i].GetE();
	      PT = theDecayRateVector[i].GetTaos();
	      PR = theDecayRateVector[i].GetDecayRateC();

	      // calculate the decay rate of the isotope
	      // decayRate is the radioactivity of isotope (PZ,PA,PE) at the 
	      // time 'theDecayTime'
	      // it will be used to calculate the statistical weight of the 
	      // decay products of this isotope

	      //  G4cout <<"PA= "<< PA << " PZ= " << PZ << " PE= "<< PE  <<G4endl;
	      decayRate = 0.;
	      for ( j = 0; j < PT.size(); j++){
		taotime = GetTaoTime(theDecayTime,PT[j]);
		decayRate -= PR[j] * taotime;
		// Eq.4.23 of of the TN
		// note the negative here is required as the rate in the
                // eqation is defined to be negative, 
		// i.e. decay away, but we need pasitive value here.

		//  G4cout << j << "\t"<< PT[j]/s <<"\t"<<PR[j]<< "\t" << decayRate << G4endl ;		
	      }

	      // now calculate the statistical weight

	      // one need to fold the the source bias function with the
              // decaytime
              // also need to include the track weight! (F.Lei, 28/10/10)
	      G4double weight = weight1*decayRate*theTrack.GetWeight(); 
		
	      // add the isotope to the radioactivity tables
	      //  G4cout <<theDecayTime/s <<"\t"<<nbin<<G4endl;
	      //  G4cout << theTrack.GetWeight() <<"\t"<<weight1<<"\t"<<decayRate<< G4endl;
	      theRadioactivityTables[decayWindows[nbin-1]]->AddIsotope(PZ,PA,PE,weight);
				
	      // decay the isotope 
	      theIonTable = (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
	      parentNucleus = theIonTable->GetIon(PZ,PA,PE);

              // Create a temprary products buffer.
              // its contents to be transfered to the products at the end of the loop
              G4DecayProducts* tempprods = 0;

              // Decide whether to apply branching ratio bias or not	     
              if (BRBias){
		G4DecayTable* theDecayTable = parentNucleus->GetDecayTable();
		ndecaych = G4int(theDecayTable->entries()*G4UniformRand());
		G4VDecayChannel* theDecayChannel = theDecayTable->GetDecayChannel(ndecaych);
		if (theDecayChannel == 0) {
		  // Decay channel not found.
#ifdef G4VERBOSE
		  if (GetVerboseLevel()>0) {
		    G4cerr << " G4RadioactiveDecay::DoIt : cannot determine decay channel ";
                    G4cerr << " for this nucleus; decay as if no biasing active ";
		    G4cerr << G4endl;
		    theDecayTable ->DumpInfo();
		  }
#endif
                  tempprods = DoDecay(*parentNucleus);  // DHW 6 Dec 2010 - do decay as if no biasing
                                                        //           to avoid deref of temppprods = 0
		} else {
		  // A decay channel has been identified, so execute the DecayIt.
		  G4double tempmass = parentNucleus->GetPDGMass();      
		  tempprods = theDecayChannel->DecayIt(tempmass);
		  weight *= (theDecayChannel->GetBR())*(theDecayTable->entries());
		}
              } else {
                tempprods = DoDecay(*parentNucleus);
              }

	      // save the secondaries for buffers

              numberOfSecondaries = tempprods->entries();
	      currentTime = finalGlobalTime + theDecayTime;
	      for (index = 0; index < numberOfSecondaries; index++) {
		asecondaryparticle = tempprods->PopProducts();
		if (asecondaryparticle->GetDefinition()->GetBaryonNumber() < 5){
		  pw.push_back(weight);
		  ptime.push_back(currentTime);
		  secondaryparticles.push_back(asecondaryparticle);
		}
	      }

	      delete tempprods;

	      //end of i loop
            }

	    // end of n loop 
	  }

	// now deal with the secondaries in the two stl containers
	// and submmit them back to the tracking manager
	//
	totalNumberOfSecondaries = pw.size();
	fParticleChangeForRadDecay.SetNumberOfSecondaries(totalNumberOfSecondaries);
	for (index=0; index < totalNumberOfSecondaries; index++) 
	  { 
	    G4Track* secondary = new G4Track(
					     secondaryparticles[index], ptime[index], currentPosition);
	    secondary->SetGoodForTrackingFlag(); 	   
	    secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
	    secondary->SetWeight(pw[index]); 	   
	    fParticleChangeForRadDecay.AddSecondary(secondary); 
	  }
	//
	// make sure the original track is set to stop and its kinematic energy collected
	// 
	//theTrack.SetTrackStatus(fStopButAlive);
	//energyDeposit += theParticle->GetKineticEnergy();

      }

      //
      // Kill the parent particle.
      //
      fParticleChangeForRadDecay.ProposeTrackStatus( fStopAndKill ) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(energyDeposit);
      // 
      fParticleChangeForRadDecay.ProposeGlobalTime( finalGlobalTime );
      //
      // Reset NumberOfInteractionLengthLeft.
      //
      ClearNumberOfInteractionLengthLeft();

      return &fParticleChangeForRadDecay ;
    }
} 

///////////////////////////////////////////////////////////////////
//
// DoDecay
//
G4DecayProducts*
G4RadioactiveDecay::DoDecay(  G4ParticleDefinition& theParticleDef )
{
  G4DecayProducts* products = 0;

  // follow the decaytable and generate the secondaries...

#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) G4cout<<"Begin of DoDecay..."<<G4endl;
#endif

  G4DecayTable* theDecayTable = theParticleDef.GetDecayTable();

  // Choose a decay channel.

#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) G4cout <<"Selecte a channel..."<<G4endl;
#endif

  G4VDecayChannel* theDecayChannel = theDecayTable->SelectADecayChannel();
  if (theDecayChannel == 0) {
    // Decay channel not found.

    G4cerr <<"G4RadioactiveDecay::DoIt : can not determine decay channel";
    G4cerr <<G4endl;
    theDecayTable ->DumpInfo();
  } else {

    // A decay channel has been identified, so execute the DecayIt.

#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cerr <<"G4RadioactiveDecay::DoIt : selected decay channel  addr:";
      G4cerr <<theDecayChannel <<G4endl;
    }
#endif

    G4double tempmass = theParticleDef.GetPDGMass();
    products = theDecayChannel->DecayIt(tempmass);
  }

  return products;
}
