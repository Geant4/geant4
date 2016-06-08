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
//      http://www.space.dera.gov.uk/space_env/rdm.html
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
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
//            2) VR: Significant efficiency inprovement
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
#include "G4AlphaDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4RIsotopeTable.hh"
#include "G4BetaFermiFunction.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4DiscreteGammaDeexcitation.hh"

#include "g4std/vector"
#include "g4std/strstream"
#include "g4std/algorithm"
#include "g4std/fstream"

const G4double   G4RadioactiveDecay::levelTolerance =2.0*keV;

///////////////////////////////////////////////////////////////////////////////
//
//
// Constructor
//
G4RadioactiveDecay::G4RadioactiveDecay
  (const G4String& processName)
  :G4VRestDiscreteProcess(processName, fDecay), HighestBinValue(10.0),
   LowestBinValue(1.0e-3), TotBin(200), verboseLevel(1)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout <<"G4RadioactiveDecay constructor    Name: ";
    G4cout <<processName << G4endl;   }
#endif

  theRadioactiveDecaymessenger = new G4RadioactiveDecaymessenger(this);
  theIsotopeTable              = new G4RIsotopeTable();
  aPhysicsTable                = NULL;
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
  SBin[1]     = 1e10 * s;
  SProfile[0] = 1.;
  SProfile[1] = 1.;
  NDecayBin   = 1;
  DBin[0]     = (1e10 -1.) * s ;
  DBin[1]     = 1e10 * s;
  DProfile[0] = 1.;
  DProfile[1] = 0.;
  NSplit      = 1;
  AnalogueMC  = true ;
  FBeta       = false ;
  BRBias      = true ;
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
  if (aPhysicsTable != NULL)
    {
      aPhysicsTable->clearAndDestroy();
      delete aPhysicsTable;
    }
  //  delete theIsotopeTable;
  delete theRadioactiveDecaymessenger;
}

///////////////////////////////////////////////////////////////////////////////
//
//
// IsApplicable
//
G4bool G4RadioactiveDecay::IsApplicable(const G4ParticleDefinition &
  aParticle)
{
  //
  //
  // All particles, other than G4Ions, are rejected by default.
  //
  if (!(aParticle.GetParticleType() == "nucleus")) {return false;}
  else if (aParticle.GetPDGLifeTime() < 0. &&
	   aParticle.GetParticleName() != "GenericIon") {
    return false;
  }
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
  return G4std::binary_search(LoadedNuclei.begin(),
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
  for (G4int i = 0; i < theLogicalVolumes->entries(); i++){
    volume=theLogicalVolumes->operator()(i);
    if (volume->GetName() == aVolume) {
      ValidVolumes.push_back(aVolume);
      G4std::sort(ValidVolumes.begin(), ValidVolumes.end());
      // sort need for performing binary_search
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	G4cout << " RDM Applies to : " << aVolume << G4endl; 
#endif
    }else if(i ==  theLogicalVolumes->entries())
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
  for (G4int i = 0; i < theLogicalVolumes->entries(); i++){
    volume=theLogicalVolumes->operator()(i);
    if (volume->GetName() == aVolume) {
      G4std::vector<G4String>::iterator location;
      location = G4std::find(ValidVolumes.begin(),ValidVolumes.end(),aVolume);
      if (location != ValidVolumes.end()) {
	ValidVolumes.erase(location);
	G4std::sort(ValidVolumes.begin(), ValidVolumes.end());
      }else{
	G4cerr << " DeselectVolume:" << aVolume << " is not in the list"<< G4endl; 
      }	  
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	G4cout << " DeselectVolume: " << aVolume << " is removed from list"<<G4endl; 
#endif
    }else if(i ==  theLogicalVolumes->entries())
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
  for (G4int i = 0; i < theLogicalVolumes->entries(); i++){
    volume=theLogicalVolumes->operator()(i);
    ValidVolumes.push_back(volume->GetName());    
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout << "         RDM Applies to Volume "  << volume->GetName() << G4endl;
#endif
  }
  G4std::sort(ValidVolumes.begin(), ValidVolumes.end());
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
  for (G4int i = 0; i < theDecayRateTableVector.size(); i++)
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

  for (G4int i = 0; i < theDecayRateTableVector.size(); i++)
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
////////////////////////////////////////////////////////////////////////////////
//
//
// GetTaoTime
//
// to perform the convilution of the sourcetimeprofile function with the 
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
	taotime += SProfile[i] * (exp(-(t-SBin[i+1])/tao)-exp(-(t-SBin[i])/tao));
      }
  }
  taotime +=  SProfile[nbin] * (1-exp(-(t-SBin[nbin])/tao));
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Tao time: " <<taotime <<G4endl;}
#endif
  return  taotime ;
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
    {G4cout <<" Decay time: " <<decaytime <<"[ns]" <<G4endl;}
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
  for (G4int i = 0; i < NDecayBin; i++) 
    {
      if ( aDecayTime > DBin[i]) return i+1;	  
    }
  return  1;
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
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"mean life time: " <<meanlife <<"[ns]" <<G4endl;}
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
       if (aPhysicsTable == NULL) BuildPhysicsTable(*aParticleDef);
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
  if (aPhysicsTable != NULL) return;

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
      beta  = sqrt((1.0 - gammainv)*(1.0 +gammainv));
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
G4DecayTable *G4RadioactiveDecay::LoadDecayTable (G4ParticleDefinition
  &theParentNucleus)
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

  G4String dirName = getenv("G4RADIOACTIVEDATA");
  LoadedNuclei.push_back(theParentNucleus.GetParticleName());
  G4std::sort( LoadedNuclei.begin(), LoadedNuclei.end() );
  // sort needed to allow binary_search

  char val[100];
  G4std::ostrstream os(val,100);
  os <<dirName <<"/z" <<Z <<".a" <<A <<'\0';
  G4String file(val);

  G4std::ifstream DecaySchemeFile(file);

  if (!DecaySchemeFile)
  //
  //
  // There is no radioactive decay data for this nucleus.  Return a null
  // decay table.
  //
  {
    G4cerr <<"G4RadoactiveDecay::LoadDecayTable() : cannot find ion radioactive decay file " <<G4endl;
    theDecayTable = NULL;
    return theDecayTable;
  }
  //
  //
  // Initialise variables used for reading in radioactive decay data.
  //
  G4int    nMode = 7;
  G4bool   modeFirstRecord[7];
  G4double modeTotalBR[7];
  G4double modeSumBR[7];
  G4int i;
  for (i=0; i<nMode; i++)
  {
    modeFirstRecord[i] = true;
    modeSumBR[i]       = 0.0;
  }

  G4bool complete(false);
  G4bool found(false);
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
  while (!complete && -DecaySchemeFile.getline(inputChars, 80).eof() != EOF)
  {
    inputLine = inputChars;
    //    G4String::stripType stripend(1);
    //    inputLine = inputLine.strip(stripend);
    inputLine = inputLine.strip(1);
    if (inputChars[0] != '#' && inputLine.length() != 0)
    {
      G4std::istrstream tmpStream(inputLine);
      if (inputChars[0] == 'P')
  //
  //
  // Nucleus is a parent type.  Check the excitation level to see if it matches
  // that of theParentNucleus
  //
      {
        tmpStream >>recordType >>a >>b;
        if (found) {complete = true;}
        else {found = (abs(a*keV - E)<levelTolerance);}
      }

      else if (found)
      {
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
	else 
	  {
	  tmpStream >>theDecayMode >>a >>b >>c;
	  a/=1000.;
	  c/=1000.;
	  
	  //	cout<< "The decay mode is [LoadTable] "<<theDecayMode<<G4endl;
	  
	  switch (theDecayMode)
	    {
	    case IT:
	      //
	      //
	      // Decay mode is isomeric transition.
	      //
	      {
		G4ITDecayChannel *anITChannel = new G4ITDecayChannel
		  (GetVerboseLevel(), &theParentNucleus, b);
		theDecayTable->Insert(anITChannel);
		break;
	      }
	    case BetaMinus:
	      //
	      //
	      // Decay mode is beta-.
	      //
	      if (modeFirstRecord[1])
		{modeFirstRecord[1] = false; modeTotalBR[1] = b;}
	      else
		{
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
		    e =e0*(ptn+1.)/102.;
		    g = e+1.;
		    f = sqrt(g*g-1)*(ee-g)*(ee-g)*g;
		    pdf[ptn] = f*aBetaFermiFunction->GetFF(e);
		  }		  
		  RandGeneral* aRandomEnergy = new RandGeneral( pdf, npti);  

		  G4BetaMinusDecayChannel *aBetaMinusChannel = new
		    G4BetaMinusDecayChannel (GetVerboseLevel(), &theParentNucleus,
					     b, c*MeV, a*MeV, n, FBeta, aRandomEnergy);
		  theDecayTable->Insert(aBetaMinusChannel);
		  modeSumBR[1] += b;


		  delete[] pdf;
		  delete aBetaFermiFunction;
		}
	      break;
	    case BetaPlus:
	      //
	      //
	      // Decay mode is beta+.
	      //
	      if (modeFirstRecord[2])
		{modeFirstRecord[2] = false; modeTotalBR[2] = b;}
	      else
		{
		  G4BetaFermiFunction* aBetaFermiFunction = new G4BetaFermiFunction (A, -(Z-1));
		  e0 = c*MeV/0.511;
		  n = aBetaFermiFunction->GetFFN(e0);

		  // now to work out the histogram and initialise the random generator
		  G4int npti = 100;				
		  G4double* pdf = new G4double[npti];
		  G4int ptn;
		  G4double g,e,ee,f;
		  ee = e0+1.;
		  for (ptn=0; ptn<npti; ptn++) {
		    e =e0*(ptn+1.)/102.;
		    g = e+1.;
		    f = sqrt(g*g-1)*(ee-g)*(ee-g)*g;
		    pdf[ptn] = f*aBetaFermiFunction->GetFF(e);
		  }		  
		  RandGeneral* aRandomEnergy = new RandGeneral( pdf, npti);  
		  G4BetaPlusDecayChannel *aBetaPlusChannel = new 
		    G4BetaPlusDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV, n, FBeta, aRandomEnergy);
		  theDecayTable->Insert(aBetaPlusChannel);
		  modeSumBR[2] += b;

		  delete[] pdf;
		  delete aBetaFermiFunction;	      
		}
	      break;
	    case KshellEC:
	      //
	      //
	      // Decay mode is K-electron capture.
	      //
	      if (modeFirstRecord[3])
		{modeFirstRecord[3] = false; modeTotalBR[3] = b;}
	      else
		{
		  G4KshellECDecayChannel *aKECChannel = new
		    G4KshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
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
	      else
		{
		  G4LshellECDecayChannel *aLECChannel = new
		    G4LshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
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
	      else
		{
		  G4LshellECDecayChannel *aLECChannel = new
		    G4LshellECDecayChannel (GetVerboseLevel(), &theParentNucleus,
					    b, c*MeV, a*MeV);
		  theDecayTable->Insert(aLECChannel);
		  //delete aLECChannel;
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
	      else
		{
		  G4AlphaDecayChannel *anAlphaChannel = new
		    G4AlphaDecayChannel (GetVerboseLevel(), &theParentNucleus,
					 b, c*MeV, a*MeV);
		  theDecayTable->Insert(anAlphaChannel);
		  //	      delete anAlphaChannel;
		  modeSumBR[6] += b;
		}
	      break;
	    }
	  }
      }
    }
  }  
  DecaySchemeFile.close();

  //
  //
  // Go through the decay table and make sure that the branching ratios are
  // correctly normalised.
  //
  G4VDecayChannel       *theChannel             = NULL;
  G4NuclearDecayChannel *theNuclearDecayChannel = NULL;
  G4String mode                     = "";
  G4int j                           = 0;
  G4double theBR                    = 0.0;
  for (i=0; i<theDecayTable->entries(); i++)
  {
    theChannel             = theDecayTable->GetDecayChannel(i);
    theNuclearDecayChannel = static_cast<G4NuclearDecayChannel *>(theChannel);
    theDecayMode           = theNuclearDecayChannel->GetDecayMode();
    j          = 0;
    if (theDecayMode != IT)
    {
      theBR = theChannel->GetBR();
      theChannel->SetBR(theBR*modeTotalBR[theDecayMode]/modeSumBR[theDecayMode]);
    }
  }  
  return theDecayTable;
}

////////////////////////////////////////////////////////////////////////
//
//
void G4RadioactiveDecay::SetDecayRate(G4int theZ, G4int theA, G4double theE, 
				       G4int theG, G4std::vector<G4double> theRates, 
				       G4std::vector<G4double> theTaos)
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
// 
void G4RadioactiveDecay::AddDecayRateTable(const G4ParticleDefinition &theParentNucleus)
{
  // 1) To calculate all the coefficiecies required to derive the radioactivities for all 
  // progeny of theParentNucleus
  //
  // 2) Add the coefficiencies to the decay rate table vector 
  //
  
  //
  // Create and initialise variables used in the method.
  //

  theDecayRateVector.clear();
  
  G4int nGeneration = 0;
  G4std::vector<G4double> rates;
  G4std::vector<G4double> taos;
  
  // start rate is -1.
  rates.push_back(-1.);
  //
  //
  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();
  G4double E = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4double tao = theParentNucleus.GetPDGLifeTime();
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
  G4int k;
  G4VDecayChannel       *theChannel             = NULL;
  G4NuclearDecayChannel *theNuclearDecayChannel = NULL;
  G4ITDecayChannel *theITChannel = NULL;
  G4BetaMinusDecayChannel *theBetaMinusChannel = NULL;
  G4BetaPlusDecayChannel *theBetaPlusChannel = NULL;
  G4AlphaDecayChannel *theAlphaChannel = NULL;
  G4RadioactiveDecayMode theDecayMode;
  //  G4NuclearLevelManager levelManager;
  //const G4NuclearLevel* level;
  G4double theBR = 0.0;
  G4int AP = 0;
  G4int ZP = 0;
  G4int AD = 0;
  G4int ZD = 0;
  G4double EP = 0.;
  G4std::vector<G4double> TP;
  G4std::vector<G4double> RP;
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
 
  while (!stable) {
    nGeneration++;
    for (j = nS; j< nT; j++) {
      ZP = theDecayRateVector[j].GetZ();
      AP = theDecayRateVector[j].GetA();
      EP = theDecayRateVector[j].GetE();
      RP = theDecayRateVector[j].GetDecayRateC();
      TP = theDecayRateVector[j].GetTaos();
      
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0){
	G4cout <<"G4RadioactiveDecay::AddDecayRateTable : "
	       << " daughters of ("<< ZP <<", "<<AP<<", "
	       << EP <<") "
	       << " are being calculated. "
	  
	       <<" generation = "
	       << nGeneration << G4endl;
      }
#endif  
      
      theIonTable = (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
      aParentNucleus = theIonTable->GetIon(ZP,AP,EP);
      if (!IsLoaded(*aParentNucleus)){
	aParentNucleus->SetDecayTable(LoadDecayTable(*aParentNucleus));
      }
      aTempDecayTable = aParentNucleus->GetDecayTable();
      
      //
      //
      // Go through the decay table and to combine the same decay channels
      //
      for (i=0; i< 7; i++) brs[i] = 0.0;
      
      G4DecayTable *theDecayTable = new G4DecayTable();
      
      for (i=0; i<aTempDecayTable->entries(); i++) {
	theChannel             = aTempDecayTable->GetDecayChannel(i);
	theNuclearDecayChannel = static_cast<G4NuclearDecayChannel *>(theChannel);
	theDecayMode           = theNuclearDecayChannel->GetDecayMode();
	daughterExcitation = theNuclearDecayChannel->GetDaughterExcitation ();
	theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus () ;
	AD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
	ZD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();  
	G4NuclearLevelManager levelManager = G4NuclearLevelManager (ZD, AD);
	if ( levelManager.NumberOfLevels() ) {
	  const G4NuclearLevel* level = levelManager.NearestLevel (daughterExcitation);

	  if (abs(daughterExcitation - level->Energy()) < levelTolerance) {
	    
	    // Level hafe life is in ns and I want to set the gate as 1 micros
	  
	    if (level->HalfLife() >= 1000.){    
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
	if (brs[i] > 0) {
	  switch ( i ) {
	  case 0:
	    //
	    //
	    // Decay mode is isomeric transition.
	    //
	    
	    theITChannel =  new G4ITDecayChannel
	      (0, aParentNucleus, brs[0]);
	    
	    theDecayTable->Insert(theITChannel);
	    break;
	    
	  case 1:
	    //
	    //
	    // Decay mode is beta-.
	    //
	    theBetaMinusChannel = new G4BetaMinusDecayChannel (0, aParentNucleus,
							       brs[1], 0.*MeV, 0.*MeV, 1, false, NULL);
	    theDecayTable->Insert(theBetaMinusChannel);
	    
	    break;
	    
	  case 2:
	    //
	    //
	    // Decay mode is beta+ + EC.
	    //
	    theBetaPlusChannel = new G4BetaPlusDecayChannel (GetVerboseLevel(), aParentNucleus,
							     brs[2], 0.*MeV, 0.*MeV, 1, false, NULL);
	    theDecayTable->Insert(theBetaPlusChannel);
	    break;		      
	    
	  case 6:
	    //
	    //
	    // Decay mode is alpha.
	    //
	    theAlphaChannel = new G4AlphaDecayChannel (GetVerboseLevel(), aParentNucleus,
						       brs[6], 0.*MeV, 0.*MeV);
	    theDecayTable->Insert(theAlphaChannel);
	    break;
	    
	  default:
	    break;
	  }
	}
      }
	
      // 
      // loop over all braches in theDecayTable
      //
      for ( i=0; i<theDecayTable->entries(); i++){
	theChannel             = theDecayTable->GetDecayChannel(i);
	theNuclearDecayChannel = static_cast<G4NuclearDecayChannel *>(theChannel);
	theBR = theChannel->GetBR();
	theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus();
	
	//
	// now test if the daughterNucleus is a valid one
	//
	if (IsApplicable(*theDaughterNucleus) && theBR ) {
	  A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
	  Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
	  E = ((const G4Ions*)(theDaughterNucleus))->GetExcitationEnergy();
	  
	  TaoPlus = theDaughterNucleus->GetPDGLifeTime();
	  //		cout << TaoPlus <<G4endl;
	  if (TaoPlus > 0.) {
	    // first set the taos, one simply need to add to the parent ones
	    taos.clear();
	    taos = TP;
	    taos.push_back(TaoPlus);
	    // now calculate the coefficiencies
	    //
	    // they are in two parts, first the les than n ones
	    rates.clear();
	    for (k = 0; k < RP.size(); k++){
	      theRate = TP[k]/(TP[k]-TaoPlus) * theBR * RP[k];
	      rates.push_back(theRate);
	    }
	    //
	    // the sencond part: the n:n coefficiency
	    theRate = 0.;
	    for (k = 0; k < RP.size(); k++){
	      theRate -=TaoPlus/(TP[k]-TaoPlus) * theBR * RP[k];
	    }
	    rates.push_back(theRate); 
	    
	    SetDecayRate (Z,A,E,nGeneration,rates,taos);
	    
	    theDecayRateVector.push_back(theDecayRate);
	    
	    nEntry++;
	    
	  }   
	}
      }
	// end of i loop( the branches) 
    }
    //end of for j loop
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
  
////////////////////////////////////////////////////////////////////////////////
//
//
// SetSourceTimeProfile
//
//  read in the source time profile function (histogram)
//
void G4RadioactiveDecay::SetSourceTimeProfile(G4String filename)
{
  G4std::ifstream infile ( filename, G4std::ios::in );
  if ( !infile ) G4Exception ( "Unable to open source data file" );
  
  float bin, flux;
  NSourceBin = -1;
  while (infile >> bin >> flux ) {
    NSourceBin++;
    if (NSourceBin > 99)  G4Exception ( "input source data file too big (>100 rows)" );
    SBin[NSourceBin] = bin * s;
    SProfile[NSourceBin] = flux;
  }
  SetAnalogueMonteCarlo(0);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Source Timeprofile Nbin = " << NSourceBin <<G4endl;}
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
//
// SetDecayBiasProfile
//
// read in the decay bias scheme function (histogram)
//
void G4RadioactiveDecay::SetDecayBias(G4String filename)
{
  G4std::ifstream infile ( filename, G4std::ios::in);
  if ( !infile ) G4Exception ( "Unable to open bias data file" );
  
  float bin, flux;
  NDecayBin = -1;
  while (infile >> bin >> flux ) {
    NDecayBin++;
    if (NDecayBin > 99)  G4Exception ( "input bias data file too big (>100 rows)" );
    DBin[NDecayBin] = bin * s;
    DProfile[NDecayBin] = flux;
  }
  G4int i ;
  for ( i = 1; i<= NDecayBin; i++) DProfile[i] += DProfile[i-1];
  for ( i = 0; i<= NDecayBin; i++) DProfile[i] /= DProfile[NDecayBin];
  SetAnalogueMonteCarlo(0);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<" Decay Bias Profile  Nbin = " << NDecayBin <<G4endl;}
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
//
// DecayIt
//
G4VParticleChange* G4RadioactiveDecay::DecayIt(const G4Track& theTrack, const G4Step& )
{
  //
  // Initialize the G4ParticleChange object. Get the particle details and the
  // decay table.
  //
  fParticleChangeForRadDecay.Initialize(theTrack);
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  G4ParticleDefinition *theParticleDef = theParticle->GetDefinition();

  // First check whether RDM applies to the current logical volume
  //
  if(!G4std::binary_search(ValidVolumes.begin(),
		    ValidVolumes.end(), 
		    theTrack.GetVolume()->GetLogicalVolume()->GetName()))
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	{
	  G4cout <<"G4RadioactiveDecay::DoIt : "
		 << theTrack.GetVolume()->GetLogicalVolume()->GetName()
		 << " is not selected for the RDM"<< G4endl;
	  G4cout << " There are " << ValidVolumes.size() << " volumes" << G4endl;
	  G4cout << " The Valid volumes are " << G4endl;
	  for (G4int i = 0; i< ValidVolumes.size(); i++)
	    G4cout << ValidVolumes[i] << G4endl;
	}
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      //
      //
      // Kill the parent particle.
      //
      fParticleChangeForRadDecay.SetStatusChange( fStopAndKill ) ;
      fParticleChangeForRadDecay.SetLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }
   
  // now check is the particle is valid for RDM
  //
  if (!(IsApplicable(*theParticleDef)))
    { 
      //
      // The particle is not a Ion or outside the nucleuslimits for decay
      //
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	{
	  G4cerr <<"G4RadioactiveDecay::DoIt : "
		 <<theParticleDef->GetParticleName() 
		 << " is not a valid nucleus for the RDM"<< G4endl;
	}
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      //
      //
      // Kill the parent particle.
      //
      fParticleChangeForRadDecay.SetStatusChange( fStopAndKill ) ;
      fParticleChangeForRadDecay.SetLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }
  
  if (!IsLoaded(*theParticleDef))
    {
      theParticleDef->SetDecayTable(LoadDecayTable(*theParticleDef));
    }
  G4DecayTable *theDecayTable = theParticleDef->GetDecayTable();
  
  if  (theDecayTable == NULL || theDecayTable->entries() == 0 )
    {
      //
      //
      // There are no data in the decay table.  Set the particle change parameters
      // to indicate this.
      //
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	{
	  G4cerr <<"G4RadioactiveDecay::DoIt : decay table not defined  for";
	  G4cerr <<theParticleDef->GetParticleName() <<G4endl;
	}
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      //
      //
      // Kill the parent particle.
      //
      fParticleChangeForRadDecay.SetStatusChange( fStopAndKill ) ;
      fParticleChangeForRadDecay.SetLocalEnergyDeposit(0.0);
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
	    //
	    finalGlobalTime += -log( G4UniformRand()) * theParticleDef->GetPDGLifeTime() ;
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
	  G4cout <<"G4RadioactiveDecay::DoIt : Decay vertex :";
	  G4cout <<" Time: " <<finalGlobalTime/ns <<"[ns]";
	  G4cout <<" X:" <<(theTrack.GetPosition()).x() /cm <<"[cm]";
	  G4cout <<" Y:" <<(theTrack.GetPosition()).y() /cm <<"[cm]";
	  G4cout <<" Z:" <<(theTrack.GetPosition()).z() /cm <<"[cm]";
	  G4cout <<G4endl;
	  G4cout <<"G4Decay::DoIt  : decay products in Lab. Frame" <<G4endl;
	  products->DumpInfo();
	}
#endif
	for (index=0; index < numberOfSecondaries; index++) 
	  {
	    G4Track* secondary = new G4Track
	      (products->PopProducts(), finalGlobalTime, currentPosition);
	    secondary->SetGoodForTrackingFlag();
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
	G4std::vector<G4double> PT;
	G4std::vector<G4double> PR;
	G4double taotime;
	G4double decayRate;
	
	G4int i;
	G4int j;
	G4int numberOfSecondaries;
	G4int totalNumberOfSecondaries = 0;
	G4double currentTime;
	G4int ndecaych;
	G4DynamicParticle* asecondaryparticle;
	G4DecayProducts* products = NULL;
	G4std::vector<G4DynamicParticle*> secondaryparticles;
	G4std::vector<G4double> pw;
	pw.clear();
	//now apply the nucleus splitting
	//
	//
	for (G4int n = 0; n < NSplit; n++)
	  {
	    //
	    // Get the decay time following the decay probability function 
	    // suppllied by user
	    //
	    G4double theDecayTime = GetDecayTime();
	    
	    G4int nbin = GetDecayTimeBin(theDecayTime);
	    
	    // claculate the first part of the weight function
	    
	    G4double weight1 =1./DProfile[nbin-1] 
	      *(DBin[nbin]-DBin[nbin-1])
	      /NSplit;
	    if (nbin > 1) {
	       weight1 = 1./(DProfile[nbin]-DProfile[nbin-2])
		 *(DBin[nbin]-DBin[nbin-1])
		 /NSplit;}
	    // it should be calculated in seconds
	    weight1 /= s ;
	    //
	    // loop over all the possible secondaries of the nucleus
	    // the first one is itself.
	    //
	    for ( i = 0; i<theDecayRateVector.size(); i++){
	      PZ = theDecayRateVector[i].GetZ();
	      PA = theDecayRateVector[i].GetA();
	      PE = theDecayRateVector[i].GetE();
	      PT = theDecayRateVector[i].GetTaos();
	      PR = theDecayRateVector[i].GetDecayRateC();
	      
	      // a temprary products buffer and its contents is transfered to 
	      // the products at the end of the loop
	      //
	      G4DecayProducts *tempprods = NULL;
	      
	      // calculate the decay rate of the isotope
	      // one need to fold the the source bias function with the decaytime
	      //
	      decayRate = 0.;
	      for ( j = 0; j < PT.size(); j++){
		taotime = GetTaoTime(theDecayTime,PT[j]);
		decayRate -= PR[j] * taotime;
	      }
	      
	      // decayRatehe radioactivity of isotope (PZ,PA,PE) at the 
	      // time 'theDecayTime'
	      // it will be used to calculate the statistical weight of the 
	      // decay products of this isotope
	      
	      
	      //
	      // now calculate the statistical weight
	      //
	      
	      G4double weight = weight1*decayRate; 
	      // decay the isotope 
	      theIonTable = (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
	      parentNucleus = theIonTable->GetIon(PZ,PA,PE);
	      
	      // decide whther to apply branching ratio bias or not
	      //
	      if (BRBias){
		G4DecayTable *theDecayTable = parentNucleus->GetDecayTable();
		ndecaych = G4int(theDecayTable->entries()*G4UniformRand());
		G4VDecayChannel *theDecayChannel = theDecayTable->GetDecayChannel(ndecaych);
		if (theDecayChannel == NULL)
		  {
		    // Decay channel not found.
#ifdef G4VERBOSE
		    if (GetVerboseLevel()>0)
		      {
			G4cerr <<"G4RadioactiveDecay::DoIt : can not determine decay channel";
			G4cerr <<G4endl;
			theDecayTable ->DumpInfo();
		      }
#endif
		  }
		else
		  {
		    // A decay channel has been identified, so execute the DecayIt.
		    G4double tempmass = parentNucleus->GetPDGMass();      
		    tempprods = theDecayChannel->DecayIt(tempmass);
		    weight *= (theDecayChannel->GetBR())*(theDecayTable->entries());
		  }
		}
	      else {
		tempprods = DoDecay(*parentNucleus);
	      }
	      //
	      // save the secondaries for buffers
	      //
	      numberOfSecondaries = tempprods->entries();
	      currentTime = finalGlobalTime + theDecayTime;
	      for (index=0; index < numberOfSecondaries; index++) 
		{
		  asecondaryparticle = tempprods->PopProducts();
		  if (asecondaryparticle->GetDefinition()->GetBaryonNumber() < 5){
		    pw.push_back(weight);
		    secondaryparticles.push_back(asecondaryparticle);
		  }
		}
	      //
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
		secondaryparticles[index], currentTime, currentPosition);
	    secondary->SetGoodForTrackingFlag(); 	   
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
      fParticleChangeForRadDecay.SetStatusChange( fStopAndKill ) ;
      fParticleChangeForRadDecay.SetLocalEnergyDeposit(energyDeposit);
      // 
      fParticleChangeForRadDecay.SetTimeChange( finalGlobalTime );
      //
      // Reset NumberOfInteractionLengthLeft.
      //
      ClearNumberOfInteractionLengthLeft();
      
      return &fParticleChangeForRadDecay ;
    }
} 

////////////////////////////////////////////////////////////////////////////////
//
//
// DoDecay
//
G4DecayProducts* G4RadioactiveDecay::DoDecay(  G4ParticleDefinition& theParticleDef )
{
  G4DecayProducts *products = 0;
  //
  //
  // follow the decaytable and generate the secondaries...
  // 
  //
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    {
      G4cout<<"Begin of DoDecay..."<<G4endl;
    }
#endif
  G4DecayTable *theDecayTable = theParticleDef.GetDecayTable();
  //
  // Choose a decay channel.
  //
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    {
      G4cout <<"Selecte a channel..."<<G4endl;
    }
#endif
  G4VDecayChannel *theDecayChannel = theDecayTable->SelectADecayChannel();
  if (theDecayChannel == 0)
    {
      // Decay channel not found.
      //
      G4cerr <<"G4RadioactiveDecay::DoIt : can not determine decay channel";
      G4cerr <<G4endl;
      theDecayTable ->DumpInfo();
    }
      else
    {
      //
      // A decay channel has been identified, so execute the DecayIt.
      //
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
	{
	  G4cerr <<"G4RadioactiveDecay::DoIt : selected decay channel  addr:";
	  G4cerr <<theDecayChannel <<G4endl;
	}
#endif
      
      G4double tempmass = theParticleDef.GetPDGMass();
      //
      
      products = theDecayChannel->DecayIt(tempmass);
      
    }
  return products;

}









