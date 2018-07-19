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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4RadioactiveDecayBase.hh                                         //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   9 August 2017                                                     //
//  Description: version the G4RadioactiveDecay process by F. Lei and         //
//               P.R. Truscott with biasing and activation calculations       //
//               removed to a derived class.  It performs alpha, beta,        //
//               electron capture and isomeric transition decays of           //
//               radioactive nuclei.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4RadioactiveDecayBase.hh"
#include "G4RadioactiveDecayBaseMessenger.hh"

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
#include "G4HadronicProcessStore.hh"
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

const G4double G4RadioactiveDecayBase::levelTolerance = 10.0*eV;
const G4ThreeVector G4RadioactiveDecayBase::origin(0.,0.,0.);

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
G4Mutex G4RadioactiveDecayBase::radioactiveDecayMutex = G4MUTEX_INITIALIZER;
DecayTableMap* G4RadioactiveDecayBase::master_dkmap = 0;
#endif

G4RadioactiveDecayBase::G4RadioactiveDecayBase(const G4String& processName)
 : G4VRestDiscreteProcess(processName, fDecay), isInitialised(false),
   forceDecayDirection(0.,0.,0.), forceDecayHalfAngle(0.*deg), dirPath(""),
   verboseLevel(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4RadioactiveDecayBase constructor: processName = " << processName
           << G4endl;
  }
#endif

  SetProcessSubType(fRadioactiveDecay);

  theRadioactiveDecayBaseMessenger = new G4RadioactiveDecayBaseMessenger(this);
  pParticleChange = &fParticleChangeForRadDecay;

  // Set up photon evaporation for use in G4ITDecay
  photonEvaporation = new G4PhotonEvaporation();
  photonEvaporation->RDMForced(true);
  photonEvaporation->SetICM(true);

  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(true);

  // Check data directory
  char* path_var = getenv("G4RADIOACTIVEDATA");
  if (!path_var) {
    G4Exception("G4RadioactiveDecay()","HAD_RDM_200",FatalException,
                "Environment variable G4RADIOACTIVEDATA is not set");
  } else {
    dirPath = path_var;   // convert to string
    std::ostringstream os;
    os << dirPath << "/z1.a3";   // used as a dummy 
    std::ifstream testFile;
    testFile.open(os.str() );
    if (!testFile.is_open() )
      G4Exception("G4RadioactiveDecay()","HAD_RDM_201",FatalException,
                  "Environment variable G4RADIOACTIVEDATA is set, but does not point to correct directory");
  }

  // Reset the list of user defined data files
  theUserRadioactiveDataFiles.clear();

  // Instantiate the map of decay tables
#ifdef G4MULTITHREADED
  G4AutoLock lk(&G4RadioactiveDecayBase::radioactiveDecayMutex);
  if(!master_dkmap) master_dkmap = new DecayTableMap;
#endif
  dkmap = new DecayTableMap;

  // Apply default values
  applyARM = true;
  applyICM = true;  // Always on; keep only for backward compatibility
 
  // RDM applies to all logical volumes by default
  isAllVolumesMode = true;
  SelectAllVolumes();
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}

void G4RadioactiveDecayBase::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The radioactive decay process (G4RadioactiveDecay) handles the\n"
          << "alpha, beta+, beta-, electron capture and isomeric transition\n"
          << "decays of nuclei (G4GenericIon) with masses A > 4.\n"
          << "The required half-lives and decay schemes are retrieved from\n"
          << "the RadioactiveDecay database which was derived from ENSDF.\n";
}


G4RadioactiveDecayBase::~G4RadioactiveDecayBase()
{
  delete theRadioactiveDecayBaseMessenger;
  delete photonEvaporation;
  for (DecayTableMap::iterator i = dkmap->begin(); i != dkmap->end(); i++) {
    delete i->second;
  }
  dkmap->clear();
  delete dkmap;
}


G4bool G4RadioactiveDecayBase::IsApplicable(const G4ParticleDefinition& aParticle)
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

G4DecayTable* G4RadioactiveDecayBase::GetDecayTable(const G4ParticleDefinition* aNucleus)
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


void G4RadioactiveDecayBase::SelectAVolume(const G4String aVolume)
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


void G4RadioactiveDecayBase::DeselectAVolume(const G4String aVolume)
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


void G4RadioactiveDecayBase::SelectAllVolumes() 
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


void G4RadioactiveDecayBase::DeselectAllVolumes() 
{
  ValidVolumes.clear();
  isAllVolumesMode=false;
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "RDM removed from all volumes" << G4endl; 
#endif
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanLifeTime (required by the base class)                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecayBase::GetMeanLifeTime(const G4Track& theTrack,
                                             G4ForceCondition*)
{
  // For variance reduction the time is set to 0 so as to force the particle
  // to decay immediately.
  // In analogueMC mode it returns the particle's mean-life.

  G4double meanlife = 0.;
//  if (AnalogueMC) {
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
//  }
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

G4double G4RadioactiveDecayBase::GetMeanFreePath (const G4Track& aTrack, G4double,
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

void G4RadioactiveDecayBase::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (!isInitialised) {
    isInitialised = true;
  }

  G4HadronicProcessStore::
   Instance()->RegisterParticleForExtraProcess(this,G4GenericIon::GenericIon());
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  LoadDecayTable loads the decay scheme from the RadioactiveDecay database  // 
//  for the parent nucleus.                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4DecayTable*
G4RadioactiveDecayBase::LoadDecayTable(const G4ParticleDefinition& theParentNucleus)
{
  // Generate input data file name using Z and A of the parent nucleus
  // file containing radioactive decay data.
  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();

  G4double levelEnergy = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4Ions::G4FloatLevelBase floatingLevel =
    ((const G4Ions*)(&theParentNucleus))->GetFloatLevelBase();

#ifdef G4MULTITHREADED
  G4AutoLock lk(&G4RadioactiveDecayBase::radioactiveDecayMutex);

  G4String key = theParentNucleus.GetParticleName();
  DecayTableMap::iterator master_table_ptr = master_dkmap->find(key);

  if (master_table_ptr != master_dkmap->end() ) {   // If table is there              
    return master_table_ptr->second;
  }
#endif

  //Check if data have been provided by the user
  G4String file = theUserRadioactiveDataFiles[1000*A+Z];

  if (file == "") {
/*
    if (!getenv("G4RADIOACTIVEDATA") ) {
      G4cout << "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files."
             << G4endl;
      throw G4HadronicException(__FILE__, __LINE__, " Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files.");
    }
    G4String dirName = getenv("G4RADIOACTIVEDATA");
*/
    std::ostringstream os;
    os << dirPath << "/z" << Z << ".a" << A << '\0';
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
//                anITChannel->SetHLThreshold(halflifethreshold);
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
//                aBetaMinusChannel->SetHLThreshold(halflifethreshold);
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
//                aBetaPlusChannel->SetHLThreshold(halflifethreshold);
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
//                aKECChannel->SetHLThreshold(halflifethreshold);
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
//                aLECChannel->SetHLThreshold(halflifethreshold);
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
//                aMECChannel->SetHLThreshold(halflifethreshold);
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
//                anAlphaChannel->SetHLThreshold(halflifethreshold);
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
//                aProtonChannel->SetHLThreshold(halflifethreshold);
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
//                aNeutronChannel->SetHLThreshold(halflifethreshold);
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
//    anITChannel->SetHLThreshold(halflifethreshold);
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
G4RadioactiveDecayBase::AddUserDecayDataFile(G4int Z, G4int A, G4String filename)
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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DecayIt                                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4VParticleChange*
G4RadioactiveDecayBase::DecayIt(const G4Track& theTrack, const G4Step&)
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

    G4DecayProducts* products = DoDecay(*theParticleDef);

    // If the product is the same as the input kill the track if
    // necessary to prevent infinite loop (11/05/10, F.Lei)
    if (products->entries() == 1) {
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill);
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }

    // Get parent particle information and boost the decay products to the
    // laboratory frame based on this information.

    // The Parent Energy used for the boost should be the total energy of
    // the nucleus of the parent ion without the energy of the shell electrons
    // (correction for bug 1359 by L. Desorgher)
    G4double ParentEnergy = theParticle->GetKineticEnergy()
                          + theParticle->GetParticleDefinition()->GetPDGMass();
    G4ThreeVector ParentDirection(theParticle->GetMomentumDirection());

    if (theTrack.GetTrackStatus() == fStopButAlive) {
      // This condition seems to be always True, further investigation is needed
      // (L.Desorgher)
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
G4RadioactiveDecayBase::DoDecay(const G4ParticleDefinition& theParticleDef)
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

void G4RadioactiveDecayBase::CollimateDecay(G4DecayProducts* products) {
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

void G4RadioactiveDecayBase::CollimateDecayProduct(G4DynamicParticle* daughter) {
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

G4ThreeVector G4RadioactiveDecayBase::ChooseCollimationDirection() const {
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

