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
//  File:   G4RadioactiveDecay.cc                                             //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   9 August 2017                                                     //
//  Description: version the G4RadioactiveDecay process by F. Lei and         //
//               P.R. Truscott with biasing and activation calculations       //
//               removed to a derived class.  It performs alpha, beta,        //
//               electron capture and isomeric transition decays of           //
//               radioactive nuclei.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NuclideTable.hh"
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
#include "G4TritonDecay.hh"
#include "G4ProtonDecay.hh"
#include "G4NeutronDecay.hh"
#include "G4SFDecay.hh"
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
#include "G4Triton.hh"
#include "G4Proton.hh"

#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicException.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhotonEvaporation.hh"
#include "G4HadronicParameters.hh"

#include "G4PhysicsModelCatalog.hh"
#include "G4AutoLock.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace CLHEP;

const G4double G4RadioactiveDecay::levelTolerance = 10.0*CLHEP::eV;
const G4ThreeVector G4RadioactiveDecay::origin(0.,0.,0.);

DecayTableMap* G4RadioactiveDecay::master_dkmap = nullptr;
std::map<G4int, G4String>* G4RadioactiveDecay::theUserRDataFiles = nullptr;
G4String G4RadioactiveDecay::dirPath = "";

namespace
{
  G4Mutex radioactiveDecayMutex = G4MUTEX_INITIALIZER;
}

G4RadioactiveDecay::G4RadioactiveDecay(const G4String& processName, 
                                       const G4double timeThreshold)
 : G4VRestDiscreteProcess(processName, fDecay),
   fThresholdForVeryLongDecayTime( 1.0*CLHEP::year )
{
  if (GetVerboseLevel() > 1) {
    G4cout << "G4RadioactiveDecay constructor: processName = " << processName
           << G4endl;
  }

  SetProcessSubType(fRadioactiveDecay);

  theRadioactiveDecayMessenger = new G4RadioactiveDecayMessenger(this);
  pParticleChange = &fParticleChangeForRadDecay;

  // Check data directory
  if (dirPath.empty()) {
    const char* path_var = G4FindDataDir("G4RADIOACTIVEDATA");
    if (nullptr == path_var) {
      G4Exception("G4RadioactiveDecay()", "HAD_RDM_200", FatalException,
                  "Environment variable G4RADIOACTIVEDATA is not set");
    } else {
      dirPath = path_var;   // convert to string
      std::ostringstream os;
      os << dirPath << "/z1.a3";   // used as a dummy 
      std::ifstream testFile;
      testFile.open(os.str() );
      if ( !testFile.is_open() )
        G4Exception("G4RadioactiveDecay()","HAD_RDM_201",FatalException,
                    "Environment variable G4RADIOACTIVEDATA is set, but does not point to correct directory");
    }
  }
  // Set up photon evaporation for use in G4ITDecay
  photonEvaporation = new G4PhotonEvaporation();
  photonEvaporation->RDMForced(true);
  photonEvaporation->SetICM(true);
  decayIT = new G4ITDecay(photonEvaporation);

  // Instantiate the map of decay tables
  if (nullptr == master_dkmap) {
    master_dkmap = new DecayTableMap();
  }
  if (nullptr == theUserRDataFiles) {
    theUserRDataFiles = new std::map<G4int, G4String>;
  }
 
  // RDM applies to all logical volumes by default
  SelectAllVolumes();
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);

  // The time threshold for radioactive decays can be set in 3 ways:
  // 1. Via C++ interface: G4HadronicParameters::Instance()->SetTimeThresholdForRadioactiveDecay(value)
  // 2. Via the second parameter of the G4RadioactiveDecay constructor
  // 3. Via UI command: /process/had/rdm/thresholdForVeryLongDecayTime value
  // If both 1. and 2. are specified (at the moment when the G4RadioactiveDecay constructor is called),
  // then we take the larger value, to be conservative.
  // If, later on (after invoking the G4RadioactiveDecay constructor) 3. is specified, 
  // then this value is used (and the eventual values 1. and/or 2. are ignored).
  G4double timeThresholdBis = G4HadronicParameters::Instance()->GetTimeThresholdForRadioactiveDecay();
  if ( timeThreshold > 0.0 || timeThresholdBis > 0.0 ) {
    if ( timeThreshold > timeThresholdBis ) timeThresholdBis = timeThreshold;
    fThresholdForVeryLongDecayTime = timeThresholdBis;
  }
}


G4VParticleChange* G4RadioactiveDecay::AtRestDoIt(const G4Track& theTrack,
                                                  const G4Step& theStep)
{
  return DecayIt(theTrack, theStep);
}


G4VParticleChange* G4RadioactiveDecay::PostStepDoIt(const G4Track& theTrack,
                                                    const G4Step& theStep)
{
  return DecayIt(theTrack, theStep);
}


void G4RadioactiveDecay::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The radioactive decay process (G4RadioactiveDecay) handles the\n"
          << "alpha, beta+, beta-, electron capture and isomeric transition\n"
          << "decays of nuclei (G4GenericIon) with masses A > 4.\n"
          << "The required half-lives and decay schemes are retrieved from\n"
          << "the RadioactiveDecay database which was derived from ENSDF.\n";
}


G4RadioactiveDecay::~G4RadioactiveDecay()
{
  delete theRadioactiveDecayMessenger;
  delete photonEvaporation;
  delete decayIT;
  if (nullptr != master_dkmap) {
    G4AutoLock lk(&radioactiveDecayMutex);
    if (nullptr != master_dkmap) {
      for (auto const & i : *master_dkmap) {
	delete i.second;
      }
      master_dkmap->clear();
      delete master_dkmap;
      master_dkmap = nullptr;
    }
    delete theUserRDataFiles;
    theUserRDataFiles = nullptr;
    lk.unlock();
  }
}


G4bool G4RadioactiveDecay::IsApplicable(const G4ParticleDefinition& aParticle)
{
  const G4String& pname = aParticle.GetParticleName();
  if (pname == "GenericIon" || pname == "triton") { return true; }
  // All particles other than G4Ions, are rejected by default
  const G4Ions* p = dynamic_cast<const G4Ions*>(&aParticle);
  if (nullptr == p) { return false; }

  // excited isomere may decay via gamma evaporation
  if (p->GetExcitationEnergy() > 0.0) { return true; }

  // Check on life time 
  G4double lifeTime = p->GetPDGLifeTime();
  if (lifeTime < 0.0 || lifeTime > fThresholdForVeryLongDecayTime) {
    return false;
  }

  // Determine whether the nuclide falls into the correct A and Z range
  G4int A = p->GetAtomicMass();
  G4int Z = p->GetAtomicNumber();

  if (A > theNucleusLimits.GetAMax() || A < theNucleusLimits.GetAMin() ||
      Z > theNucleusLimits.GetZMax() || Z < theNucleusLimits.GetZMin()) {
    return false;
  }

  return true;
}


G4DecayTable* G4RadioactiveDecay::GetDecayTable(const G4ParticleDefinition* aNucleus)
{
  G4String key = aNucleus->GetParticleName();
  auto ptr = master_dkmap->find(key);

  G4DecayTable* theDecayTable = nullptr;
  if ( ptr == master_dkmap->end() ) {
    // Load new file if table not there
    const G4Ions* ion = dynamic_cast<const G4Ions*>(aNucleus);
    if (nullptr != ion) {
      theDecayTable = LoadDecayTable(ion);
    }
  } else {
    theDecayTable = ptr->second;
  }
  return theDecayTable;
}


void G4RadioactiveDecay::SelectAVolume(const G4String& aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* volume = nullptr;
  volume = theLogicalVolumes->GetVolume(aVolume);
  if (volume != nullptr)
  {
    ValidVolumes.push_back(aVolume);
    std::sort(ValidVolumes.begin(), ValidVolumes.end());
    // sort need for performing binary_search

    if (GetVerboseLevel() > 0)
      G4cout << " Radioactive decay applied to " << aVolume << G4endl;
  }
  else
  {
    G4ExceptionDescription ed;
    ed << aVolume << " is not a valid logical volume name."
       << " Decay not activated for it."
       << G4endl;
    G4Exception("G4RadioactiveDecay::SelectAVolume()", "HAD_RDM_300",
                JustWarning, ed);
  }
}


void G4RadioactiveDecay::DeselectAVolume(const G4String& aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* volume = nullptr;
  volume = theLogicalVolumes->GetVolume(aVolume);
  if (volume != nullptr)
  {
    auto location= std::find(ValidVolumes.cbegin(),ValidVolumes.cend(),aVolume);
    if (location != ValidVolumes.cend() )
    {
      ValidVolumes.erase(location);
      std::sort(ValidVolumes.begin(), ValidVolumes.end());
      isAllVolumesMode = false;
      if (GetVerboseLevel() > 0)
        G4cout << " G4RadioactiveDecay::DeselectAVolume: " << aVolume
               << " is removed from list " << G4endl;
    }
    else
    {
      G4ExceptionDescription ed;
      ed << aVolume << " is not in the list.  No action taken." << G4endl;
      G4Exception("G4RadioactiveDecay::DeselectAVolume()", "HAD_RDM_300",
                  JustWarning, ed);
    }
  }
  else
  {
    G4ExceptionDescription ed;
    ed << aVolume << " is not a valid logical volume name.  No action taken." 
       << G4endl;
    G4Exception("G4RadioactiveDecay::DeselectAVolume()", "HAD_RDM_300",
                JustWarning, ed);
  }
}


void G4RadioactiveDecay::SelectAllVolumes() 
{
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* volume = nullptr;
  ValidVolumes.clear();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    G4cout << " RDM Applies to all Volumes"  << G4endl;
#endif
  for (std::size_t i = 0; i < theLogicalVolumes->size(); ++i){
    volume = (*theLogicalVolumes)[i];
    ValidVolumes.push_back(volume->GetName());    
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
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
  if (GetVerboseLevel() > 1) G4cout << "RDM removed from all volumes" << G4endl; 
#endif
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanLifeTime (required by the base class)                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecay::GetMeanLifeTime(const G4Track& theTrack,
                                             G4ForceCondition*)
{
  G4double meanlife = DBL_MAX;
  const G4ParticleDefinition* theParticleDef = theTrack.GetParticleDefinition();
  if (!IsApplicable(*theParticleDef)) { return meanlife; }
  G4double theLife = theParticleDef->GetPDGLifeTime();
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "G4RadioactiveDecay::GetMeanLifeTime() for " 
	   << theParticleDef->GetParticleName() << G4endl;
    G4cout << "KineticEnergy(GeV)=" << theTrack.GetKineticEnergy()/CLHEP::GeV
           << " Mass(GeV)=" << theParticleDef->GetPDGMass()/CLHEP::GeV
           << " LifeTime(ns)=" << theLife/CLHEP::ns << G4endl;
  }
#endif
  if (theLife >= 0.0 && theLife <= fThresholdForVeryLongDecayTime) {
    meanlife = theLife;
  } 

  if (meanlife == DBL_MAX) {
    const G4Ions* ion = dynamic_cast<const G4Ions*>(theParticleDef);
    if (nullptr != ion && ion->GetExcitationEnergy() > 0.0) {
      meanlife = 0.0;
    }
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2)
    G4cout << "G4RadioactiveDecay::GetMeanLifeTime: " 
	   << meanlife/CLHEP::s << " second " << G4endl;
#endif

  return meanlife;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanFreePath for decay in flight                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4RadioactiveDecay::GetMeanFreePath(const G4Track& aTrack, G4double,
                                             G4ForceCondition* fc)
{
  G4double res = DBL_MAX;
  G4double lifeTime = GetMeanLifeTime(aTrack, fc);
  if (lifeTime > 0.0 && lifeTime < DBL_MAX) {
    auto dParticle = aTrack.GetDynamicParticle();
    res = lifeTime*dParticle->GetTotalEnergy()*aTrack.GetVelocity()/dParticle->GetMass(); 
  } else {
    res = lifeTime;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "G4RadioactiveDecay::GetMeanFreePath() for "
	   << aTrack.GetDefinition()->GetParticleName() << G4endl;
    G4cout << "  kinEnergy(GeV)=" << aTrack.GetKineticEnergy()/CLHEP::GeV
           << " lifeTime(ns)=" << lifeTime
           << " mean free path(cm)=" << res/CLHEP::cm << G4endl;
  }
#endif
  return res;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  BuildPhysicsTable - initialization of atomic de-excitation                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void G4RadioactiveDecay::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if (isInitialised) { return; }
  isInitialised = true;
  if (G4HadronicParameters::Instance()->GetVerboseLevel() > 0  &&
      G4Threading::IsMasterThread() && "GenericIon" == p.GetParticleName()) {
    StreamInfo(G4cout, "\n");
  }
  photonEvaporation->Initialise();
  photonEvaporation->RDMForced(true);
  photonEvaporation->SetICM(true);
  decayIT->SetARM(applyARM);

  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  StreamInfo - stream out parameters                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void
G4RadioactiveDecay::StreamInfo(std::ostream& os, const G4String& endline)
{
  G4DeexPrecoParameters* deex =
    G4NuclearLevelData::GetInstance()->GetParameters();
  G4EmParameters* emparam = G4EmParameters::Instance();
  G4double minMeanLife = G4NuclideTable::GetInstance()->GetMeanLifeThreshold();

  G4long prec = os.precision(5);
  os << "======================================================================"
     << endline;
  os << "======          Radioactive Decay Physics Parameters           ======="
     << endline;
  os << "======================================================================"
     << endline;
  os << "min MeanLife (from G4NuclideTable)                "
     << G4BestUnit(minMeanLife, "Time") << endline;     
  os << "Max life time (from G4DeexPrecoParameters)        "
     << G4BestUnit(deex->GetMaxLifeTime(), "Time") << endline;
  os << "Internal e- conversion flag                       "
     << deex->GetInternalConversionFlag() << endline;
  os << "Stored internal conversion coefficients           "
     << deex->StoreICLevelData() << endline;
  os << "Enabled atomic relaxation mode                    "
     << applyARM << endline;
  os << "Enable correlated gamma emission                  "
     << deex->CorrelatedGamma() << endline;
  os << "Max 2J for sampling of angular correlations       "
     << deex->GetTwoJMAX() << endline;
  os << "Atomic de-excitation enabled                      "
     << emparam->Fluo() << endline;
  os << "Auger electron emission enabled                   "
     << emparam->Auger() << endline;
  os << "Check EM cuts disabled for atomic de-excitation   "
     << emparam->DeexcitationIgnoreCut() << endline;
  os << "Use Bearden atomic level energies                 "
     << emparam->BeardenFluoDir() << endline;
  os << "Use ANSTO fluorescence model                      "
     << emparam->ANSTOFluoDir() << endline;
  os << "Threshold for very long decay time at rest        "
     << G4BestUnit(fThresholdForVeryLongDecayTime, "Time") << endline;
  os << "======================================================================"
     << G4endl;
  os.precision(prec);
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  LoadDecayTable loads the decay scheme from the RadioactiveDecay database  //
//  for the parent nucleus.                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4DecayTable* G4RadioactiveDecay::LoadDecayTable(const G4Ions* theIon)
{
  G4AutoLock lk(&radioactiveDecayMutex);
  const G4String key = theIon->GetParticleName();
  auto dtptr = master_dkmap->find(key);
  if (dtptr != master_dkmap->end()) {
    lk.unlock();
    return dtptr->second;
  }

  // Generate input data file name using Z and A of the parent nucleus
  // file containing radioactive decay data.
  G4int A = theIon->GetAtomicMass();
  G4int Z = theIon->GetAtomicNumber();

  //G4cout << "LoadDecayTable for " << key << " Z=" << Z << " A=" << A << G4endl;

  G4double levelEnergy = theIon->GetExcitationEnergy();
  G4Ions::G4FloatLevelBase floatingLevel = theIon->GetFloatLevelBase();

  //Check if data have been provided by the user
  G4String file;
  G4int ke = 1000*A + Z;
  auto ptr = theUserRDataFiles->find(ke);
  if (ptr != theUserRDataFiles->end()) {
    file = ptr->second;
  } else {
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
    const G4int nMode = G4RadioactiveDecayModeSize;
    G4double modeTotalBR[nMode] = {0.0};
    G4double modeSumBR[nMode];
    for (G4int i = 0; i < nMode; ++i) {
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
    /* Loop checking, 01.09.2015, D.Wright */
    while (!complete && !DecaySchemeFile.getline(inputChars, 120).eof()) {
      loop++;
      if (loop > 100000) {
        G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_100",
                    JustWarning, "While loop count exceeded");
        break;
      }
 
      inputLine = inputChars;
      G4StrUtil::rstrip(inputLine);
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
		  G4ITDecay* anITChannel = new G4ITDecay(theIon, decayModeTotal, 0.0, 0.0);
		  theDecayTable->Insert(anITChannel);
                }
                break;
              case BetaMinus:
                modeTotalBR[BetaMinus] = decayModeTotal; break;
              case BetaPlus:
                modeTotalBR[BetaPlus] = decayModeTotal; break;
              case KshellEC:
                modeTotalBR[KshellEC] = decayModeTotal; break;
              case LshellEC:
                modeTotalBR[LshellEC] = decayModeTotal; break;
              case MshellEC:
                modeTotalBR[MshellEC] = decayModeTotal; break;
              case NshellEC:
                modeTotalBR[NshellEC] = decayModeTotal; break;
              case Alpha:
                modeTotalBR[Alpha] = decayModeTotal; break;
              case Proton:
                modeTotalBR[Proton] = decayModeTotal; break;
              case Neutron:
                modeTotalBR[Neutron] = decayModeTotal; break;
              case SpFission:
                modeTotalBR[SpFission] = decayModeTotal; break;
              case BDProton:
                /* Not yet implemented */  break;
              case BDNeutron:
                /* Not yet implemented */  break;
              case Beta2Minus:
                /* Not yet implemented */  break;
              case Beta2Plus:
                /* Not yet implemented */  break;
              case Proton2:
                /* Not yet implemented */  break;
              case Neutron2:
                /* Not yet implemented */  break;
              case Triton:
                modeTotalBR[Triton] = decayModeTotal; break;
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
            b /= 100.;
            daughterFloatLevel = G4Ions::FloatLevelBase(daughterFloatFlag.back());

            switch (theDecayMode) {
              case BetaMinus:
              {
                G4BetaMinusDecay* aBetaMinusChannel =
                  new G4BetaMinusDecay(theIon, b, c*MeV, a*MeV,
                                       daughterFloatLevel, betaType);
		//aBetaMinusChannel->DumpNuclearInfo();
                theDecayTable->Insert(aBetaMinusChannel);
		modeSumBR[BetaMinus] += b;
              }
              break;

              case BetaPlus:
              {
                G4BetaPlusDecay* aBetaPlusChannel =
                  new G4BetaPlusDecay(theIon, b, c*MeV, a*MeV,
                                      daughterFloatLevel, betaType);
		//aBetaPlusChannel->DumpNuclearInfo();
                theDecayTable->Insert(aBetaPlusChannel);
		modeSumBR[BetaPlus] += b;
              }
              break;

              case KshellEC:  // K-shell electron capture
              {
                G4ECDecay* aKECChannel =
                  new G4ECDecay(theIon, b, c*MeV, a*MeV,
                                daughterFloatLevel, KshellEC);
		//aKECChannel->DumpNuclearInfo();
                aKECChannel->SetARM(applyARM);
                theDecayTable->Insert(aKECChannel);
                modeSumBR[KshellEC] += b;
              }
              break;

              case LshellEC:  // L-shell electron capture
              {
                G4ECDecay* aLECChannel =
                  new G4ECDecay(theIon, b, c*MeV, a*MeV,
                                daughterFloatLevel, LshellEC);
//              aLECChannel->DumpNuclearInfo();
                aLECChannel->SetARM(applyARM);
                theDecayTable->Insert(aLECChannel);
                modeSumBR[LshellEC] += b;
              }
              break;

              case MshellEC:  // M-shell electron capture
              {
                G4ECDecay* aMECChannel =
                  new G4ECDecay(theIon, b, c*MeV, a*MeV,
                                daughterFloatLevel, MshellEC);
//              aMECChannel->DumpNuclearInfo();
                aMECChannel->SetARM(applyARM);
                theDecayTable->Insert(aMECChannel);
                modeSumBR[MshellEC] += b;
              }
              break;

              case NshellEC:  // N-shell electron capture
              {
                G4ECDecay* aNECChannel =
                  new G4ECDecay(theIon, b, c*MeV, a*MeV,
                                daughterFloatLevel, NshellEC);
//              aNECChannel->DumpNuclearInfo();
                aNECChannel->SetARM(applyARM);
                theDecayTable->Insert(aNECChannel);
                modeSumBR[NshellEC] += b;
              }
              break;

              case Alpha:
              {
                G4AlphaDecay* anAlphaChannel =
                  new G4AlphaDecay(theIon, b, c*MeV, a*MeV,
                                   daughterFloatLevel);
//              anAlphaChannel->DumpNuclearInfo();
                theDecayTable->Insert(anAlphaChannel);
                modeSumBR[Alpha] += b;
              }
              break;

	      case Proton:
              {
                G4ProtonDecay* aProtonChannel =
                  new G4ProtonDecay(theIon, b, c*MeV, a*MeV,
                                    daughterFloatLevel);
//              aProtonChannel->DumpNuclearInfo();
                theDecayTable->Insert(aProtonChannel);
                modeSumBR[Proton] += b;
              }
              break;

              case Neutron:
              {
                G4NeutronDecay* aNeutronChannel =
                  new G4NeutronDecay(theIon, b, c*MeV, a*MeV,
                                     daughterFloatLevel);
//              aNeutronChannel->DumpNuclearInfo();
                theDecayTable->Insert(aNeutronChannel);
                modeSumBR[Neutron] += b;
              }
              break;

              case SpFission:
              {
                G4SFDecay* aSpontFissChannel =
                  new G4SFDecay(theIon, b, c*MeV, a*MeV, daughterFloatLevel);
                theDecayTable->Insert(aSpontFissChannel);
                modeSumBR[SpFission] += b;
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

              case Triton:
                {
		  G4TritonDecay* aTritonChannel =
                    new G4TritonDecay(theIon, b, c*MeV, a*MeV, daughterFloatLevel);
		  theDecayTable->Insert(aTritonChannel);
		  modeSumBR[Triton] += b;
                }
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
    for (G4int i = 0; i < theDecayTable->entries(); ++i) {
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
    G4ITDecay* anITChannel = new G4ITDecay(theIon, 1.0, 0.0, 0.0);
    theDecayTable->Insert(anITChannel);
  }

  if (GetVerboseLevel() > 1) {
    theDecayTable->DumpInfo();
  }

  // store in master library 
  (*master_dkmap)[theIon->GetParticleName()] = theDecayTable;
  lk.unlock();
  return theDecayTable;
}

void G4RadioactiveDecay::AddUserDecayDataFile(G4int Z, G4int A,
                                              const G4String& filename)
{
  if (Z < 1 || A < 2) G4cout << "Z and A not valid!" << G4endl;

  std::ifstream DecaySchemeFile(filename);
  if (DecaySchemeFile) {
    G4int ID_ion = A*1000 + Z;
    (*theUserRDataFiles)[ID_ion] = filename;
  } else {
    G4ExceptionDescription ed;
    ed << filename << " does not exist! " << G4endl;
    G4Exception("G4RadioactiveDecay::AddUserDecayDataFile()", "HAD_RDM_001",
                FatalException, ed);
  }
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
      if (GetVerboseLevel()>1) {
        G4cout <<"G4RadioactiveDecay::DecayIt : "
               << theTrack.GetVolume()->GetLogicalVolume()->GetName()
               << " is not selected for the RDM"<< G4endl;
        G4cout << " There are " << ValidVolumes.size() << " volumes" << G4endl;
        G4cout << " The Valid volumes are: ";
	for (auto const & vol : ValidVolumes) { G4cout << vol << " " << G4endl; }
        G4cout << G4endl;
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
  G4DecayTable* theDecayTable = GetDecayTable(theParticleDef);
  if ( theDecayTable == nullptr || theDecayTable->entries() == 0) { 
    // Particle is not an ion or is outside the nucleuslimits for decay
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "G4RadioactiveDecay::DecayIt : "
             << theParticleDef->GetParticleName() 
             << " is outside (Z,A) limits set for the decay or has no decays." 
             << G4endl;
    }
#endif
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;
  }
  //G4cout << "DecayIt for " << theParticleDef->GetParticleName() 
  //	 << " isAllVolumesMode:" << isAllVolumesMode
  //	 << " decayTable=" << theDecayTable << G4endl;

  // Data found. Decay nucleus without variance reduction. 
  DecayAnalog(theTrack, theDecayTable);
  return &fParticleChangeForRadDecay;
} 

void G4RadioactiveDecay::DecayAnalog(const G4Track& theTrack,
                                     G4DecayTable* decayTable)
{
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();
  //G4cout << "DecayIt for " << theParticleDef->GetParticleName() << G4endl;
  G4DecayProducts* products = DoDecay(*theParticleDef, decayTable);

  // Check if the product is the same as input and kill the track if
  // necessary to prevent infinite loop (11/05/10, F.Lei)
  if (nullptr == products || products->entries() == 1) {
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill);
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    delete products;
    return;
  }

  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = theTrack.GetGlobalTime();
  G4double finalLocalTime = theTrack.GetLocalTime();

  // Get parent particle information and boost the decay products to the
  // laboratory frame

  // ParentEnergy used for the boost should be the total energy of the nucleus
  // of the parent ion without the energy of the shell electrons
  // (correction for bug 1359 by L. Desorgher)
  G4double ParentEnergy = theParticle->GetKineticEnergy()
                        + theParticle->GetParticleDefinition()->GetPDGMass();
  G4ThreeVector ParentDirection(theParticle->GetMomentumDirection());

  if (theTrack.GetTrackStatus() == fStopButAlive) {
    // this condition seems to be always True, further investigation is needed (L.Desorgher)

    // The particle is decayed at rest
    // Since the time is for the particle at rest, need to add additional time
    // lapsed between particle coming to rest and the actual decay.  This time
    // is sampled with the mean-life of the particle.  Need to protect the case 
    // PDGTime < 0.  (F.Lei 11/05/10)
    G4double temptime = -std::log(G4UniformRand() ) *
                        theParticleDef->GetPDGLifeTime();
    if (temptime < 0.) temptime = 0.;
    finalGlobalTime += temptime;
    finalLocalTime += temptime;
    energyDeposit += theParticle->GetKineticEnergy();
    
    // Kill the parent particle, and ignore its decay, if it decays later than the
    // threshold fThresholdForVeryLongDecayTime (whose default value corresponds
    // to more than twice the age of the universe).
    // This kind of cut has been introduced (in April 2021) in order to avoid to
    // account energy depositions happening after many billions of years in
    // ordinary materials used in calorimetry, in particular Tungsten and Lead
    // (via their natural unstable, but very long lived, isotopes, such as
    // W183, W180 and Pb204).
    // Note that the cut is not on the average, mean lifetime, but on the actual
    // sampled global decay time.
    if ( finalGlobalTime > fThresholdForVeryLongDecayTime ) {
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
      fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      delete products;
      return;
    }     
  }
  products->Boost(ParentEnergy, ParentDirection);

  // Add products in theParticleChangeForRadDecay.
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForRadDecay.SetNumberOfSecondaries(numberOfSecondaries);

  if (GetVerboseLevel() > 1) {
    G4cout << "G4RadioactiveDecay::DecayAnalog: Decay vertex :";
    G4cout << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cout << " X:" << (theTrack.GetPosition()).x() /cm << "[cm]";
    G4cout << " Y:" << (theTrack.GetPosition()).y() /cm << "[cm]";
    G4cout << " Z:" << (theTrack.GetPosition()).z() /cm << "[cm]";
    G4cout << G4endl;
    G4cout << "G4Decay::DecayIt : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
    products->IsChecked();
  }

  const G4int modelID_forIT = G4PhysicsModelCatalog::GetModelID( "model_RDM_IT" );
  G4int modelID = modelID_forIT + 10*theRadDecayMode;
  const G4int modelID_forAtomicRelaxation =
    G4PhysicsModelCatalog::GetModelID( "model_RDM_AtomicRelaxation" );
  for ( G4int index = 0; index < numberOfSecondaries; ++index ) {
    G4Track* secondary = new G4Track( products->PopProducts(), finalGlobalTime,
                                      theTrack.GetPosition() );
    secondary->SetWeight( theTrack.GetWeight() );
    secondary->SetCreatorModelID( modelID );
    // Change for atomics relaxation
    if ( theRadDecayMode == IT  &&  index > 0 ) {
      if ( index == numberOfSecondaries-1 ) {
	secondary->SetCreatorModelID( modelID_forIT );
      } else {
	secondary->SetCreatorModelID( modelID_forAtomicRelaxation) ;
      }
    } else if ( theRadDecayMode >= KshellEC  &&  theRadDecayMode <= NshellEC  &&
		index < numberOfSecondaries-1 ) {
      secondary->SetCreatorModelID( modelID_forAtomicRelaxation );
    }
    secondary->SetGoodForTrackingFlag();
    secondary->SetTouchableHandle( theTrack.GetTouchableHandle() );
    fParticleChangeForRadDecay.AddSecondary( secondary );
  }

  delete products;

  // Kill the parent particle
  fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
  fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(energyDeposit);
  fParticleChangeForRadDecay.ProposeLocalTime(finalLocalTime);

  // Reset NumberOfInteractionLengthLeft.
  ClearNumberOfInteractionLengthLeft();
}


G4DecayProducts*
G4RadioactiveDecay::DoDecay(const G4ParticleDefinition& theParticleDef,
                            G4DecayTable* theDecayTable)
{
  G4DecayProducts* products = nullptr;

  // Choose a decay channel.
  // G4DecayTable::SelectADecayChannel checks to see if sum of daughter masses
  // exceeds parent mass. Pass it the parent mass + maximum Q value to account
  // for difference in mass defect.
  G4double parentPlusQ = theParticleDef.GetPDGMass() + 30.*MeV;
  G4VDecayChannel* theDecayChannel = theDecayTable->SelectADecayChannel(parentPlusQ);

  if (theDecayChannel == nullptr) {
    // Decay channel not found.
    G4ExceptionDescription ed;
    ed << " Cannot determine decay channel for " << theParticleDef.GetParticleName() << G4endl;
    G4Exception("G4RadioactiveDecay::DoDecay", "HAD_RDM_013",
                FatalException, ed);
  } else {
    // A decay channel has been identified, so execute the DecayIt.
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "G4RadioactiveDecay::DoIt : selected decay channel addr: "
             << theDecayChannel << G4endl;
    }
#endif
    theRadDecayMode = static_cast<G4NuclearDecay*>(theDecayChannel)->GetDecayMode();

    // for IT decay use local G4ITDecay class
    if (theRadDecayMode == IT) {
      decayIT->SetupDecay(&theParticleDef);
      products = decayIT->DecayIt(0.0);
    } else { 
      // for others decayes use shared  class
      products = theDecayChannel->DecayIt(theParticleDef.GetPDGMass());
    }

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
  if (GetVerboseLevel() > 1) G4cout << "Begin of CollimateDecay..." << G4endl;
#endif

  // Particles suitable for directional biasing (for if-blocks below)
  static const G4ParticleDefinition* electron = G4Electron::Definition();
  static const G4ParticleDefinition* positron = G4Positron::Definition();
  static const G4ParticleDefinition* neutron  = G4Neutron::Definition();
  static const G4ParticleDefinition* gamma    = G4Gamma::Definition();
  static const G4ParticleDefinition* alpha    = G4Alpha::Definition();
  static const G4ParticleDefinition* triton   = G4Triton::Definition();
  static const G4ParticleDefinition* proton   = G4Proton::Definition();

  G4ThreeVector newDirection;		// Re-use to avoid memory churn
  for (G4int i=0; i<products->entries(); ++i) {
    G4DynamicParticle* daughter = (*products)[i];
    const G4ParticleDefinition* daughterType =
                                  daughter->GetParticleDefinition();
    if (daughterType == electron || daughterType == positron ||
	daughterType == neutron || daughterType == gamma ||
	daughterType == alpha || daughterType == triton || daughterType == proton) {
      CollimateDecayProduct(daughter);
    }
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

