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
//  File:   G4Radioactivation.cc                                              //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   29 August 2017                                                    //
//  Description: activation process derived from the original                 //
//               G4RadioactiveDecay of F. Lei and P.R. Truscott in which      //
//               biasing and activation calculations are separated from the   //
//               unbiased decay chain calculation performed in the base       //
//               class.                                                       //  
//                                                                            //                         
////////////////////////////////////////////////////////////////////////////////

#include "G4Radioactivation.hh"
#include "G4RadioactivationMessenger.hh"

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

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace CLHEP;

G4Radioactivation::G4Radioactivation(const G4String& processName)
 : G4RadioactiveDecay(processName)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4Radioactivation constructor: processName = " << processName
           << G4endl;
  }
#endif

// DHW  SetProcessSubType(fRadioactiveDecay);
  theRadioactivationMessenger = new G4RadioactivationMessenger(this);

  // Apply default values.
  NSourceBin  = 1;
  SBin[0]     = 0.* s;
  SBin[1]     = 1.* s;    // Convert to ns
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
  AnalogueMC = true;
  BRBias = true;
  halflifethreshold = 1000.*nanosecond;
}


void G4Radioactivation::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The G4Radioactivation process performs radioactive decay of\n"
          << "nuclides (G4GenericIon) in biased mode which includes nucleus\n"
          << "duplication, branching ratio biasing, source time convolution\n"
          << "and detector time convolution.  It is designed for use in\n"
          << "activation physics.\n"
          << "The required half-lives and decay schemes are retrieved from\n"
          << "the RadioactiveDecay database which was derived from ENSDF.\n";
}


G4Radioactivation::~G4Radioactivation()
{
  delete theRadioactivationMessenger;
}

G4DecayTable* G4Radioactivation::GetDecayTable1(const G4ParticleDefinition* aNucleus)
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

G4bool
G4Radioactivation::IsRateTableReady(const G4ParticleDefinition& aParticle)
{
  // Check whether the radioactive decay rates table for the ion has already
  // been calculated.
  G4String aParticleName = aParticle.GetParticleName();
  for (std::size_t i = 0; i < theParentChainTable.size(); ++i) {
    if (theParentChainTable[i].GetIonName() == aParticleName) return true;
  }
  return false;
}


void
G4Radioactivation::GetChainsFromParent(const G4ParticleDefinition& aParticle)
{
  // Retrieve the decay rate table for the specified aParticle
  G4String aParticleName = aParticle.GetParticleName();

  for (std::size_t i = 0; i < theParentChainTable.size(); ++i) {
    if (theParentChainTable[i].GetIonName() == aParticleName) {
      theDecayRateVector = theParentChainTable[i].GetItsRates();
    }
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "The DecayRate Table for " << aParticleName << " is selected."
           <<  G4endl;
  }
#endif
}

// ConvolveSourceTimeProfile performs the convolution of the source time profile
// function with a single exponential characterized by a decay constant in the 
// decay chain.  The time profile is treated as a step function so that the 
// convolution integral can be done bin-by-bin.
// This implements Eq. 4.13 of DERA technical note, with SProfile[i] = F(t')

G4double
G4Radioactivation::ConvolveSourceTimeProfile(const G4double t, const G4double tau)
{
  G4double convolvedTime = 0.0;
  G4int nbin;
  if ( t > SBin[NSourceBin]) {
    nbin  = NSourceBin;
  } else {
    nbin = 0;

    G4int loop = 0;
    while (t > SBin[nbin]) {  // Loop checking, 01.09.2015, D.Wright
      loop++;
      if (loop > 1000) {
        G4Exception("G4Radioactivation::ConvolveSourceTimeProfile()",
                    "HAD_RDM_100", JustWarning, "While loop count exceeded");
        break;
      }
      nbin++;
    }
    nbin--;
  }

  // Use expm1 wherever possible to avoid large cancellation errors in
  // 1 - exp(x) for small x
  G4double earg = 0.0;
  if (nbin > 0) {
    for (G4int i = 0; i < nbin; i++) {
      earg = (SBin[i+1] - SBin[i])/tau;
      if (earg < 100.) {
        convolvedTime += SProfile[i] * std::exp((SBin[i] - t)/tau) *
                         std::expm1(earg);
      } else {
        convolvedTime += SProfile[i] *
          (std::exp(-(t-SBin[i+1])/tau)-std::exp(-(t-SBin[i])/tau));
      }
    }
  }
  convolvedTime -= SProfile[nbin] * std::expm1((SBin[nbin] - t)/tau);
  // tau divided out of final result to provide probability of decay in window

  if (convolvedTime < 0.)  {
    G4cout << " Convolved time =: " << convolvedTime << " reset to zero! " << G4endl;
    G4cout << " t = " << t << " tau = " << tau << G4endl;
    G4cout << SBin[nbin] << " " << SBin[0] << G4endl;
    convolvedTime = 0.;
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2)
    G4cout << " Convolved time: " << convolvedTime << G4endl;
#endif
  return convolvedTime;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetDecayTime                                                              //
//    Randomly select a decay time for the decay process, following the       //
//    supplied decay time bias scheme.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double G4Radioactivation::GetDecayTime()
{
  G4double decaytime = 0.;
  G4double rand = G4UniformRand();
  G4int i = 0;

  G4int loop = 0;
  while (DProfile[i] < rand) {  /* Loop checking, 01.09.2015, D.Wright */
    // Entries in DProfile[i] are all between 0 and 1 and arranged in inreaseing order
    // Comparison with rand chooses which time bin to sample  
    i++;
    loop++;
    if (loop > 100000) {
      G4Exception("G4Radioactivation::GetDecayTime()", "HAD_RDM_100",
                  JustWarning, "While loop count exceeded");
      break;
    }
  }

  rand = G4UniformRand();
  decaytime = DBin[i] + rand*(DBin[i+1]-DBin[i]);
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2)
    G4cout <<" Decay time: " <<decaytime/s <<"[s]" <<G4endl;
#endif
  return  decaytime;	    
}


G4int G4Radioactivation::GetDecayTimeBin(const G4double aDecayTime)
{
  G4int i = 0;

  G4int loop = 0;
  while (aDecayTime > DBin[i] ) {   /* Loop checking, 01.09.2015, D.Wright */
    i++;
    loop++;
    if (loop > 100000) {
      G4Exception("G4Radioactivation::GetDecayTimeBin()", "HAD_RDM_100", 
                  JustWarning, "While loop count exceeded");
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

G4double G4Radioactivation::GetMeanLifeTime(const G4Track& theTrack,
                                            G4ForceCondition*)
{
  // For variance reduction time is set to 0 so as to force the particle
  // to decay immediately.
  // In analogue mode it returns the particle's mean-life.
  G4double meanlife = 0.;
  if (AnalogueMC) meanlife = G4RadioactiveDecay::GetMeanLifeTime(theTrack, 0); 
  return meanlife;
}


void
G4Radioactivation::SetDecayRate(G4int theZ, G4int theA, G4double theE, 
                                 G4int theG, std::vector<G4double> theCoefficients, 
                                 std::vector<G4double> theTaos)
//  Why not make this a method of G4RadioactiveDecayRate? (e.g. SetParameters)
{ 
  //fill the decay rate vector 
  ratesToDaughter.SetZ(theZ);
  ratesToDaughter.SetA(theA);
  ratesToDaughter.SetE(theE);
  ratesToDaughter.SetGeneration(theG);
  ratesToDaughter.SetDecayRateC(theCoefficients);
  ratesToDaughter.SetTaos(theTaos);
}


void G4Radioactivation::
CalculateChainsFromParent(const G4ParticleDefinition& theParentNucleus)
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
  theDecayRateVector.push_back(ratesToDaughter);
  nEntry++;

  // Now start treating the secondary generations.
  G4bool stable = false;
//  G4int i;
  G4int j;
  G4VDecayChannel* theChannel = 0;
  G4NuclearDecay* theNuclearDecayChannel = 0;

  G4ITDecay* theITChannel = 0;
  G4BetaMinusDecay* theBetaMinusChannel = 0;
  G4BetaPlusDecay* theBetaPlusChannel = 0;
  G4AlphaDecay* theAlphaChannel = 0;
  G4ProtonDecay* theProtonChannel = 0;
  G4TritonDecay* theTritonChannel = 0;
  G4NeutronDecay* theNeutronChannel = 0;
  G4SFDecay* theFissionChannel = 0;

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
  const G4int nMode = G4RadioactiveDecayModeSize;
  G4double brs[nMode];
  //
  theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4int loop = 0;
  while (!stable) {   /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 10000) {
      G4Exception("G4Radioactivation::CalculateChainsFromParent()", "HAD_RDM_100",
                  JustWarning, "While loop count exceeded");
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
      if (GetVerboseLevel() > 1) {
        G4cout << "G4RadioactiveDecay::CalculateChainsFromParent: daughters of ("
               << ZP << ", " << AP << ", " << EP
               << ") are being calculated,  generation = " << nGeneration
               << G4endl;
      }
//      G4cout << " Taus = " << G4endl;
//      for (G4int ii = 0; ii < TP.size(); ii++) G4cout << TP[ii] << ", " ;
//      G4cout << G4endl;

      aParentNucleus = theIonTable->GetIon(ZP,AP,EP);
      parentDecayTable = GetDecayTable1(aParentNucleus);

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
      for (G4int i = 0; i < parentDecayTable->entries(); i++) {
        theChannel = parentDecayTable->GetDecayChannel(i);
        theNuclearDecayChannel = static_cast<G4NuclearDecay*>(theChannel);
        theDecayMode = theNuclearDecayChannel->GetDecayMode();
        daughterExcitation = theNuclearDecayChannel->GetDaughterExcitation();
        theDaughterNucleus = theNuclearDecayChannel->GetDaughterNucleus() ;
        AD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
        ZD = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();  
        const G4LevelManager* levelManager =
          G4NuclearLevelData::GetInstance()->GetLevelManager(ZD,AD);

        // Check each nuclide to see if it is metastable (lifetime > 1 usec) 
        // If so, add it to the decay chain by inserting its decay channel in
        // summedDecayTable.  If not, just add its BR to sum for that decay mode. 
        if (levelManager->NumberOfTransitions() ) {
          nearestEnergy = levelManager->NearestLevelEnergy(daughterExcitation);
          if (std::abs(daughterExcitation - nearestEnergy) < levelTolerance) {
            // Level half-life is in ns and the threshold is set to 1 micros
            // by default, user can set it via the UI command
            nearestLevelIndex = (G4int)levelManager->NearestLevelIndex(daughterExcitation);
            if (levelManager->LifeTime(nearestLevelIndex)*ns >= halflifethreshold){
              // save the metastable decay channel 
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

      brs[BetaPlus] = brs[BetaPlus]+brs[KshellEC]+brs[LshellEC]+brs[MshellEC]+brs[NshellEC];  // Combine beta+ and EC 
      brs[KshellEC] = brs[LshellEC] = brs[MshellEC] = brs[NshellEC] = 0.0;
      for (G4int i = 0; i < nMode; i++) {                 // loop over decay modes
        if (brs[i] > 0.) {
          switch (i) {
          case IT:
            // Decay mode is isomeric transition
            theITChannel = new G4ITDecay(aParentNucleus, brs[IT], 0.0, 0.0,
                                         photonEvaporation);

            summedDecayTable->Insert(theITChannel);
            break;

          case BetaMinus:
            // Decay mode is beta-
            theBetaMinusChannel = new G4BetaMinusDecay(aParentNucleus, brs[BetaMinus],
                                                       0.*MeV, 0.*MeV,
                                                       noFloat, allowed);
            summedDecayTable->Insert(theBetaMinusChannel);
            break;

          case BetaPlus:
            // Decay mode is beta+ + EC.
            theBetaPlusChannel = new G4BetaPlusDecay(aParentNucleus, brs[BetaPlus],
                                                     0.*MeV, 0.*MeV,
                                                     noFloat, allowed);
            summedDecayTable->Insert(theBetaPlusChannel);
            break;

          case Alpha:
            // Decay mode is alpha.
            theAlphaChannel = new G4AlphaDecay(aParentNucleus, brs[Alpha], 0.*MeV,
                                               0.*MeV, noFloat);
            summedDecayTable->Insert(theAlphaChannel);
            break;

          case Proton:
            // Decay mode is proton.
            theProtonChannel = new G4ProtonDecay(aParentNucleus, brs[Proton], 0.*MeV,
                                                 0.*MeV, noFloat);
            summedDecayTable->Insert(theProtonChannel);
            break;

          case Neutron:
            // Decay mode is neutron.
            theNeutronChannel = new G4NeutronDecay(aParentNucleus, brs[Neutron], 0.*MeV,
                                                   0.*MeV, noFloat);
            summedDecayTable->Insert(theNeutronChannel);
            break;

          case SpFission:
            // Decay mode is spontaneous fission
            theFissionChannel = new G4SFDecay(aParentNucleus, brs[SpFission], 0.*MeV,
                                              0.*MeV, noFloat);
            summedDecayTable->Insert(theFissionChannel);
            break;
	    
          case BDProton:
            // Not yet implemented
            break;

          case BDNeutron:
            // Not yet implemented
            break;

          case Beta2Minus:
            // Not yet implemented
            break;

          case Beta2Plus:
            // Not yet implemented
            break;

          case Proton2:
            // Not yet implemented
            break;

          case Neutron2:
            // Not yet implemented
            break;

          case Triton:
            // Decay mode is Triton.
            theTritonChannel = new G4TritonDecay(aParentNucleus, brs[Triton], 0.*MeV,
                                                0.*MeV, noFloat);
            summedDecayTable->Insert(theTritonChannel);
            break;

          default:
            break;
          }
        }
      }

      // loop over all branches in summedDecayTable
      //
      for (G4int i = 0; i < summedDecayTable->entries(); i++){
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
        if (IsApplicable(*theDaughterNucleus) && theBR > 0.0 && 
            aParentNucleus != theDaughterNucleus) { 
          // need to make sure daughter has decay table
          parentDecayTable = GetDecayTable1(theDaughterNucleus);
          if (parentDecayTable->entries() ) {
            A = ((const G4Ions*)(theDaughterNucleus))->GetAtomicMass();
            Z = ((const G4Ions*)(theDaughterNucleus))->GetAtomicNumber();
            E = ((const G4Ions*)(theDaughterNucleus))->GetExcitationEnergy();

            TaoPlus = theDaughterNucleus->GetPDGLifeTime();
            if (TaoPlus <= 0.)  TaoPlus = 1e-100;

            // first set the taos, one simply need to add to the parent ones
            taos.clear();
            taos = TP;   // load lifetimes of all previous generations 
            std::size_t k;
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
            theDecayRateVector.push_back(ratesToDaughter);
            nEntry++;
          } // there are entries in the table
        } // nuclide is OK to decay
      } // end of loop (i) over decay table branches

      delete summedDecayTable;

    } // Getting contents of decay rate vector (end loop on j)
    nS = nT;
    nT = nEntry;
    if (nS == nT) stable = true;
  } // while nuclide is not stable

  // end of while loop
  // the calculation completed here


  // fill the first part of the decay rate table
  // which is the name of the original particle (isotope)
  chainsFromParent.SetIonName(theParentNucleus.GetParticleName()); 

  // now fill the decay table with the newly completed decay rate vector
  chainsFromParent.SetItsRates(theDecayRateVector);

  // finally add the decayratetable to the tablevector
  theParentChainTable.push_back(chainsFromParent);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  SetSourceTimeProfile                                                      //
//     read in the source time profile function (histogram)                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void G4Radioactivation::SetSourceTimeProfile(G4String filename)
{
  std::ifstream infile ( filename, std::ios::in );
  if (!infile) {
    G4ExceptionDescription ed;
    ed << " Could not open file " << filename << G4endl; 
    G4Exception("G4Radioactivation::SetSourceTimeProfile()", "HAD_RDM_001",
                FatalException, ed);
  }

  G4double bin, flux;
  NSourceBin = -1;

  G4int loop = 0;
  while (infile >> bin >> flux) {  /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 10000) {
      G4Exception("G4Radioactivation::SetSourceTimeProfile()", "HAD_RDM_100",
                  JustWarning, "While loop count exceeded");
      break;
    }
 
    NSourceBin++;
    if (NSourceBin > 99) {
      G4Exception("G4Radioactivation::SetSourceTimeProfile()", "HAD_RDM_002",
                  FatalException, "Input source time file too big (>100 rows)");

    } else {
      SBin[NSourceBin] = bin * s;    // Convert read-in time to ns
      SProfile[NSourceBin] = flux;   // Dimensionless
    }
  }

  AnalogueMC = false;
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2)
    G4cout <<" Source Timeprofile Nbin = " << NSourceBin <<G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  SetDecayBiasProfile                                                       //
//    read in the decay bias scheme function (histogram)                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void G4Radioactivation::SetDecayBias(G4String filename)
{
  std::ifstream infile(filename, std::ios::in);
  if (!infile) G4Exception("G4Radioactivation::SetDecayBias()", "HAD_RDM_001",
                           FatalException, "Unable to open bias data file" );

  G4double bin, flux;
  G4int dWindows = 0;
  G4int i ;

  theRadioactivityTables.clear();

  NDecayBin = -1;

  G4int loop = 0;
  while (infile >> bin >> flux ) {  /* Loop checking, 01.09.2015, D.Wright */
    NDecayBin++;
    loop++;
    if (loop > 10000) {
      G4Exception("G4Radioactivation::SetDecayBias()", "HAD_RDM_100",
                  JustWarning, "While loop count exceeded");
      break;
    }
 
    if (NDecayBin > 99) {
      G4Exception("G4Radioactivation::SetDecayBias()", "HAD_RDM_002",
                  FatalException, "Input bias file too big (>100 rows)" );
    } else {
      DBin[NDecayBin] = bin * s;     // Convert read-in time to ns
      DProfile[NDecayBin] = flux;    // Dimensionless
      if (flux > 0.) {
        decayWindows[NDecayBin] = dWindows;
        dWindows++;
        G4RadioactivityTable *rTable = new G4RadioactivityTable() ;
        theRadioactivityTables.push_back(rTable);
      }
    }
  }
  for ( i = 1; i<= NDecayBin; i++) DProfile[i] += DProfile[i-1];  // Cumulative flux vs i 
  for ( i = 0; i<= NDecayBin; i++) DProfile[i] /= DProfile[NDecayBin];
                             // Normalize so entries increase from 0 to 1
  // converted to accumulated probabilities

  AnalogueMC = false;
  infile.close();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2)
    G4cout <<" Decay Bias Profile  Nbin = " << NDecayBin <<G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DecayIt                                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4VParticleChange*
G4Radioactivation::DecayIt(const G4Track& theTrack, const G4Step&)
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
        for (std::size_t i = 0; i< ValidVolumes.size(); ++i)
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
    if (GetVerboseLevel() > 1) {
      G4cout << "G4RadioactiveDecay::DecayIt : "
             << theParticleDef->GetParticleName()
             << " is not an ion or is outside (Z,A) limits set for the decay. "
             << " Set particle change accordingly. "
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

  G4DecayTable* theDecayTable = GetDecayTable1(theParticleDef);

  if (theDecayTable == 0 || theDecayTable->entries() == 0) {
    // No data in the decay table.  Set particle change parameters
    // to indicate this.
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "G4RadioactiveDecay::DecayIt : "
             << "decay table not defined for "
             << theParticleDef->GetParticleName()
             << ". Set particle change accordingly. "
             << G4endl;
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
    if (AnalogueMC) {
      G4RadioactiveDecay::DecayAnalog(theTrack);

    } else {
      // Proceed with decay using variance reduction 
      G4double energyDeposit = 0.0;
      G4double finalGlobalTime = theTrack.GetGlobalTime();
      G4double finalLocalTime = theTrack.GetLocalTime(); 
      G4int index;
      G4ThreeVector currentPosition;
      currentPosition = theTrack.GetPosition();

      G4IonTable* theIonTable;
      G4ParticleDefinition* parentNucleus;

      // Get decay chains for the given nuclide
      if (!IsRateTableReady(*theParticleDef)) CalculateChainsFromParent(*theParticleDef);
      GetChainsFromParent(*theParticleDef);

      // Declare some of the variables required in the implementation
      G4int PZ;
      G4int PA;
      G4double PE;
      G4String keyName;
      std::vector<G4double> PT;
      std::vector<G4double> PR;
      G4double taotime;
      long double decayRate;

      std::size_t i;
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

      // Now apply the nucleus splitting
      for (G4int n = 0; n < NSplit; n++) {
        // Get the decay time following the decay probability function 
        // supplied by user  
        G4double theDecayTime = GetDecayTime();
        G4int nbin = GetDecayTimeBin(theDecayTime);

        // calculate the first part of the weight function  
        G4double weight1 = 1.; 
        if (nbin == 1) {
          weight1 = 1./DProfile[nbin-1] 
                    *(DBin[nbin]-DBin[nbin-1])/NSplit;   // width of window in ns
        } else if (nbin > 1) {
          // Go from nbin to nbin-2 because flux entries in file alternate between 0 and 1 
          weight1 = 1./(DProfile[nbin]-DProfile[nbin-2])
                    *(DBin[nbin]-DBin[nbin-1])/NSplit;
          // weight1 = (probability of choosing one of the bins)*(time width of bin)/NSplit
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

          // The array of arrays theDecayRateVector contains all possible decay
          // chains of a given parent nucleus (ZP,AP,EP) to a given descendant
          // nuclide (Z,A,E).
          //
          // theDecayRateVector[0] contains the decay parameters of the parent
          // nucleus
          //           PZ = ZP
          //           PA = AP
          //           PE = EP
          //           PT[] = {TP}
          //           PR[] = {RP}
          //
          // theDecayRateVector[1] contains the decay of the parent to the first
          // generation daughter (Z1,A1,E1).
          //           PZ = Z1
          //           PA = A1
          //           PE = E1
          //           PT[] = {TP, T1}
          //           PR[] = {RP, R1}
          //
          // theDecayRateVector[2] contains the decay of the parent to the first
          // generation daughter (Z1,A1,E1) and the decay of the first
          // generation daughter to the second generation daughter (Z2,A2,E2).
          //           PZ = Z2
          //           PA = A2
          //           PE = E2
          //           PT[] = {TP, T1, T2}
          //           PR[] = {RP, R1, R2}
          //
          // theDecayRateVector[3] may contain a branch chain
          //           PZ = Z2a
          //           PA = A2a
          //           PE = E2a
          //           PT[] = {TP, T1, T2a}
          //           PR[] = {RP, R1, R2a}
          //
          // and so on.

          // Calculate the decay rate of the isotope.  decayRate is the
          // radioactivity of isotope (PZ,PA,PE) at 'theDecayTime'
          // it will be used to calculate the statistical weight of the 
          // decay products of this isotope

          // For each nuclide, calculate all the decay chains which can reach
          // the parent nuclide
          decayRate = 0.L;
          for (G4int j = 0; j < G4int(PT.size() ); j++) {
            taotime = ConvolveSourceTimeProfile(theDecayTime,PT[j]);
            decayRate -= PR[j] * (long double)taotime;  // PRs are Acoeffs, taotime is inverse time
            // Eq.4.23 of of the TN
            // note the negative here is required as the rate in the
            // equation is defined to be negative, 
            // i.e. decay away, but we need positive value here.

            // G4cout << j << "\t"<< PT[j]/s << "\t" << PR[j] << "\t" << decayRate << G4endl;		
          }

          // At this point any negative decay rates are probably small enough
          // (order 10**-30) that negative values are likely due to cancellation
          // errors.  Set them to zero.
          if (decayRate < 0.0) decayRate = 0.0;

          //  G4cout <<theDecayTime/s <<"\t"<<nbin<<G4endl;
          //  G4cout << theTrack.GetWeight() <<"\t"<<weight1<<"\t"<<decayRate<< G4endl;

          // Add isotope to the radioactivity tables
          // One table for each observation time window specified in
          // SetDecayBias(G4String filename)

          theRadioactivityTables[decayWindows[nbin-1]]
                ->AddIsotope(PZ,PA,PE,weight1*decayRate,theTrack.GetWeight());

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
            G4DecayTable* decayTable = GetDecayTable1(parentNucleus);
            ndecaych = G4int(decayTable->entries()*G4UniformRand());
            G4VDecayChannel* theDecayChannel = decayTable->GetDecayChannel(ndecaych);

            if (theDecayChannel == 0) {
              // Decay channel not found.

              if (GetVerboseLevel() > 0) {
                G4cout << " G4RadioactiveDecay::DoIt : cannot determine decay channel ";
                G4cout << " for this nucleus; decay as if no biasing active. ";
                G4cout << G4endl;
                decayTable ->DumpInfo();
              }

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
          for (index = 0; index < numberOfSecondaries; ++index) {
            asecondaryparticle = tempprods->PopProducts();
            if (asecondaryparticle->GetDefinition()->GetPDGStable() ) {
              pw.push_back(weight);
              ptime.push_back(currentTime);
              secondaryparticles.push_back(asecondaryparticle);
            }
            // Generate gammas and Xrays from excited nucleus, added by L.Desorgher
            else if (((const G4Ions*)(asecondaryparticle->GetDefinition()))
                      ->GetExcitationEnergy() > 0. && weight > 0.) {  //Compute the gamma
              G4ParticleDefinition* apartDef = asecondaryparticle->GetDefinition();
              AddDeexcitationSpectrumForBiasMode(apartDef,weight,currentTime,pw,
                                                 ptime,secondaryparticles);
            }
          }

          delete tempprods;
        } // end of i loop
      } // end of n loop 

      // now deal with the secondaries in the two stl containers
      // and submmit them back to the tracking manager
      totalNumberOfSecondaries = (G4int)pw.size();
      fParticleChangeForRadDecay.SetNumberOfSecondaries(totalNumberOfSecondaries);
      for (index=0; index < totalNumberOfSecondaries; ++index) { 
        G4Track* secondary = new G4Track(secondaryparticles[index],
                                         ptime[index], currentPosition);
        secondary->SetGoodForTrackingFlag(); 	   
        secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
        secondary->SetWeight(pw[index]); 	   
        fParticleChangeForRadDecay.AddSecondary(secondary); 
      }

      // Kill the parent particle
      fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(energyDeposit); 
      fParticleChangeForRadDecay.ProposeLocalTime(finalLocalTime);
      // Reset NumberOfInteractionLengthLeft.
      ClearNumberOfInteractionLengthLeft();
    } // end VR decay

    return &fParticleChangeForRadDecay;
  }  // end of data found branch 
} 


// Add gamma, X-ray, conversion and auger electrons for bias mode
void 
G4Radioactivation::AddDeexcitationSpectrumForBiasMode(G4ParticleDefinition* apartDef,
                                        G4double weight,G4double currentTime,
                                        std::vector<double>& weights_v,
                                        std::vector<double>& times_v,
                                        std::vector<G4DynamicParticle*>& secondaries_v)
{
  G4double elevel=((const G4Ions*)(apartDef))->GetExcitationEnergy();
  G4double life_time=apartDef->GetPDGLifeTime();
  G4ITDecay* anITChannel = 0;

  while (life_time < halflifethreshold && elevel > 0.) {
    anITChannel = new G4ITDecay(apartDef, 100., elevel, elevel, photonEvaporation);
    G4DecayProducts* pevap_products = anITChannel->DecayIt(0.);
    G4int nb_pevapSecondaries = pevap_products->entries();

    G4DynamicParticle* a_pevap_secondary = 0;
    G4ParticleDefinition* secDef = 0;
    for (G4int ind = 0; ind < nb_pevapSecondaries; ind++) {
      a_pevap_secondary= pevap_products->PopProducts();
      secDef = a_pevap_secondary->GetDefinition();

      if (secDef->GetBaryonNumber() > 4) {
        elevel = ((const G4Ions*)(secDef))->GetExcitationEnergy();
        life_time = secDef->GetPDGLifeTime();
        apartDef = secDef;
        if (secDef->GetPDGStable() ) {
          weights_v.push_back(weight);
          times_v.push_back(currentTime);
          secondaries_v.push_back(a_pevap_secondary);
        }
      } else {
        weights_v.push_back(weight);
        times_v.push_back(currentTime);
        secondaries_v.push_back(a_pevap_secondary);
      }
    }

    delete anITChannel;
    delete pevap_products;
  }
}

