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

#include "G4Types.hh"

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "G4HadronicProcess.hh"
// #include "G4EffectiveCharge.hh"
#include "G4HadProjectile.hh"
#include "G4ElementVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Element.hh"
#include "G4ParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4StableIsotopes.hh"
#include "G4HadTmpUtil.hh"

#include "G4HadLeadBias.hh"
#include "G4HadronicException.hh"
#include "G4HadReentrentException.hh"
#include "G4HadronicInteractionWrapper.hh"

#include "G4HadSignalHandler.hh"

#include <typeinfo>

namespace G4HadronicProcess_local
{
  extern "C" void G4HadronicProcessHandler_1(int)
  {
    G4HadronicWhiteBoard::Instance().Dump();
  }
} 

G4IsoParticleChange * G4HadronicProcess::theIsoResult = 0;
G4IsoParticleChange * G4HadronicProcess::theOldIsoResult = 0;
G4bool G4HadronicProcess::isoIsEnabled = true;

void G4HadronicProcess::
EnableIsotopeProductionGlobally()  {isoIsEnabled = true;}

void G4HadronicProcess::
DisableIsotopeProductionGlobally() {isoIsEnabled = false;}

G4HadronicProcess::G4HadronicProcess( const G4String &processName,
                                      G4ProcessType   aType ) :
G4VDiscreteProcess( processName, aType)
{ 
  ModelingState = 0;
  isoIsOnAnyway = -1;
  theTotalResult = new G4ParticleChange();
  theCrossSectionDataStore = new G4CrossSectionDataStore();
  aScaleFactor = 1;
  xBiasOn = false;
  if(getenv("SwitchLeadBiasOn")) theBias.push_back(new G4HadLeadBias());
}

G4HadronicProcess::~G4HadronicProcess()
{ 
  delete theTotalResult;

  std::for_each(theProductionModels.begin(),
                theProductionModels.end(), G4Delete());
  std::for_each(theBias.begin(), theBias.end(), G4Delete());
 
  delete theOldIsoResult; delete theIsoResult;
  delete theCrossSectionDataStore;
}

void G4HadronicProcess::RegisterMe( G4HadronicInteraction *a )
{ 
  try{GetManagerPointer()->RegisterMe( a );}   
  catch(G4HadronicException & aE)
  {
    aE.Report(std::cout);
    G4Exception("G4HadronicProcess", "007", FatalException,
    "Could not register G4HadronicInteraction");
  }
}

G4double G4HadronicProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{ 
  G4double sigma = 0.0;
  try
  {
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    if( !IsApplicable(*aParticle->GetDefinition()))
    {
      G4cout << "Unrecoverable error: "<<G4endl;
      G4ProcessManager * it = aParticle->GetDefinition()->GetProcessManager();
      G4ProcessVector * itv = it->GetProcessList();
      G4cout <<aParticle->GetDefinition()->GetParticleName()<< 
      " has the following processes:"<<G4endl;
      for(G4int i=0; i<itv->size(); i++)
      {
        G4cout <<"  "<<(*itv)[i]->GetProcessName()<<G4endl;		 
      }
      G4cout << "for kinetic energy "<<aParticle->GetKineticEnergy()<<G4endl;
      G4cout << "and material "<<aTrack.GetMaterial()->GetName()<<G4endl;
      G4Exception("G4HadronicProcess", "007", FatalException,
      std::string(this->GetProcessName()+
      " was called for "+
      aParticle->GetDefinition()->GetParticleName()).c_str() );
    }
    G4Material *aMaterial = aTrack.GetMaterial();
    ModelingState = 1;
    
    sigma = theCrossSectionDataStore->GetCrossSection(aParticle, aMaterial);

    sigma *= aScaleFactor;
    theLastCrossSection = sigma;
  }
  catch(G4HadronicException aR)
  { 
    aR.Report(G4cout); 
    G4Exception("G4HadronicProcess", "007", FatalException,
    "G4HadronicProcess::GetMeanFreePath failed");
  } 
  if( sigma > 0.0 )
    return 1.0/sigma;
  else
    return DBL_MAX;
}


G4Element* G4HadronicProcess::ChooseAandZ(
const G4DynamicParticle *aParticle, const G4Material *aMaterial )
{
  std::pair<G4double, G4double> ZA = 
      theCrossSectionDataStore->SelectRandomIsotope(aParticle, aMaterial);
  G4double ZZ = ZA.first;
  G4double AA = ZA.second;

  targetNucleus.SetParameters(AA, ZZ);

  const G4int numberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  G4Element* chosen = 0;
  for (G4int i = 0; i < numberOfElements; i++) {
    chosen = (*theElementVector)[i];
    if (chosen->GetZ() == ZZ) break;
  }
  return chosen;
}


struct G4Nancheck{ bool operator()(G4double aV){return (!(aV<1))&&(!(aV>-1));}};

G4VParticleChange *G4HadronicProcess::GeneralPostStepDoIt(
const G4Track &aTrack, const G4Step &)
{
  // Debugging stuff

  bool G4HadronicProcess_debug_flag = false;
  if(getenv("G4HadronicProcess_debug")) G4HadronicProcess_debug_flag = true;
  if(G4HadronicProcess_debug_flag) 
         std::cout << "@@@@ hadronic process start "<< std::endl;
  // G4cout << theNumberOfInteractionLengthLeft<<G4endl;
  #ifndef G4HadSignalHandler_off
  G4HadSignalHandler aHandler(G4HadronicProcess_local::G4HadronicProcessHandler_1);
  #endif

  if(aTrack.GetTrackStatus() != fAlive && aTrack.GetTrackStatus() != fSuspend) 
  {
    G4cerr << "G4HadronicProcess: track in unusable state - "
    <<aTrack.GetTrackStatus()<<G4endl;
    G4cerr << "G4HadronicProcess: returning unchanged track "<<G4endl;
    G4Exception("G4HadronicProcess", "001", JustWarning, "bailing out");
    theTotalResult->Clear();
    theTotalResult->Initialize(aTrack);
    return theTotalResult;
  }

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  G4Material *aMaterial = aTrack.GetMaterial();
  G4double originalEnergy = aParticle->GetKineticEnergy();
  G4double kineticEnergy = originalEnergy;

  // More debugging

  G4Nancheck go_wild;
  if(go_wild(originalEnergy) ||
    go_wild(aParticle->Get4Momentum().x()) ||
  go_wild(aParticle->Get4Momentum().y()) ||
  go_wild(aParticle->Get4Momentum().z()) ||
  go_wild(aParticle->Get4Momentum().t())
  )
  {
    G4Exception("G4HadronicProcess", "001", JustWarning, "NaN in input energy or momentum - bailing out.");
    theTotalResult->Clear();
    theTotalResult->Initialize(aTrack);
    return theTotalResult;
  }

  // Get kinetic energy per nucleon for ions

  if(aParticle->GetDefinition()->GetBaryonNumber() > 1.5)
          kineticEnergy/=aParticle->GetDefinition()->GetBaryonNumber();

  G4Element* anElement = 0;
  try
  {
    anElement = ChooseAandZ( aParticle, aMaterial );
  }
  catch(G4HadronicException & aR)
  {
    aR.Report(G4cout);
    G4cout << "Unrecoverable error for:"<<G4endl;
    G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
    G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
    G4cout << " - Particle type = "
    <<aParticle->GetDefinition()->GetParticleName()<<G4endl;
    G4Exception("G4HadronicProcess", "007", FatalException,
    "GeneralPostStepDoIt failed on element selection.");
  }

  try
  {
    theInteraction = ChooseHadronicInteraction( kineticEnergy,
    aMaterial, anElement );
  }
  catch(G4HadronicException & aE)
  {
    aE.Report(std::cout);
    G4cout << "Unrecoverable error for:"<<G4endl;
    G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
    G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
    G4cout << " - Particle type = " 
           << aParticle->GetDefinition()->GetParticleName()<<G4endl;
    G4Exception("G4HadronicProcess", "007", FatalException,
    "ChooseHadronicInteraction failed.");
  }

  // Initialize the hadronic projectile from the track

  G4HadProjectile thePro(aTrack);
  
  G4HadFinalState* result = 0;
  G4int reentryCount = 0;

  do
  {
    try
    {
      // Call the interaction

      G4HadronicInteractionWrapper aW;
      result = aW.ApplyInteraction(thePro, targetNucleus, theInteraction,
                                   GetProcessName(),
                                   theInteraction->GetModelName());
    }
    catch(G4HadReentrentException aR)
    {
      aR.Report(G4cout);
      G4cout << " G4HadronicProcess re-entering the ApplyYourself call for "
             <<G4endl;
      G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
      G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
      G4cout << " - Particle type = " 
             << aParticle->GetDefinition()->GetParticleName() << G4endl;
      result = 0; // here would still be leaking...
      if(reentryCount>100)
      {
        G4Exception("G4HadronicProcess", "007", FatalException,
        "GetHadronicProcess: Reentering ApplyYourself too often - GeneralPostStepDoIt failed.");  
      }
      G4Exception("G4HadronicProcess", "007", FatalException,
      "GetHadronicProcess: GeneralPostStepDoIt failed (Reentering ApplyYourself not yet supported.)");  
    }
    catch(G4HadronicException aR)
    {
      aR.Report(G4cout);
      G4cout << " G4HadronicProcess failed in ApplyYourself call for" 
             << G4endl;
      G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
      G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
      G4cout << " - Particle type = " 
             << aParticle->GetDefinition()->GetParticleName() << G4endl;
      G4Exception("G4HadronicProcess", "007", FatalException,
      "GeneralPostStepDoIt failed.");
    }
  }
  while(!result);

  if(!ModelingState && !getenv("BypassAllSafetyChecks") ) 
  {
    G4cout << "ERROR IN EXECUTION -- HADRONIC PROCESS STATE NOT VALID"<<G4endl;
    G4cout << "Result will be of undefined quality."<<G4endl;
  }

  // NOT USED ?? Projectile particle has changed character during interaction
  if(result->GetStatusChange() == isAlive && 
     thePro.GetDefinition() != aTrack.GetDefinition())
  {
    G4DynamicParticle * aP = 
             const_cast<G4DynamicParticle *>(aTrack.GetDynamicParticle());
    aP->SetDefinition(const_cast<G4ParticleDefinition *>(thePro.GetDefinition()));
  }

  result->SetTrafoToLab(thePro.GetTrafoToLab());

  /*
  // Loop over charged ion secondaries

  for(G4int i=0; i<result->GetNumberOfSecondaries(); i++)
  {
    G4DynamicParticle* aSecTrack = result->GetSecondary(i)->GetParticle();
    if(aSecTrack->GetDefinition()->GetPDGCharge()>1.5)
    {
      G4EffectiveCharge aCalculator;
      G4double charge = 
           aCalculator.GetCharge(aMaterial, aSecTrack->GetKineticEnergy(),
                                 aSecTrack->GetDefinition()->GetPDGMass(),
                                 aSecTrack->GetDefinition()->GetPDGCharge());
      if(getenv("GHADChargeDebug")) 
      {
        std::cout << "Recoil fractional charge is "
        << charge/aSecTrack->GetDefinition()->GetPDGCharge()<<" "
        << charge <<" "<<aSecTrack->GetDefinition()->GetPDGCharge()<<std::endl;
      }
      aSecTrack->SetCharge(charge);
    }
  }
  */

  if(getenv("HadronicDoitLogging") )
  {
    G4cout << "HadronicDoitLogging "
    << GetProcessName() <<" "
    << aParticle->GetDefinition()->GetPDGEncoding()<<" "
    << originalEnergy<<" "
    << aParticle->GetMomentum()<<" "
    << targetNucleus.GetN()<<" "
    << targetNucleus.GetZ()<<" "
    << G4endl;
  }

  ClearNumberOfInteractionLengthLeft();
  if(isoIsOnAnyway!=-1)
  {
    if(isoIsEnabled||isoIsOnAnyway)
    {
      result = DoIsotopeCounting(result, aTrack, targetNucleus);
    }
  }
  
  G4double e=aTrack.GetKineticEnergy();
  ModelingState = 0;
  if(e<5*GeV)
  {
    for(size_t i=0; i<theBias.size(); i++)
    {
      result = theBias[i]->Bias(result);
    }
  }

  // Put hadronic final state particles into G4ParticleChange

  FillTotalResult(result, aTrack);
  if(G4HadronicProcess_debug_flag) 
    std::cout << "@@@@ hadronic process end "<< std::endl;
    
  return theTotalResult;
}


G4HadFinalState* 
G4HadronicProcess::DoIsotopeCounting(G4HadFinalState * aResult,
                                     const G4Track & aTrack,
                                     const G4Nucleus & aNucleus)
{
  // get the PC from iso-production
  delete theOldIsoResult;
  theOldIsoResult = 0;
  delete theIsoResult;
  theIsoResult = new G4IsoParticleChange;
  G4bool done = false;
  G4IsoResult * anIsoResult = 0;
  for(unsigned int i=0; i<theProductionModels.size(); i++)
  {
    anIsoResult = theProductionModels[i]->GetIsotope(aTrack, aNucleus);
    if(anIsoResult!=0)
    {
      done = true;
      break;
    }
  }

  // If no production models active, use default iso production
  if(!done) anIsoResult = ExtractResidualNucleus(aTrack, aNucleus, aResult); 

  // Add all info explicitely and add typename from model called.
  theIsoResult->SetIsotope(anIsoResult->GetIsotope());
  theIsoResult->SetProductionPosition(aTrack.GetPosition());
  theIsoResult->SetProductionTime(aTrack.GetGlobalTime());
  theIsoResult->SetParentParticle(*aTrack.GetDynamicParticle());
  theIsoResult->SetMotherNucleus(anIsoResult->GetMotherNucleus());
  theIsoResult->SetProducer(typeid(*theInteraction).name());
  
  delete anIsoResult;

  // If isotope production is enabled the GetIsotopeProductionInfo() 
  // method must be called or else a memory leak will result
  //
  // The following code will fix the memory leak, but remove the 
  // isotope information:
  //
  //  if(theIsoResult) {
  //    delete theIsoResult;
  //    theIsoResult = 0;
  //  }
  
  return aResult;
}

G4IsoResult* 
G4HadronicProcess::ExtractResidualNucleus(const G4Track&,
                                          const G4Nucleus& aNucleus,
                                          G4HadFinalState* aResult)
{
  G4double A = aNucleus.GetN();
  G4double Z = aNucleus.GetZ();
  G4double bufferA = 0;
  G4double bufferZ = 0;
  
  // loop over aResult, and decrement A, Z accordingly
  // cash the max
  for(G4int i=0; i<aResult->GetNumberOfSecondaries(); i++)
  {
    G4HadSecondary* aSecTrack = aResult->GetSecondary(i);
    if(bufferA<aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber())
    {
      bufferA = aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber();
      bufferZ = aSecTrack->GetParticle()->GetDefinition()->GetPDGCharge();
    }
    Z-=aSecTrack->GetParticle()->GetDefinition()->GetPDGCharge();
    A-=aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber();
  }
  
  // if the fragment was part of the final state, it is 
  // assumed to be the heaviest secondary.
  if(A<0.1)
  {
    A = bufferA;
    Z = bufferZ;
  }
  
  // prepare the IsoResult.

  std::ostringstream ost1;
  ost1 <<Z<<"_"<<A;
  G4String biff = ost1.str();
  G4IsoResult * theResult = new G4IsoResult(biff, aNucleus);

  return theResult;
}

G4double G4HadronicProcess::XBiasSurvivalProbability()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double biasedProbability = 1.-std::exp(-nLTraversed);
  G4double realProbability = 1-std::exp(-nLTraversed/aScaleFactor);
  result = (biasedProbability-realProbability)/biasedProbability;
  return result;
}

G4double G4HadronicProcess::XBiasSecondaryWeight()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  result = 
     1./aScaleFactor*std::exp(-nLTraversed/aScaleFactor*(1-1./aScaleFactor));
  return result;
}

void 
G4HadronicProcess::FillTotalResult(G4HadFinalState * aR, const G4Track & aT)
{
  G4Nancheck go_wild;
  theTotalResult->Clear();
  theTotalResult->ProposeLocalEnergyDeposit(0.);
  theTotalResult->Initialize(aT);
  theTotalResult->SetSecondaryWeightByProcess(true);
  theTotalResult->ProposeTrackStatus(fAlive);
  G4double rotation = 2.*pi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);

  /*
  if(xBiasOn)
  {
    G4cout << "BiasDebug "<<GetProcessName()<<" "
    <<aScaleFactor<<" "
    <<XBiasSurvivalProbability()<<" "
    <<XBiasSecondaryWeight()<<" "
    <<G4endl;
  }
  */
  // if(GetProcessName() != "LElastic") std::cout << "Debug -1 "<<aR->GetStatusChange()<<std::endl;
  if(aR->GetStatusChange()==stopAndKill)
  {
    if( xBiasOn && G4UniformRand()<XBiasSurvivalProbability() )
    {
      theTotalResult->ProposeParentWeight( XBiasSurvivalProbability()*aT.GetWeight() );
    }
    else
    {
      theTotalResult->ProposeTrackStatus(fStopAndKill);
      theTotalResult->ProposeEnergy( 0.0 );
    }
  }
  else if(aR->GetStatusChange()!=stopAndKill )
  {
    if(aR->GetStatusChange()==suspend)
    {
      theTotalResult->ProposeTrackStatus(fSuspend);
      if(xBiasOn)
      {
        G4Exception("G4HadronicProcess", "007", FatalException,
        "Cannot cross-section bias a process that suspends tracks.");
      }
    } else if (aT.GetKineticEnergy() == 0) {
      theTotalResult->ProposeTrackStatus(fStopButAlive);
    }

    if(xBiasOn && G4UniformRand()<XBiasSurvivalProbability())
    {
      theTotalResult->ProposeParentWeight( XBiasSurvivalProbability()*aT.GetWeight() );
      G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
      if(go_wild(aR->GetEnergyChange()))
      {
        G4Exception("G4HadronicProcess", "007", FatalException,
        "surviving track received NaN energy.");  
      }
      if(go_wild(aR->GetMomentumChange().x()) || 
        go_wild(aR->GetMomentumChange().y()) || 
      go_wild(aR->GetMomentumChange().z()))
      {
        G4Exception("G4HadronicProcess", "007", FatalException,
        "surviving track received NaN momentum.");  
      }
      G4double newM=aT.GetDefinition()->GetPDGMass();
      G4double newE=aR->GetEnergyChange() + newM;
      G4double newP=std::sqrt(newE*newE - newM*newM);
      G4DynamicParticle * aNew = 
      new G4DynamicParticle(aT.GetDefinition(), newE, newP*aR->GetMomentumChange());
      G4HadSecondary * theSec = new G4HadSecondary(aNew, newWeight);
      aR->AddSecondary(theSec);
    }
    else
    {
      G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
      theTotalResult->ProposeParentWeight(newWeight); // This is multiplicative
      if(aR->GetEnergyChange()>-.5) 
      {
        if(go_wild(aR->GetEnergyChange()))
        {
          G4Exception("G4HadronicProcess", "007", FatalException,
          "track received NaN energy.");  
        }
        theTotalResult->ProposeEnergy(aR->GetEnergyChange());
      }
      G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
      newDirection*=aR->GetTrafoToLab();
      theTotalResult->ProposeMomentumDirection(newDirection.vect());
    }
  }
  else
  {
    G4cerr << "Track status is "<< aR->GetStatusChange()<<G4endl;
    G4Exception("G4HadronicProcess", "007", FatalException,
    "use of unsupported track-status.");
  }

  if(GetProcessName() != "hElastic" && GetProcessName() != "HadronElastic"
     &&  theTotalResult->GetTrackStatus()==fAlive
     && aR->GetStatusChange()==isAlive
    )
  {
    // Use for debugging:   G4double newWeight = theTotalResult->GetParentWeight();

    G4double newKE = std::max(DBL_MIN, aR->GetEnergyChange());
    G4DynamicParticle* aNew = new G4DynamicParticle(aT.GetDefinition(), 
                                                    aR->GetMomentumChange(), 
                                                    newKE);
    G4HadSecondary* theSec = new G4HadSecondary(aNew, 1.0);
    aR->AddSecondary(theSec);
    aR->SetStatusChange(stopAndKill);

    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeEnergy( 0.0 );

  }
  theTotalResult->ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
  theTotalResult->SetNumberOfSecondaries(aR->GetNumberOfSecondaries());

  if(aR->GetStatusChange() != stopAndKill)
  {
    G4double newM=aT.GetDefinition()->GetPDGMass();
    G4double newE=aR->GetEnergyChange() + newM;
    G4double newP=std::sqrt(newE*newE - newM*newM);
    G4ThreeVector newPV = newP*aR->GetMomentumChange();
    G4LorentzVector newP4(newE, newPV);
    newP4.rotate(rotation, it);
    newP4*=aR->GetTrafoToLab();
    theTotalResult->ProposeMomentumDirection(newP4.vect().unit());
  }

  for(G4int i=0; i<aR->GetNumberOfSecondaries(); i++)
  {
    G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
    theM.rotate(rotation, it);
    theM*=aR->GetTrafoToLab();

    if(go_wild(theM.e()))
    {
      G4Exception("G4HadronicProcess", "007", FatalException,
      "secondary track received NaN energy.");  
    }
    if(go_wild(theM.x()) || 
      go_wild(theM.y()) || 
    go_wild(theM.z()))
    {
      G4Exception("G4HadronicProcess", "007", FatalException,
      "secondary track received NaN momentum.");  
    }

    aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);
    G4double time = aR->GetSecondary(i)->GetTime();
    if(time<0) time = aT.GetGlobalTime();

    G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
    aT.GetGlobalTime(),
    aT.GetPosition());

    G4double newWeight = aT.GetWeight()*aR->GetSecondary(i)->GetWeight();
    //static G4double pinelcount=0;
    if(xBiasOn) newWeight *= XBiasSecondaryWeight();
    // G4cout << "#### ParticleDebug "
    // <<GetProcessName()<<" "
    // <<aR->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()<<" "
    // <<aScaleFactor<<" "
    // <<XBiasSurvivalProbability()<<" "
    // <<XBiasSecondaryWeight()<<" "
    // <<aT.GetWeight()<<" "
    // <<aR->GetSecondary(i)->GetWeight()<<" "
    // <<aR->GetSecondary(i)->GetParticle()->Get4Momentum()<<" "
    // <<G4endl;
    track->SetWeight(newWeight);
    G4double trackDeb = track->GetKineticEnergy();
    if( (  trackDeb<0 
      || (trackDeb>aT.GetKineticEnergy()+1*GeV) ) && getenv("GHADEnergyBalanceDebug") )
      {
        G4cout << "Debugging hadronic processes: "<<track->GetKineticEnergy()
        <<" "<<aT.GetKineticEnergy()
        <<" "<<GetProcessName()
        <<" "<<aT.GetDefinition()->GetParticleName() 
        <<G4endl;
      }
      track->SetTouchableHandle(aT.GetTouchableHandle());
      theTotalResult->AddSecondary(track);
  }

  aR->Clear();
  return;
}
/* end of file */
