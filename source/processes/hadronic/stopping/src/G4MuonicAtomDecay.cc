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
//---------------------------------------------------------------------
//
// GEANT4 Class
//
// File name: G4MuonicAtomDecay
//
// 20170522 K L Genser first implementation based on code by
// V.Ivantchenko & G4HadronicProcess & G4Decay
//
// Class Description:
//
// MuonicAtom Process where Muon either decays in orbit or is captured by the nucleus
//
// Modifications:
//
//
//------------------------------------------------------------------------

#include "G4MuonicAtomDecay.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicProcessType.hh"
#include "G4Nucleus.hh"
#include "G4ProcessManager.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4ForceCondition.hh"
#include "G4MuonicAtom.hh"
#include "G4MuonicAtomHelper.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4CascadeInterface.hh"
#include "G4MuMinusCapturePrecompound.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonicAtomDecay::G4MuonicAtomDecay(G4HadronicInteraction* hiptr,
                                     const G4String& name)
  : G4VRestDiscreteProcess(name, fDecay),
    fMuMass(G4MuonMinus::MuonMinus()->GetPDGMass()),
    cmptr(hiptr),
    verboseLevel(0)
{
  // This is not a hadronic process; assume it is a kind of decay
  enableAtRestDoIt = true;
  enablePostStepDoIt = true; // it is a streach; fixme
  theProcessSubType = 221; // (see G4DecayProcessType.hh) fixme
  if (!cmptr) {
    // cmptr = new G4CascadeInterface(); // Bertini - Pointer owned by InteractionRegistry
    cmptr = new G4MuMinusCapturePrecompound(); // Precompound
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonicAtomDecay::~G4MuonicAtomDecay()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuonicAtomDecay::IsApplicable(const G4ParticleDefinition& a)
{
  return ( a.GetParticleType() == "MuonicAtom" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void
// G4MuonicAtomDecay::PreparePhysicsTable(const G4ParticleDefinition& p)
// {
//   G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,&p);
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void G4MuonicAtomDecay::BuildPhysicsTable(const G4ParticleDefinition& p)
// {
//   G4HadronicProcessStore::Instance()->PrintInfo(&p);
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonicAtomDecay::AtRestGetPhysicalInteractionLength(
    const G4Track& aTrack, G4ForceCondition* condition)
{
  *condition = NotForced;
   // check if this is the beginning of tracking
  if (theNumberOfInteractionLengthLeft < 0.) {
    ResetNumberOfInteractionLengthLeft();
  }
  return theNumberOfInteractionLengthLeft*GetMeanLifeTime(aTrack, condition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonicAtomDecay::PostStepGetPhysicalInteractionLength(
    const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return DBL_MAX; // this will need to be changed in future; fixme
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonicAtomDecay::GetMeanLifeTime(const G4Track& aTrack,
                                            G4ForceCondition*)
{
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4MuonicAtom* muatom = static_cast<G4MuonicAtom*>(aParticleDef);
   G4double meanlife = muatom->GetPDGLifeTime();
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
     G4cout << "mean life time: "<< meanlife/ns << "[ns]" << G4endl;
   }
#endif
   return  meanlife;
}


G4VParticleChange* G4MuonicAtomDecay::DecayIt(const G4Track& aTrack,
                                              const G4Step&)
{

  // mainly based on G4HadronStoppingProcess & G4Decay
  // if primary is not Alive then do nothing
  theTotalResult.Clear(); //   G4ParticleChange*
  theTotalResult.Initialize(aTrack);
  theTotalResult.ProposeWeight(aTrack.GetWeight());
  if(aTrack.GetTrackStatus() != fAlive &&
     aTrack.GetTrackStatus() != fStopButAlive) {
    return &theTotalResult;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
  G4MuonicAtom const* muatom = static_cast<G4MuonicAtom const*>(aParticleDef);
  G4Ions const* baseion = muatom->GetBaseIon();
  G4int Z = baseion->GetAtomicNumber();
  G4double Zd = Z;
  G4double KEnergy = G4MuonicAtomHelper::GetKShellEnergy(Zd); // fixme check
  G4HadProjectile thePro(aTrack); // G4HadProjectile, here the muonic atom
  thePro.SetBoundEnergy(KEnergy);

  G4ForceCondition* condition = nullptr; // it is unused in the following call anyway
  G4double meanlife = GetMeanLifeTime(aTrack, condition);

  G4HadFinalState* result = nullptr;  // result before converting to G4VParticleChange*
  // G4int nSecondaries = 0;
  // save track time and start from zero time
  //  G4double time0 = aTrack.GetGlobalTime(); FillResult does it
  // see G4Decay DecayIt
  // see time0 down below
  thePro.SetGlobalTime(0.0);

  // do we need G4double fRemainderLifeTime; ???

  G4double maDTime =  theNumberOfInteractionLengthLeft*meanlife; //fixme check
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4MuonicAtomDecay::DecayIt time set to: "<< maDTime/ns << "[ns]" << G4endl;
  }
#endif

  // decide on DIO or Capture

  G4double lambdac = 1./muatom->GetDIOLifeTime();
  G4double lambdad = 1./muatom->GetNCLifeTime();
  G4double lambda  = lambdac + lambdad;

  if ( G4UniformRand()*lambda < lambdac) {
    //  if ( false ) { // force NC for testing

    // DIO
    // result = dmptr->ApplyYourself(thePro, *nucleus); // not quite the reaction;
    // using G4PhaseSpaceDecayChannel

#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
	G4cout << "G4MuonicAtomDecay::DecayIt: selected DIO mode" << G4endl;
      }
#endif

    // decay table; we use it only for the DIO which is more of a decay
    // code mostly copied from G4Decay

    G4DecayProducts* products = nullptr;
    G4DecayTable   *decaytable = aParticleDef->GetDecayTable();
    G4VDecayChannel* decaychannel = nullptr;
    G4double massParent = aParticle->GetMass();
    decaychannel = decaytable->SelectADecayChannel(massParent);
    if (decaychannel == nullptr) {
      // decay channel not found
      G4ExceptionDescription ed;
      ed << "Can not determine decay channel for "
         << aParticleDef->GetParticleName() << G4endl
         << "  mass of dynamic particle: " << massParent/GeV << " (GEV)" << G4endl
         << "  dacay table has " << decaytable->entries() << " entries" << G4endl;
      G4double checkedmass=massParent;
      if (massParent < 0.) {
        checkedmass=aParticleDef->GetPDGMass();
        ed << "Using PDG mass ("<<checkedmass/GeV << "(GeV)) in IsOKWithParentMass" << G4endl;	
      }
      for (G4int ic =0;ic <decaytable->entries();++ic) {
        G4VDecayChannel * dc= decaytable->GetDecayChannel(ic);
        ed << ic << ": BR " << dc->GetBR() << ", IsOK? "
           << dc->IsOKWithParentMass(checkedmass)
           << ", --> ";
        G4int ndaughters=dc->GetNumberOfDaughters();
        for (G4int id=0;id<ndaughters;++id) {
          if (id>0) ed << " + ";   // seperator, except for first
          ed << dc->GetDaughterName(id);
        }
        ed << G4endl;
      }
      G4Exception("G4MuonicAtomDecay::DecayIt", "DECAY003", FatalException,ed);
      return &theTotalResult;
    } else {
      // execute DecayIt()
#ifdef G4VERBOSE
      G4int temp = decaychannel->GetVerboseLevel();
      if (GetVerboseLevel()>1) {
	G4cout << "G4MuonicAtomDecay::DecayIt  : selected decay channel  addr:"
	       << decaychannel <<G4endl;
	decaychannel->SetVerboseLevel(GetVerboseLevel());
      }
#endif
      products = decaychannel->DecayIt(aParticle->GetMass());
      if(products == nullptr) {
	G4ExceptionDescription ed;
	ed << "No products are generated for "
           << aParticleDef->GetParticleName();
	G4Exception("G4MuonicAtomDecay::DecayIt","DECAY003",FatalException,ed);
        return &theTotalResult;
      }
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
	decaychannel->SetVerboseLevel(temp);
      }
      if (GetVerboseLevel()>2) {
	if (! products->IsChecked() ) products->DumpInfo();
      }
#endif
    }

    // get parent particle information ...................................
    G4double   ParentEnergy  = aParticle->GetTotalEnergy();
    G4double   ParentMass    = aParticle->GetMass();
    if (ParentEnergy < ParentMass) {
      if (GetVerboseLevel()>0) {
        G4cout << "G4MuonicAtomDecay::DecayIt  : Total Energy is less than its mass" << G4endl;
        G4cout << " Particle: " << aParticle->GetDefinition()->GetParticleName();
        G4cout << " Energy:"    << ParentEnergy/MeV << "[MeV]";
        G4cout << " Mass:"    << ParentMass/MeV << "[MeV]";
        G4cout << G4endl;
      }
      G4Exception( "G4MuonicAtomDecay::DecayIt ",
                   "DECAY102",JustWarning,
                   "Total Energy is less than its mass");
      ParentEnergy = ParentMass;
    }

    G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

    //boost all decay products to laboratory frame
    G4double energyDeposit = 0.0;
    G4double finalGlobalTime = aTrack.GetGlobalTime();
    G4double finalLocalTime = aTrack.GetLocalTime();
    if (aTrack.GetTrackStatus() == fStopButAlive ){
      // AtRest case
      finalGlobalTime += maDTime;
      finalLocalTime += maDTime;
      energyDeposit += aParticle->GetKineticEnergy();
    } else {
      // PostStep case
      products->Boost( ParentEnergy, ParentDirection);
    }

    //add products in theTotalResult
    G4int numberOfSecondaries = products->entries();
    theTotalResult.SetNumberOfSecondaries(numberOfSecondaries);
   
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4MuonicAtomDecay::DecayIt  : Decay vertex :";
      G4cout << " Time: " << finalGlobalTime/ns << "[ns]";
      G4cout << " X:" << (aTrack.GetPosition()).x() /cm << "[cm]";
      G4cout << " Y:" << (aTrack.GetPosition()).y() /cm << "[cm]";
      G4cout << " Z:" << (aTrack.GetPosition()).z() /cm << "[cm]";
      G4cout << G4endl;
      G4cout << "G4MuonicAtomDecay::DecayIt  : decay products in Lab. Frame" << G4endl;
      products->DumpInfo();
    }
#endif
    G4int index;
    G4ThreeVector currentPosition;
    const G4TouchableHandle thand = aTrack.GetTouchableHandle();
    for (index=0; index < numberOfSecondaries; index++)
      {
        // get current position of the track
        currentPosition = aTrack.GetPosition();
        // create a new track object
        G4Track* secondary = new G4Track( products->PopProducts(),
                                          finalGlobalTime ,
                                          currentPosition );
        // switch on good for tracking flag
        secondary->SetGoodForTrackingFlag();
        secondary->SetTouchableHandle(thand);
        // add the secondary track in the List
        theTotalResult.AddSecondary(secondary);
      }
    delete products;

    // Kill the parent particle
    theTotalResult.ProposeTrackStatus( fStopAndKill ) ;
    theTotalResult.ProposeLocalEnergyDeposit(energyDeposit);
    theTotalResult.ProposeLocalTime( finalLocalTime );

    // Clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return &theTotalResult ;

  } else { //either or

    // nuclearCapture

    // model
    // need to be able to choose between preco or bertini; no good way to do it?
    // hardcoded in the constructor for now

#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
	G4cout << "G4MuonicAtomDecay::DecayIt: selected NC  mode" << G4endl;
      }
#endif

    G4int A = baseion->GetAtomicMass();
    //    G4Nucleus* nucleus = GetTargetNucleusPointer(); // from G4HadronicProcess
    G4Nucleus nucleus;
    nucleus.SetParameters(A, Z);

    // we define a local projectile here which will be the orbiting muon
    // we shall assume it is at rest; fixme

    // G4HadProjectile, here the muon
    G4HadProjectile theMuPro(G4DynamicParticle(G4MuonMinus::MuonMinus(),
                                               G4ThreeVector(0.,0.,0.)));
    theMuPro.SetBoundEnergy(KEnergy);
    theMuPro.SetGlobalTime(0.0);

    G4int reentryCount = 0; // this may be in the model already; check fixme <---
    do {
      // sample final state
      // nuclear interaction should keep G4HadFinalState object
      // model should define time of each secondary particle
      try {
        result = cmptr->ApplyYourself(theMuPro, nucleus); // muon and muonic atom nucleus
        ++reentryCount;
      }
      catch(G4HadronicException & aR) {
        G4ExceptionDescription ed;
        ed << "Call for " << cmptr->GetModelName() << G4endl;
        ed << "  Z= "
           << nucleus.GetZ_asInt()
           << "  A= " << nucleus.GetA_asInt() << G4endl;
        DumpState(aTrack,"ApplyYourself",ed);
        ed << " ApplyYourself failed" << G4endl;
        G4Exception("G4MuonicAtomDecay::DecayIt", "HAD_MAD_101",
                    FatalException, ed);
      }

      // Check the result for catastrophic energy non-conservation
      // result = CheckResult(theMuPro, nucleus, result);

      if(reentryCount>100) {
        G4ExceptionDescription ed;
        ed << "Call for " << cmptr->GetModelName() << G4endl;
        ed << "  Z= "
           << nucleus.GetZ_asInt()
           << "  A= " << nucleus.GetA_asInt() << G4endl;
        DumpState(aTrack,"ApplyYourself",ed);
        ed << " ApplyYourself does not completed after 100 attempts" << G4endl;
        G4Exception("G4MuonicAtomDecay::DecayIt", "HAD_MAD_102",
                    FatalException, ed);
      }
      // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
    } while(result == nullptr);

    // add delay time of capture (inter + intra)
    G4int nsec = (G4int)result->GetNumberOfSecondaries();
    for(G4int i=0; i<nsec; ++i) {
      G4HadSecondary* sec = result->GetSecondary(i);
      G4double ctime = sec->GetTime();
      sec->SetTime(maDTime + ctime); // we add time0 in the next stage
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "G4MuonicAtomDecay::DecayIt time set to: "
               << (maDTime + ctime)/ns << "[ns]" << G4endl;
      }
#endif
    }

    FillResult(result,aTrack);

    ClearNumberOfInteractionLengthLeft();
    return &theTotalResult;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonicAtomDecay::ProcessDescription(std::ostream& outFile) const
{
  outFile << "MuonicAtom process where Muon decays in orbit or is captured by the nucleus." <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonicAtomDecay::FillResult(G4HadFinalState * aR, const G4Track & aT)
{
  // based on G4HadronicProcess::FillResult

  theTotalResult.ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());

  G4double rotation = CLHEP::twopi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);

  G4double efinal = aR->GetEnergyChange();
  if(efinal < 0.0) { efinal = 0.0; }

  // check status of primary
  if(aR->GetStatusChange() == stopAndKill) {
    theTotalResult.ProposeTrackStatus(fStopAndKill);
    theTotalResult.ProposeEnergy( 0.0 );

    // check its final energy
  } else if(0.0 == efinal) {
    theTotalResult.ProposeEnergy( 0.0 );
    if(aT.GetParticleDefinition()->GetProcessManager()
       ->GetAtRestProcessVector()->size() > 0)
         { theTotalResult.ProposeTrackStatus(fStopButAlive); }
    else { theTotalResult.ProposeTrackStatus(fStopAndKill); } // check fixme

    // primary is not killed apply rotation and Lorentz transformation
  } else  {
    theTotalResult.ProposeTrackStatus(fAlive);
    G4double mass = aT.GetParticleDefinition()->GetPDGMass();
    G4double newE = efinal + mass;
    G4double newP = std::sqrt(efinal*(efinal + 2*mass));
    G4ThreeVector newPV = newP*aR->GetMomentumChange();
    G4LorentzVector newP4(newE, newPV);
    newP4.rotate(rotation, it);
    newP4 *= aR->GetTrafoToLab();
    theTotalResult.ProposeMomentumDirection(newP4.vect().unit());
    newE = newP4.e() - mass;
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1 && newE <= 0.0) {
      G4ExceptionDescription ed;
      DumpState(aT,"Primary has zero energy after interaction",ed);
      G4Exception("G4MuonicAtomDecay::FillResults", "HAD_MAD_103", JustWarning, ed);
    }
#endif
    if(newE < 0.0) { newE = 0.0; }
    theTotalResult.ProposeEnergy( newE );
  }
  //G4cout << "FillResult: Efinal= " << efinal << " status= "
  //	 << theTotalResult.GetTrackStatus()
  //	 << "  fKill= " << fStopAndKill << G4endl;

  // check secondaries: apply rotation and Lorentz transformation
  G4int nSec = (G4int)aR->GetNumberOfSecondaries();
  theTotalResult.SetNumberOfSecondaries(nSec);
  G4double weight = aT.GetWeight();

  if (nSec > 0) {
    G4double time0 = aT.GetGlobalTime();
    for (G4int i = 0; i < nSec; ++i) {
      G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
      theM.rotate(rotation, it);
      theM *= aR->GetTrafoToLab();
      aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);

      // time of interaction starts from zero
      G4double time = aR->GetSecondary(i)->GetTime();
      if (time < 0.0) { time = 0.0; }

      // take into account global time
      time += time0;

      G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
                                   time, aT.GetPosition());
      track->SetCreatorModelID(aR->GetSecondary(i)->GetCreatorModelID());
      G4double newWeight = weight*aR->GetSecondary(i)->GetWeight();
	// G4cout << "#### ParticleDebug "
	// <<GetProcessName()<<" "
	//<<aR->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()<<" "
	// <<aScaleFactor<<" "
	// <<XBiasSurvivalProbability()<<" "
	// <<XBiasSecondaryWeight()<<" "
	// <<aT.GetWeight()<<" "
	// <<aR->GetSecondary(i)->GetWeight()<<" "
	// <<aR->GetSecondary(i)->GetParticle()->Get4Momentum()<<" "
	// <<G4endl;
      track->SetWeight(newWeight);
      track->SetTouchableHandle(aT.GetTouchableHandle());
      theTotalResult.AddSecondary(track);
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4double e = track->GetKineticEnergy();
        if (e <= 0.0) {
          G4ExceptionDescription ed;
          DumpState(aT,"Secondary has zero energy",ed);
          ed << "Secondary " << track->GetDefinition()->GetParticleName()
             << G4endl;
          G4Exception("G4MuonicAtomDecay::FillResults", "HAD_MAD_103",
		      JustWarning,ed);
        }
      }
#endif
    }
  }
  aR->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonicAtomDecay::DumpState(const G4Track& aTrack,
				  const G4String& method,
				  G4ExceptionDescription& ed)
{
  ed << "Unrecoverable error in the method " << method << " of "
     << GetProcessName() << G4endl;
  ed << "TrackID= "<< aTrack.GetTrackID() << "  ParentID= "
     << aTrack.GetParentID()
     << "  " << aTrack.GetParticleDefinition()->GetParticleName()
     << G4endl;
  ed << "Ekin(GeV)= " << aTrack.GetKineticEnergy()/CLHEP::GeV
     << ";  direction= " << aTrack.GetMomentumDirection() << G4endl;
  ed << "Position(mm)= " << aTrack.GetPosition()/CLHEP::mm << ";";

  if (aTrack.GetMaterial()) {
    ed << "  material " << aTrack.GetMaterial()->GetName();
  }
  ed << G4endl;

  if (aTrack.GetVolume()) {
    ed << "PhysicalVolume  <" << aTrack.GetVolume()->GetName()
       << ">" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonicAtomDecay::GetMeanFreePath(const G4Track& aTrack,G4double, G4ForceCondition*)
{
  // based on G4Decay::GetMeanFreePath check; fixme

  // get particle
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
  G4double aMass = aParticle->GetMass();
  G4double aLife = aParticleDef->GetPDGLifeTime();

  // returns the mean free path in GEANT4 internal units
  G4double pathlength;
  G4double aCtau = c_light * aLife;

  // check if the particle is stable?
  if (aParticleDef->GetPDGStable()) {
    pathlength = DBL_MAX;

    //check if the particle has very short life time ?
  } else if (aCtau < DBL_MIN) {
    pathlength =  DBL_MIN;

  } else {
    //calculate the mean free path
    // by using normalized kinetic energy (= Ekin/mass)
    G4double   rKineticEnergy = aParticle->GetKineticEnergy()/aMass;
    const G4double HighestValue = 20.0; //
    if ( rKineticEnergy > HighestValue) {
      // gamma >>  1
      pathlength = ( rKineticEnergy + 1.0)* aCtau;
    } else if ( rKineticEnergy < DBL_MIN ) {
      // too slow particle
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "G4MuonicAtomDecay::GetMeanFreePath()   !!particle stops!!";
        G4cout << aParticleDef->GetParticleName() << G4endl;
        G4cout << "KineticEnergy:" << aParticle->GetKineticEnergy()/GeV <<"[GeV]";
      }
#endif
      pathlength = DBL_MIN;
    } else {
      // beta <1
      pathlength = (aParticle->GetTotalMomentum())/aMass*aCtau ;
    }
  }
  return  pathlength;
}
