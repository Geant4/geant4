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
//    G4AntiNeutronAnnihilationAtRest physics process
//    Larry Felawka (TRIUMF), April 1998
//---------------------------------------------------------------------

#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicDeprecate.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
 
#define MAX_SECONDARIES 100

// constructor
 
G4AntiNeutronAnnihilationAtRest::G4AntiNeutronAnnihilationAtRest(const G4String& processName,
                                                                 G4ProcessType aType) :
  G4VRestProcess (processName, aType),       // initialization
  massPionMinus(G4PionMinus::PionMinus()->GetPDGMass()/GeV),
  massPionZero(G4PionZero::PionZero()->GetPDGMass()/GeV),
  massPionPlus(G4PionPlus::PionPlus()->GetPDGMass()/GeV),
  massGamma(G4Gamma::Gamma()->GetPDGMass()/GeV),
  massAntiNeutron(G4AntiNeutron::AntiNeutron()->GetPDGMass()/GeV),
  massNeutron(G4Neutron::Neutron()->GetPDGMass()/GeV),
  pdefGamma(G4Gamma::Gamma()),
  pdefPionPlus(G4PionPlus::PionPlus()),
  pdefPionZero(G4PionZero::PionZero()),
  pdefPionMinus(G4PionMinus::PionMinus()),
  pdefProton(G4Proton::Proton()),
  pdefNeutron(G4Neutron::Neutron()),
  pdefAntiNeutron(G4AntiNeutron::AntiNeutron()),
  pdefDeuteron(G4Deuteron::Deuteron()),
  pdefTriton(G4Triton::Triton()),
  pdefAlpha(G4Alpha::Alpha())
{
  G4HadronicDeprecate("G4AntiNeutronAnnihilationAtRest");
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessSubType(fHadronAtRest);
  pv   = new G4GHEKinematicsVector [MAX_SECONDARIES+1];
  eve  = new G4GHEKinematicsVector [MAX_SECONDARIES];
  gkin = new G4GHEKinematicsVector [MAX_SECONDARIES];

  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
  globalTime = targetAtomicMass = targetCharge = evapEnergy1 
    = evapEnergy3 = 0.0;
  ngkine = ntot = 0;
}
 
// destructor
 
G4AntiNeutronAnnihilationAtRest::~G4AntiNeutronAnnihilationAtRest()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
  delete [] pv;
  delete [] eve;
  delete [] gkin;
}
 
void G4AntiNeutronAnnihilationAtRest::PreparePhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4AntiNeutronAnnihilationAtRest::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}
 
// methods.............................................................................
 
G4bool G4AntiNeutronAnnihilationAtRest::IsApplicable(
				 const G4ParticleDefinition& particle
				 )
{
   return ( &particle == pdefAntiNeutron );

}

// Warning - this method may be optimized away if made "inline"
G4int G4AntiNeutronAnnihilationAtRest::GetNumberOfSecondaries()
{
  return ( ngkine );

}

// Warning - this method may be optimized away if made "inline"
G4GHEKinematicsVector* G4AntiNeutronAnnihilationAtRest::GetSecondaryKinematics()
{
  return ( &gkin[0] );

}

G4double G4AntiNeutronAnnihilationAtRest::AtRestGetPhysicalInteractionLength(
				   const G4Track& track,
				   G4ForceCondition* condition
				   )
{
  // beggining of tracking 
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean life time
  currentInteractionLength = GetMeanLifeTime(track, condition);

  if ((currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4AntiNeutronAnnihilationAtRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }

  return theNumberOfInteractionLengthLeft * currentInteractionLength;

}

G4VParticleChange* G4AntiNeutronAnnihilationAtRest::AtRestDoIt(
                                            const G4Track& track,
					    const G4Step& 
					    )
//
// Handles AntiNeutrons at rest; an AntiNeutron can either create secondaries
// or do nothing (in which case it should be sent back to decay-handling
// section
//
{

//   Initialize ParticleChange
//     all members of G4VParticleChange are set to equal to 
//     corresponding member in G4Track

  aParticleChange.Initialize(track);

//   Store some global quantities that depend on current material and particle

  globalTime = track.GetGlobalTime()/s;
  G4Material * aMaterial = track.GetMaterial();
  const G4int numberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  const G4double* theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
  G4double normalization = 0;
  for ( G4int i1=0; i1 < numberOfElements; i1++ )
  {
    normalization += theAtomicNumberDensity[i1] ; // change when nucleon specific
                                                  // probabilities are included.
  }
  G4double runningSum= 0.;
  G4double random = G4UniformRand()*normalization;
  for ( G4int i2=0; i2 < numberOfElements; i2++ )
  {
    runningSum += theAtomicNumberDensity[i2]; // change when nucleon specific
                                              // probabilities are included.
    if (random<=runningSum)
    {
      targetCharge = G4double( ((*theElementVector)[i2])->GetZ());
      targetAtomicMass = (*theElementVector)[i2]->GetN();
    }
  }
  if (random>runningSum)
  {
    targetCharge = G4double( ((*theElementVector)[numberOfElements-1])->GetZ());
    targetAtomicMass = (*theElementVector)[numberOfElements-1]->GetN();
  }

  if (verboseLevel>1) {
    G4cout << "G4AntiNeutronAnnihilationAtRest::AtRestDoIt is invoked " <<G4endl;
    }

  G4ParticleMomentum momentum;
  G4float localtime;

  G4ThreeVector   position = track.GetPosition();

  GenerateSecondaries(); // Generate secondaries

  aParticleChange.SetNumberOfSecondaries( ngkine ); 

  for ( G4int isec = 0; isec < ngkine; isec++ ) {
    G4DynamicParticle* aNewParticle = new G4DynamicParticle;
    aNewParticle->SetDefinition( gkin[isec].GetParticleDef() );
    aNewParticle->SetMomentum( gkin[isec].GetMomentum() * GeV );

    localtime = globalTime + gkin[isec].GetTOF();

    G4Track* aNewTrack = new G4Track( aNewParticle, localtime*s, position );
		aNewTrack->SetTouchableHandle(track.GetTouchableHandle());
    aParticleChange.AddSecondary( aNewTrack ); 

  }

  aParticleChange.ProposeLocalEnergyDeposit( 0.0*GeV );

  aParticleChange.ProposeTrackStatus(fStopAndKill); // Kill the incident AntiNeutron

//   clear InteractionLengthLeft

  ResetNumberOfInteractionLengthLeft();

  return &aParticleChange;

}


void G4AntiNeutronAnnihilationAtRest::GenerateSecondaries()
{
  G4int index;
  G4int l;
  G4int nopt;
  G4int i;
  // DHW 15 May 2011: unused: static G4ParticleDefinition* jnd;

  for (i = 1; i <= MAX_SECONDARIES; ++i) {
    pv[i].SetZero();
  }


  ngkine = 0;            // number of generated secondary particles
  ntot = 0;
  result.SetZero();
  result.SetMass( massAntiNeutron );
  result.SetKineticEnergyAndUpdate( 0. );
  result.SetTOF( 0. );
  result.SetParticleDef( pdefAntiNeutron );

  // *** SELECT PROCESS FOR CURRENT PARTICLE ***

  AntiNeutronAnnihilation(&nopt);

  // *** CHECK WHETHER THERE ARE NEW PARTICLES GENERATED ***
  if (ntot != 0 || result.GetParticleDef() != pdefAntiNeutron) {
    // *** CURRENT PARTICLE IS NOT THE SAME AS IN THE BEGINNING OR/AND ***
    // *** ONE OR MORE SECONDARIES HAVE BEEN GENERATED ***

    // --- INITIAL PARTICLE TYPE HAS BEEN CHANGED ==> PUT NEW TYPE ON ---
    // --- THE GEANT TEMPORARY STACK ---

    // --- PUT PARTICLE ON THE STACK ---
    gkin[0] = result;
    gkin[0].SetTOF( result.GetTOF() * 5e-11 );
    ngkine = 1;

    // --- ALL QUANTITIES ARE TAKEN FROM THE GHEISHA STACK WHERE THE ---
    // --- CONVENTION IS THE FOLLOWING ---

    // --- ONE OR MORE SECONDARIES HAVE BEEN GENERATED ---
    for (l = 1; l <= ntot; ++l) {
      index = l - 1;
      // DHW  15 May 2011: unused: jnd = eve[index].GetParticleDef();

      // --- ADD PARTICLE TO THE STACK IF STACK NOT YET FULL ---
      if (ngkine < MAX_SECONDARIES) {
	gkin[ngkine] = eve[index];
	gkin[ngkine].SetTOF( eve[index].GetTOF() * 5e-11 );
	++ngkine;
      }
    }
  }
  else {
    // --- NO SECONDARIES GENERATED AND PARTICLE IS STILL THE SAME ---
    // --- ==> COPY EVERYTHING BACK IN THE CURRENT GEANT STACK ---
    ngkine = 0;
    ntot = 0;
    globalTime += result.GetTOF() * G4float(5e-11);
  }

  // --- LIMIT THE VALUE OF NGKINE IN CASE OF OVERFLOW ---
  ngkine = G4int(std::min(ngkine,G4int(MAX_SECONDARIES)));

} // GenerateSecondaries


void G4AntiNeutronAnnihilationAtRest::Poisso(G4float xav, G4int *iran)
{
  G4int i;
  G4float r, p1, p2, p3;
  G4int fivex;
  G4float rr, ran, rrr, ran1;

  // *** GENERATION OF POISSON DISTRIBUTION ***
  // *** NVE 16-MAR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (27-OCT-1983)

  // --- USE NORMAL DISTRIBUTION FOR <X> > 9.9 ---
  if (xav > G4float(9.9)) {
    // ** NORMAL DISTRIBUTION WITH SIGMA**2 = <X>
    Normal(&ran1);
    ran1 = xav + ran1 * std::sqrt(xav);
    *iran = G4int(ran1);
    if (*iran < 0) {
      *iran = 0;
    }
  }
  else {
    fivex = G4int(xav * G4float(5.));
    *iran = 0;
    if (fivex > 0) {
      r = G4Exp(-G4double(xav));
      ran1 = G4UniformRand();
      if (ran1 > r) {
	rr = r;
	for (i = 1; i <= fivex; ++i) {
	  ++(*iran);
	  if (i <= 5) {
	    rrr = G4Pow::GetInstance()->powN(xav, i) / NFac(i);
	  }
	  // ** STIRLING' S FORMULA FOR LARGE NUMBERS
	  if (i > 5) {
	    rrr = G4Exp(i * G4Log(xav) -
		      (i + G4float(.5)) * G4Log(i * G4float(1.)) +
		      i - G4float(.9189385));
	  }
	  rr += r * rrr;
	  if (ran1 <= rr) {
	    break;
	  }
	}
      }
    }
    else {
      // ** FOR VERY SMALL XAV TRY IRAN=1,2,3
      p1 = xav * G4Exp(-G4double(xav));
      p2 = xav * p1 / G4float(2.);
      p3 = xav * p2 / G4float(3.);
      ran = G4UniformRand();
      if (ran >= p3) {
	if (ran >= p2) {
	  if (ran >= p1) {
	    *iran = 0;
	  }
	  else {
	    *iran = 1;
	  }
	}
	else {
	  *iran = 2;
	}
      }
      else {
	*iran = 3;
      }
    }
  }

} // Poisso


G4int G4AntiNeutronAnnihilationAtRest::NFac(G4int n)
{
  G4int ret_val;

  G4int i, j;

  // *** NVE 16-MAR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (27-OCT-1983)

  ret_val = 1;
  j = n;
  if (j > 1) {
    if (j > 10) {
      j = 10;
    }
    for (i = 2; i <= j; ++i) {
      ret_val *= i;
    }
  }
  return ret_val;

} // NFac


void G4AntiNeutronAnnihilationAtRest::Normal(G4float *ran)
{
  // *** NVE 14-APR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (27-OCT-1983)

  *ran = (G4float)(-6. + 12.*G4UniformRand());
} // Normal


void G4AntiNeutronAnnihilationAtRest::AntiNeutronAnnihilation(G4int *nopt)
{
  G4float brr[3] = { G4float(.125),G4float(.25),G4float(.5) };

  G4float r__1;

  G4int i, ii, kk;
  G4int nt;
  G4float cfa, eka;
  G4int ika, nbl;
  G4float ran, pcm;
  G4int isw;
  G4float tex;
  G4ParticleDefinition* ipa1;
  G4float ran1, ran2, ekin, tkin;
  G4float targ;
  G4ParticleDefinition* inve;
  G4float ekin1, ekin2, black;
  G4float pnrat, rmnve1, rmnve2;
  G4float ek, en;

  // *** ANTI NEUTRON ANNIHILATION AT REST ***
  // *** NVE 04-MAR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (09-JULY-1987)

  // NOPT=0    NO ANNIHILATION
  // NOPT=1    ANNIH.IN PI+ PI-
  // NOPT=2    ANNIH.IN PI0 PI0
  // NOPT=3    ANNIH.IN PI+ PI0
  // NOPT=4    ANNIH.IN GAMMA GAMMA

  pv[1].SetZero();
  pv[1].SetMass( massAntiNeutron );
  pv[1].SetKineticEnergyAndUpdate( 0. );
  pv[1].SetTOF( result.GetTOF() );
  pv[1].SetParticleDef( result.GetParticleDef() );
  isw = 1;
  ran = G4UniformRand();
  if (ran > brr[0]) {
    isw = 2;
  }
  if (ran > brr[1]) {
    isw = 3;
  }
  if (ran > brr[2]) {
    isw = 4;
  }
  *nopt = isw;
  // **
  // **  EVAPORATION
  // **
  rmnve1 = massPionPlus;
  rmnve2 = massPionMinus;
  if (isw == 2) {
    rmnve1 = massPionZero;
    rmnve2 = massPionZero;
  }
  if (isw == 3) {
    rmnve2 = massPionZero;
  }
  if (isw == 4) {
    rmnve1 = massGamma;
    rmnve2 = massGamma;
  }
  ek = massNeutron + massAntiNeutron - rmnve1 - rmnve2;
  tkin = ExNu(ek);
  ek -= tkin;
  if (ek < G4float(1e-4)) {
    ek = G4float(1e-4);
  }
  ek /= G4float(2.);
  en = ek + (rmnve1 + rmnve2) / G4float(2.);
  r__1 = en * en - rmnve1 * rmnve2;
  pcm = r__1 > 0 ? std::sqrt(r__1) : 0;
  pv[2].SetZero();
  pv[2].SetMass( rmnve1 );
  pv[3].SetZero();
  pv[3].SetMass( rmnve2 );
  if (isw > 3) {
    pv[2].SetMass( 0. );
    pv[3].SetMass( 0. );
  }
  pv[2].SetEnergyAndUpdate( std::sqrt(pv[2].GetMass()*pv[2].GetMass()+pcm*pcm) );
  pv[2].SetTOF( result.GetTOF() );
  pv[3].SetEnergy( std::sqrt(pv[3].GetMass()*pv[3].GetMass()+pcm*pcm) );
  pv[3].SetMomentumAndUpdate( -pv[2].GetMomentum().x(), -pv[2].GetMomentum().y(), -pv[2].GetMomentum().z() );
  pv[3].SetTOF( result.GetTOF() );
  switch ((int)isw) {
    case 1:
      pv[2].SetParticleDef( pdefPionPlus );
      pv[3].SetParticleDef( pdefPionMinus );
      break;
    case 2:
      pv[2].SetParticleDef( pdefPionZero );
      pv[3].SetParticleDef( pdefPionZero );
      break;
    case 3:
      pv[2].SetParticleDef( pdefPionPlus );
      pv[3].SetParticleDef( pdefPionZero );
      break;
    case 4:
      pv[2].SetParticleDef( pdefGamma );
      pv[3].SetParticleDef( pdefGamma );
      break;
  }
  nt = 3;
  if (targetAtomicMass >= G4float(1.5)) {
    cfa = (targetAtomicMass - G4float(1.)) / G4float(120.) *
      G4float(.025) * G4Exp(-G4double(targetAtomicMass - G4float(1.)) /
			  G4float(120.));
    targ = G4float(1.);
    tex = evapEnergy1;
    if (tex >= G4float(.001)) {
      black = (targ * G4float(1.25) +
	       G4float(1.5)) * evapEnergy1 / (evapEnergy1 + evapEnergy3);
      Poisso(black, &nbl);
      if (G4float(G4int(targ) + nbl) > targetAtomicMass) {
	nbl = G4int(targetAtomicMass - targ);
      }
      if (nt + nbl > (MAX_SECONDARIES - 2)) {
	nbl = (MAX_SECONDARIES - 2) - nt;
      }
      if (nbl > 0) {
	ekin = tex / nbl;
	ekin2 = 0.0f;
	for (i = 1; i <= nbl; ++i) {
	  if (nt == (MAX_SECONDARIES - 2)) {
	    continue;
	  }
	  if (ekin2 > tex) {
	    break;
	  }
	  ran1 = G4UniformRand();
	  Normal(&ran2);
	  ekin1 = -G4double(ekin) * G4Log(ran1) -
	    cfa * (ran2 * G4float(.5) + G4float(1.));
	  if (ekin1 < 0.0f) {
	    ekin1 = G4Log(ran1) * G4float(-.01);
	  }
	  ekin1 *= G4float(1.);
	  ekin2 += ekin1;
	  if (ekin2 > tex) {
	    ekin1 = tex - (ekin2 - ekin1);
	  }
	  if (ekin1 < 0.0f) {
	    ekin1 = G4float(.001);
	  }
	  ipa1 = pdefNeutron;
	  pnrat = G4float(1.) - targetCharge / targetAtomicMass;
	  if (G4UniformRand() > pnrat) {
	    ipa1 = pdefProton;
	  }
	  ++nt;
	  pv[nt].SetZero();
	  pv[nt].SetMass( ipa1->GetPDGMass()/GeV );
	  pv[nt].SetKineticEnergyAndUpdate( ekin1 );
	  pv[nt].SetTOF( result.GetTOF() );
	  pv[nt].SetParticleDef( ipa1 );
	}
	if (targetAtomicMass >= G4float(230.) && ek <= G4float(2.)) {
	  ii = nt + 1;
	  kk = 0;
	  eka = ek;
	  if (eka > G4float(1.)) {
	    eka *= eka;
	  }
	  if (eka < 0.1f) {
	    eka = 0.1f;
	  }
	  ika = G4int(G4float(3.6) / eka);
	  for (i = 1; i <= nt; ++i) {
	    --ii;
	    if (pv[ii].GetParticleDef() != pdefProton) {
	      continue;
	    }
	    ipa1 = pdefNeutron;
	    pv[ii].SetMass( ipa1->GetPDGMass()/GeV );
	    pv[ii].SetParticleDef( ipa1 );
	    ++kk;
	    if (kk > ika) {
	      break;
	    }
	  }
	}
      }
    }
    // **
    // ** THEN ALSO DEUTERONS, TRITONS AND ALPHAS
    // **
    tex = evapEnergy3;
    if (tex >= G4float(.001)) {
      black = (targ * G4float(1.25) + G4float(1.5)) * evapEnergy3 /
	(evapEnergy1 + evapEnergy3);
      Poisso(black, &nbl);
      if (nt + nbl > (MAX_SECONDARIES - 2)) {
	nbl = (MAX_SECONDARIES - 2) - nt;
      }
      if (nbl > 0) {
	ekin = tex / nbl;
	ekin2 = 0.0f;
	for (i = 1; i <= nbl; ++i) {
	  if (nt == (MAX_SECONDARIES - 2)) {
	    continue;
	  }
	  if (ekin2 > tex) {
	    break;
	  }
	  ran1 = G4UniformRand();
	  Normal(&ran2);
	  ekin1 = -G4double(ekin) * G4Log(ran1) -
	    cfa * (ran2 * G4float(.5) + G4float(1.));
	  if (ekin1 < 0.0f) {
	    ekin1 = G4Log(ran1) * G4float(-.01);
	  }
	  ekin1 *= G4float(1.);
	  ekin2 += ekin1;
	  if (ekin2 > tex) {
	    ekin1 = tex - (ekin2 - ekin1);
	  }
	  if (ekin1 < 0.0f) {
	    ekin1 = G4float(.001);
	  }
	  ran = G4UniformRand();
	  inve = pdefDeuteron;
	  if (ran > G4float(.6)) {
	    inve = pdefTriton;
	  }
	  if (ran > G4float(.9)) {
	    inve = pdefAlpha;
	  }
	  ++nt;
	  pv[nt].SetZero();
	  pv[nt].SetMass( inve->GetPDGMass()/GeV );
	  pv[nt].SetKineticEnergyAndUpdate( ekin1 );
	  pv[nt].SetTOF( result.GetTOF() );
	  pv[nt].SetParticleDef( inve );
	}
      }
    }
  }
  result = pv[2];
  if (nt == 2) {
    return;
  }
  for (i = 3; i <= nt; ++i) {
    if (ntot >= MAX_SECONDARIES) {
      return;
    }
    eve[ntot++] = pv[i];
  }

} // AntiNeutronAnnihilation


G4double G4AntiNeutronAnnihilationAtRest::ExNu(G4float ek1)
{
  G4float ret_val, r__1;

  G4float cfa, gfa, ran1, ran2, ekin1, atno3;
  G4int magic;
  G4float fpdiv;

  // *** NUCLEAR EVAPORATION AS FUNCTION OF ATOMIC NUMBER ATNO ***
  // *** AND KINETIC ENERGY EKIN OF PRIMARY PARTICLE ***
  // *** NVE 04-MAR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (10-DEC-1986)

  ret_val = 0.f;
  if (targetAtomicMass >= G4float(1.5)) {
    magic = 0;
    if (G4int(targetCharge + 0.1f) == 82) {
      magic = 1;
    }
    ekin1 = ek1;
    if (ekin1 < 0.1f) {
      ekin1 = 0.1f;
    }
    if (ekin1 > G4float(4.)) {
      ekin1 = G4float(4.);
    }
    // **   0.35 VALUE AT 1 GEV
    // **   0.05 VALUE AT 0.1 GEV
    cfa = G4float(.13043478260869565);
    cfa = cfa * G4Log(ekin1) + G4float(.35);
    if (cfa < G4float(.15)) {
      cfa = G4float(.15);
    }
    ret_val = cfa * G4float(7.716) * G4Exp(-G4double(cfa));
    atno3 = targetAtomicMass;
    if (atno3 > G4float(120.)) {
      atno3 = G4float(120.);
    }
    cfa = (atno3 - G4float(1.)) /
      G4float(120.) * G4Exp(-G4double(atno3 - G4float(1.)) / G4float(120.));
    ret_val *= cfa;
    r__1 = ekin1;
    fpdiv = G4float(1.) - r__1 * r__1 * G4float(.25);
    if (fpdiv < G4float(.5)) {
      fpdiv = G4float(.5);
    }
    gfa = (targetAtomicMass - G4float(1.)) /
      G4float(70.) * G4float(2.) *
      G4Exp(-G4double(targetAtomicMass - G4float(1.)) / G4float(70.));
    evapEnergy1 = ret_val * fpdiv;
    evapEnergy3 = ret_val - evapEnergy1;
    Normal(&ran1);
    Normal(&ran2);
    if (magic == 1) {
      ran1 = 0.0f;
      ran2 = 0.0f;
    }
    evapEnergy1 *= ran1 * gfa + G4float(1.);
    if (evapEnergy1 < 0.0f) {
      evapEnergy1 = 0.0f;
    }
    evapEnergy3 *= ran2 * gfa + G4float(1.);
    if (evapEnergy3 < 0.0f) {
      evapEnergy3 = 0.0f;
    }

    // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
    while ((ret_val = evapEnergy1 + evapEnergy3) >= ek1) {
      evapEnergy1 *= G4float(1.) - G4UniformRand() * G4float(.5);
      evapEnergy3 *= G4float(1.) - G4UniformRand() * G4float(.5);
    }
  }
  return ret_val;

} // ExNu
