// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronCaptureAtRest.cc,v 1.2 1999-12-15 14:53:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4NeutronCaptureAtRest physics process --------
//                   by Larry Felawka (TRIUMF), April 1998
//                     E-mail: felawka@alph04.triumf.ca
// **************************************************************
//-----------------------------------------------------------------------------

#include "G4NeutronCaptureAtRest.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh" 
#include <string.h>
#include <math.h>
#include <stdio.h>
 
#define MAX_SECONDARIES 100

// constructor
 
G4NeutronCaptureAtRest::G4NeutronCaptureAtRest(const G4String& processName)
  : G4VRestProcess (processName),       // initialization
  massProton(G4Proton::Proton()->GetPDGMass()/GeV),
  massNeutron(G4Neutron::Neutron()->GetPDGMass()/GeV),
  massElectron(G4Electron::Electron()->GetPDGMass()/GeV),
  massDeuteron(G4Deuteron::Deuteron()->GetPDGMass()/GeV),
  massAlpha(G4Alpha::Alpha()->GetPDGMass()/GeV),
  pdefGamma(G4Gamma::Gamma()),
  pdefNeutron(G4Neutron::Neutron())
{
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }

  pv   = new G4GHEKinematicsVector [MAX_SECONDARIES+1];
  eve  = new G4GHEKinematicsVector [MAX_SECONDARIES];
  gkin = new G4GHEKinematicsVector [MAX_SECONDARIES];

}
 
// destructor
 
G4NeutronCaptureAtRest::~G4NeutronCaptureAtRest(){;}
 
 
// methods.............................................................................
 
G4bool G4NeutronCaptureAtRest::IsApplicable(
				 const G4ParticleDefinition& particle
				 )
{
   return ( &particle == pdefNeutron );

}
 
// Warning - this method may be optimized away if made "inline"
G4int G4NeutronCaptureAtRest::GetNumberOfSecondaries()
{
  return ( ngkine );

}

// Warning - this method may be optimized away if made "inline"
G4GHEKinematicsVector* G4NeutronCaptureAtRest::GetSecondaryKinematics()
{
  return ( &gkin[0] );

}

G4double G4NeutronCaptureAtRest::AtRestGetPhysicalInteractionLength(
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
    G4cout << "G4NeutronCaptureAtRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }

  return theNumberOfInteractionLengthLeft * currentInteractionLength;

}

G4VParticleChange* G4NeutronCaptureAtRest::AtRestDoIt(
                                            const G4Track& track,
					    const G4Step& stepData
					    )
//
// Handles Neutrons at rest; a Neutron can either create secondaries or
// do nothing (in which case it should be sent back to decay-handling
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
      targetCharge = G4double((*theElementVector)(i2)->GetZ());
      targetAtomicMass = (*theElementVector)(i2)->GetN();
    }
  }
  if (random>runningSum)
  {
    targetCharge = G4double((*theElementVector)(numberOfElements-1)->GetZ());
    targetAtomicMass = (*theElementVector)(numberOfElements-1)->GetN();

  }

  if (verboseLevel>1) {
    G4cout << "G4NeutronCaptureAtRest::AtRestDoIt is invoked " <<G4endl;
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
    aParticleChange.AddSecondary( aNewTrack );

  }

  aParticleChange.SetLocalEnergyDeposit( 0.0*GeV );

  aParticleChange.SetStatusChange(fStopAndKill); // Kill the incident Neutron

//   clear InteractionLengthLeft

  ResetNumberOfInteractionLengthLeft();

  return &aParticleChange;

}


void G4NeutronCaptureAtRest::GenerateSecondaries()
{
  static G4int index;
  static G4int l;
  static G4int nopt;
  static G4int i;
  static G4ParticleDefinition* jnd;

  for (i = 1; i <= MAX_SECONDARIES; ++i) {
    pv[i].SetZero();
  }

  ngkine = 0;            // number of generated secondary particles
  ntot = 0;
  result.SetZero();
  result.SetMass( massNeutron );
  result.SetKineticEnergyAndUpdate( 0. );
  result.SetTOF( 0. );
  result.SetParticleDef( pdefNeutron );

  NeutronCapture(&nopt);

  // *** CHECK WHETHER THERE ARE NEW PARTICLES GENERATED ***
  if (ntot != 0 || result.GetParticleDef() != pdefNeutron) {
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
      jnd = eve[index].GetParticleDef();

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
  ngkine = G4int(G4std::min(ngkine,G4int(MAX_SECONDARIES)));

} // GenerateSecondaries


void G4NeutronCaptureAtRest::Normal(G4float *ran)
{
  static G4int i;

  // *** NVE 14-APR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (27-OCT-1983)

  *ran = G4float(-6.);
  for (i = 1; i <= 12; ++i) {
    *ran += G4UniformRand();
  }

} // Normal


void G4NeutronCaptureAtRest::NeutronCapture(G4int *nopt)
{
  static G4int nt;
  static G4float xp, pcm;
  static G4float ran;

  // *** ROUTINE FOR CAPTURE OF NEUTRAL BARYONS ***
  // *** NVE 04-MAR-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (02-DEC-1986)

  *nopt = 1;
  pv[1] = result;
  pv[2].SetZero();
  pv[2].SetMass( AtomAs(targetAtomicMass, targetCharge) );
  pv[2].SetMomentumAndUpdate( 0., 0., 0. );
  pv[2].SetTOF( result.GetTOF() );
  pv[2].SetParticleDef( NULL );
  pv[MAX_SECONDARIES].Add( pv[1], pv[2] );
  pv[MAX_SECONDARIES].SetMomentum( -pv[MAX_SECONDARIES].GetMomentum().x(), -pv[MAX_SECONDARIES].GetMomentum().y(), -pv[MAX_SECONDARIES].GetMomentum().z() );
  pv[MAX_SECONDARIES].SetParticleDef( NULL );
  Normal(&ran);
  pcm = ran * G4float(.001) + G4float(.0065);
  ran = G4UniformRand();
  result.SetTOF( result.GetTOF() - log(ran) * G4float(480.) );
  pv[3].SetZero();
  pv[3].SetMass( 0. );
  pv[3].SetKineticEnergyAndUpdate( pcm );
  pv[3].SetTOF( result.GetTOF() );
  pv[3].SetParticleDef( pdefGamma );
  pv[3].Lor( pv[3], pv[MAX_SECONDARIES] );
  nt = 3;
  xp = G4float(.008) - pcm;
  if (xp >= G4float(0.)) {
    nt = 4;
    pv[4].SetZero();
    pv[4].SetMass( 0. );
    pv[4].SetKineticEnergyAndUpdate( xp );
    pv[4].SetTOF( result.GetTOF() );
    pv[4].SetParticleDef( pdefGamma );
    pv[4].Lor( pv[4], pv[MAX_SECONDARIES] );
  }
  result = pv[3];
  if (nt == 4) {
    if (ntot < MAX_SECONDARIES-1) {
      eve[ntot++] = pv[4];
    }
  }

} // NeutronCapture


G4double G4NeutronCaptureAtRest::AtomAs(G4float a, G4float z)
{
  G4float ret_val;
  G4double d__1, d__2;

  static G4double aa;
  static G4int ia, iz;
  static G4double zz;
  static G4float rma, rmd;
  static G4int ipp;
  static G4float rmn, rmp;
  static G4int izz;
  static G4float rmel;
  static G4double mass;

  // *** DETERMINATION OF THE ATOMIC MASS ***
  // *** NVE 19-MAY-1988 CERN GENEVA ***
  // ORIGIN : H.FESEFELDT (02-DEC-1986)

  // --- GET ATOMIC (= ELECTRONS INCL.) MASSES (IN MEV) FROM RMASS ARRAY ---
  // --- ELECTRON ---
  rmel = massElectron * G4float(1e3);
  // --- PROTON ---
  rmp = massProton * G4float(1e3);
  // --- NEUTRON ---
  rmn = massNeutron * G4float(1e3);
  // --- DEUTERON ---
  rmd = massDeuteron * G4float(1e3) + rmel;
  // --- ALPHA ---
  rma = massAlpha * G4float(1e3) + rmel * G4float(2.);

  ret_val = G4float(0.);
  aa = a * 1.;
  zz = z * 1.;
  ia = G4int(a + G4float(.5));
  if (ia < 1) {
    return ret_val;
  }
  iz = G4int(z + G4float(.5));
  if (iz < 0 || iz > ia) {
    return ret_val;
  }
  mass = 0.;
  if (ia == 1) {
    if (iz == 0) {
      mass = rmn;
    }
    else if (iz == 1) {
      mass = rmp + rmel;
    }
  }
  else if (ia == 2 && iz == 1) {
    mass = rmd;
  }
  else if (ia == 4 && iz == 2) {
    mass = rma;
  }
  else if (ia == 2 && iz != 1 || ia == 3 || ia == 4 && iz != 2 || ia > 4) {
    d__1 = aa / G4float(2.) - zz;
    d__2 = zz;
    mass = (aa - zz) * rmn + zz * rmp + zz * rmel - aa * G4float(15.67) +
      pow(aa, .6666667) * G4float(17.23) + d__1 * d__1 * G4float(93.15) / aa +
      d__2 * d__2 * G4float(.6984523) / pow(aa, .3333333);
    ipp = (ia - iz) % 2;
    izz = iz % 2;
    if (ipp == izz) {
      mass += (ipp + izz - 1) * G4float(12.) * pow(aa, -.5);
    }
  }
  ret_val = mass * G4float(.001);
  return ret_val;

} // AtomAs
