// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuonMinusCaptureAtRest.cc,v 1.4 2000-07-12 09:17:11 hpw Exp $
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
//      ------------ G4MuonMinusCaptureAtRest physics process --------
//                   by Larry Felawka (TRIUMF)
//                     E-mail: felawka@alph04.triumf.ca
//                   and Art Olin (TRIUMF)
//                     E-mail: olin@triumf.ca
//                            April 1998
// **************************************************************
//      V.Ivanchenko   7 Apr 2000 Advance model for electromagnetic
//                                capture and cascade
//-----------------------------------------------------------------------------

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "G4He3.hh"


#define NINT(x) ((x)>=0 ? (G4int)floor((x) + .5) : -(G4int)floor(.5 - (x)))

#define Y0  1.5
#define B0  8.
#define AMPROT .93827231
#define AMNEUT .93956563
#define AMMUON .105658389	// Muon mass (GeV)
#define AMMRED ( AMNEUT / ( AMPROT + AMMUON ) )
#define FERTHO .00000001433
#define BEXC12 ( FERTHO * 72.40715579499394 )
#define AMUAMU .93149432
#define AMELCT 5.1099906e-4	// Electron mass (GeV)
#define AMELEC AMELCT
#define AMUC12 ( AMUAMU - AMELEC / 2. - BEXC12 / 12. )
#define EBNDAV ( 0.5 * (AMPROT + AMNEUT) - AMUC12 )
#define R0NUCL 1.12
#define PLABRC .197327053 // Reduced Planck constant Times the light velocity
#define FSCTO2 5.3251361962113614e-5 // (Fine structure constant)^2
#define O16OLD 931.145
#define O16NEW 931.19826
#define O16RAT ( O16NEW / O16OLD )
#define C12NEW 931.49432
#define ADJUST -.08322737768178909
#define C0M1E1 .306
#define C0E2E1 .71
#define HNDFM1 .05
#define HNDFE2 10.
//  Gammin : threshold for deexcitation gammas production, set to 1 keV
//  (this means that up to 1 keV of energy unbalancing can occur
//   during an event)
#define GAMMIN 1.0e-6
//  Tvepsi : "epsilon" for excitation energy, set to gammin / 100
#define TVEPSI ( GAMMIN / 100. )
//  Anglgb = this parameter should be set equal to the machine
//           "zero" with respect to unit
#define ANGLGB 5.0e-16
//  Onepls = 1+ of the machine, it is 1 + 2 x Anglgb
#define ONEPLS 1.000000000000001
#define ASMTOG ( 6.0 / M_PI*M_PI )
#define COULPR 0.001439965
#define AMUMEV ( 1.0e3 * AMUAMU )
#define HNDFE1 5.e-03
#define RCCOUL 1.7
#define EXPMIN -88.
#define EXPMAX 88.
#define IDMAX8 183
#define MXFRAG 100
#define MXGKIN 500
#define MXSECS 999
#define MXEVAP 100
#define NALLWP 39


// constructor

G4MuonMinusCaptureAtRest::G4MuonMinusCaptureAtRest(const G4String& processName)
  : G4VRestProcess (processName),       // initialization
    massGamma(G4Gamma::Gamma()->GetPDGMass()/GeV),
    massNeutrinoE(G4NeutrinoE::NeutrinoE()->GetPDGMass()/GeV),
    massProton(G4Proton::Proton()->GetPDGMass()/GeV),
    massNeutron(G4Neutron::Neutron()->GetPDGMass()/GeV),
    pdefGamma(G4Gamma::Gamma()),
    pdefNeutrinoE(G4NeutrinoE::NeutrinoE()),
    pdefMuonMinus(G4MuonMinus::MuonMinus()),
    pdefNeutron(G4Neutron::Neutron())
{
  massFragm[0] = G4Neutron::Neutron()->GetPDGMass()/GeV;
  massFragm[1] = G4Proton::Proton()->GetPDGMass()/GeV;
  massFragm[2] = G4Deuteron::Deuteron()->GetPDGMass()/GeV;
  massFragm[3] = G4Triton::Triton()->GetPDGMass()/GeV;
  massFragm[4] = G4He3::He3()->GetPDGMass()/GeV;
  massFragm[5] = G4Alpha::Alpha()->GetPDGMass()/GeV;

  pdefFragm[0] = G4Neutron::Neutron();
  pdefFragm[1] = G4Proton::Proton();
  pdefFragm[2] = G4Deuteron::Deuteron();
  pdefFragm[3] = G4Triton::Triton();
  pdefFragm[4] = G4He3::He3();
  pdefFragm[5] = G4Alpha::Alpha();

  nGkine = 0;

  if (verboseLevel>=0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }

  Fragments   = new G4GHEKinematicsVector [MXFRAG];
  Secondaries = new G4GHEKinematicsVector [MXSECS];
  Evaporates  = new G4GHEKinematicsVector [MXEVAP*6];
  Gkin        = new G4GHEKinematicsVector [MXGKIN];
  Cascade     = new G4GHEKinematicsVector [17];

  InitializeMuCapture();

  pSelector  = new G4StopElementSelector();
  pEMCascade = new G4MuMinusCaptureCascade();
}

// destructor

G4MuonMinusCaptureAtRest::~G4MuonMinusCaptureAtRest()
{
  delete [] Fragments;
  delete [] Secondaries;
  delete [] Evaporates;
  delete [] Gkin;
  delete [] Cascade;

}


// methods.............................................................................

G4bool G4MuonMinusCaptureAtRest::IsApplicable(
					      const G4ParticleDefinition& particle
					      )
{
  return ( &particle == pdefMuonMinus );

}

G4double G4MuonMinusCaptureAtRest::GetMeanLifeTime(const G4Track& track,
						   G4ForceCondition* condition)
{
  GetCaptureIsotope( track );
  return (tDelay );

}

// Warning - this method may be optimized away if made "inline"
G4int G4MuonMinusCaptureAtRest::GetNumberOfSecondaries()
{
  return ( nGkine );

}

// Warning - this method may be optimized away if made "inline"
G4GHEKinematicsVector* G4MuonMinusCaptureAtRest::GetSecondaryKinematics()
{
  return ( &Gkin[0] );

}

G4double G4MuonMinusCaptureAtRest::AtRestGetPhysicalInteractionLength(
								      const G4Track& track,
								      G4ForceCondition* condition
								      )
{
  // beggining of tracking
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  //  *condition = NotForced;
  // condition is set to "ExclusivelyForced" by V.Ivanchenko
       *condition = ExclusivelyForced;

  // get mean life time
  currentInteractionLength = GetMeanLifeTime(track, condition);

  if ((currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4MuonMinusCaptureAtRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }

  // Return 0 interaction length to get for this process
  // the 100% probability 

  return 0.0;
  //  return theNumberOfInteractionLengthLeft * currentInteractionLength;

}

G4VParticleChange* G4MuonMinusCaptureAtRest::AtRestDoIt(
							const G4Track& track,
							const G4Step& stepData
							)
//
// Handles MuonMinuss at rest; a MuonMinus can either create secondaries or
// do nothing (in which case it should be sent back to decay-handling
// section
//
{

//   Initialize ParticleChange
//     all members of G4VParticleChange are set to equal to
//     corresponding member in G4Track

  aParticleChange.Initialize(track);

//   Store some global quantities that depend on current material and particle

  G4double globalTime = track.GetGlobalTime()/s;

  if (verboseLevel > 1) {
    G4cout << "G4MuonMinusCaptureAtRest::AtRestDoIt is invoked " << G4endl;
  }

  G4ParticleMomentum momentum;
  G4double localtime;

  G4ThreeVector   position = track.GetPosition();

  // Generate secondaries from electromagnetic cascade
  nCascade = 0;
  G4double mass = GetIsotopicMass(targetAtomicMass,targetCharge) * GeV;
  nCascade      = pEMCascade->DoCascade(targetCharge, mass, Cascade);

  // Decay or Capture?
  G4double lambdac  = pSelector->GetMuonCaptureRate(targetCharge,targetAtomicMass);
  G4double lambdad  = pSelector->GetMuonDecayRate(targetCharge,targetAtomicMass);
  
  if( G4UniformRand()*(lambdac + lambdad) > lambdac) {

    // Decay
    pEMCascade->DoBoundMuonMinusDecay(targetCharge,mass,&nCascade,Cascade);

    // Generate secondaries from capture
  } else {

    DoMuCapture();
  }
    
  aParticleChange.SetNumberOfSecondaries( nGkine + nCascade );

  // Store nuclear cascade
  if(nGkine > 0) {
    for ( G4int isec = 0; isec < nGkine; isec++ ) {
      G4DynamicParticle* aNewParticle = new G4DynamicParticle;
      aNewParticle->SetDefinition( Gkin[isec].GetParticleDef() );
      aNewParticle->SetMomentum( Gkin[isec].GetMomentum() * GeV );

      localtime = globalTime + Gkin[isec].GetTOF();

      G4Track* aNewTrack = new G4Track( aNewParticle, localtime*s, position );
      aParticleChange.AddSecondary( aNewTrack );
    }
  }

  // Store electromagnetic cascade

  if(nCascade > 0) {
    localtime = globalTime*s + tDelay;

    for ( G4int isec = 0; isec < nCascade; isec++ ) {
      G4DynamicParticle* aNewParticle = new G4DynamicParticle;
      aNewParticle->SetDefinition( Cascade[isec].GetParticleDef() );
      aNewParticle->SetMomentum( Cascade[isec].GetMomentum() );

      G4Track* aNewTrack = new G4Track( aNewParticle, localtime, position );
      aParticleChange.AddSecondary( aNewTrack );
    }
  }

  aParticleChange.SetLocalEnergyDeposit(0.0);

  aParticleChange.SetStatusChange(fStopAndKill); // Kill the incident MuonMinus

//   clear InteractionLengthLeft

  ResetNumberOfInteractionLengthLeft();

  return &aParticleChange;

}


void G4MuonMinusCaptureAtRest::CascadeCorrection(G4double zztar,
						 G4double bbtar)
{
  // System generated locals
  G4double d__1;

  // Local variables
  static G4int i;
  static G4double hkap;
  static G4double p2help;
  static G4double efrmmx, ato1o3;
  static G4double hhlp[2];

  const G4double onethird = 1./3.;
  const G4double apfrmx = pow(9.*M_PI/8.,onethird) * PLABRC / R0NUCL;
  static G4double bbold = -1e10;
  static G4double zzold = -1e10;

  //  Reduction factors for intran. cascade energy, taken from Alsmiller
  //  Incoming baryons
  //  Incoming mesons

  if (bbtar != bbold || zztar != zzold) {
    //  Supply the fraction of the total kinetic energy to be
    //  used for intranuclear cascade nucleons
    ato1o3 = pow(bbtar, onethird);
    d__1 = bbtar - zztar;
    hkap = bbtar * bbtar / (zztar * zztar + d__1 * d__1);
    hhlp[0] = pow(hkap * zztar, onethird) / ato1o3;
    hhlp[1] = pow(hkap * (bbtar - zztar), onethird) / ato1o3;
    for (i = 1; i <= 2; ++i) {
      nucleonMaxFermiEn[i - 1] = hhlp[i - 1] * apfrmx;
      d__1 = nucleonMaxFermiEn[i - 1];
      p2help = d__1 * d__1;
      efrmmx = sqrt(p2help + nucleonMassSquared[i - 1]) - nucleonMass[i - 1];
      nucleonPotWell[i - 1] = efrmmx + nucleonBindingEn[i - 1];
    }
    bbold = bbtar;
    zzold = zztar;
  }
  return;

} // CascadeCorrection


G4double G4MuonMinusCaptureAtRest::CoulombBarrier(G4int i, G4double z)
{
  // System generated locals
  G4double ret_val;

  // Local variables
  static G4int n;
  static G4double x;

  static G4double t[28] = {
     0.36, 0.77, 0.08, 0.00, 0.51, 0.81, 0.00,
     0.00, 0.60, 0.85,-0.06, 0.00, 0.66, 0.89,
    -0.10, 0.00, 0.68, 0.93,-0.10, 0.00, 0.69,
     0.97,-0.10, 0.00, 0.69, 1.00,-0.10, 0.00
  };

  if (z >= 70.) {
    ret_val = t[i + 23];
  } else {
    if (z <= 10.) {
      ret_val = t[i - 1];
    } else {
      n = (G4int) (z * .1 + 1.);
      x = (G4double) (n * 10);
      x = (x - z) * .1;
      ret_val = x * t[i + (n - 1 << 2) - 5] + (1. - x) * t[i + (n << 2) - 5];
    }
  }
  return ret_val;

} // CoulombBarrier


void G4MuonMinusCaptureAtRest::ExcitationEnergyLevel(G4int ja, G4int jz,
						     G4double *eex1st,
						     G4double *eex2nd,
						     G4double *eexcon)
{
  // Initialized data

  static G4int jaold = 0;

  // Local variables
  static G4double dja;
  static G4int iodd, inodd, izodd;
  static G4double sqatar;

  static G4double cam4[130] = {
    0.00, 5.44, 0.00, 2.76, 0.00, 3.34, 0.00, 2.70, 0.00, 1.90,
    0.00, 2.12, 0.00, 2.13, 0.00, 1.54, 0.00, 1.42, 0.00, 1.51,
    0.00, 1.73, 0.00, 1.44, 0.00, 1.45, 0.00, 1.37, 0.00, 1.09,
    0.00, 1.36, 0.00, 1.42, 0.00, 1.33, 0.00, 1.20, 0.00, 1.00,
    0.00, 1.16, 0.00, 1.28, 0.00, 1.38, 0.00, 1.38, 0.00, 1.32,
    0.00, 1.04, 0.00, 1.11, 0.00, 1.13, 0.00, 1.21, 0.00, 1.43,
    0.00, 1.15, 0.00, 0.99, 0.00, 0.91, 0.00, 0.92, 0.00, 1.00,
    0.00, 1.11, 0.00, 1.23, 0.00, 0.85, 0.00, 0.98, 0.00, 0.72,
    0.00, 0.80, 0.00, 0.77, 0.00, 0.89, 0.00, 0.92, 0.00, 0.80,
    0.00, 0.81, 0.00, 0.69, 0.00, 0.70, 0.00, 0.76, 0.00, 0.73,
    0.00, 0.80, 0.00, 0.74, 0.00, 0.73, 0.00, 0.72, 0.00, 0.72,
    0.00, 0.72, 0.00, 0.71, 0.00, 0.69, 0.00, 0.68, 0.00, 0.66,
    0.00, 0.61, 0.00, 0.42, 0.00, 0.36, 0.00, 0.41, 0.00, 0.49
  };
  static G4double cam5[200] = {
    0.00, 5.98, 0.00, 2.77, 0.00, 3.16, 0.00, 3.01, 0.00, 1.68,
    0.00, 1.73, 0.00, 2.17, 0.00, 1.74, 0.00, 1.75, 0.00, 1.72,
    0.00, 1.63, 0.00, 1.41, 0.00, 1.29, 0.00, 1.47, 0.00, 1.32,
    0.00, 1.46, 0.00, 1.44, 0.00, 1.46, 0.00, 1.52, 0.00, 1.51,
    0.00, 1.47, 0.00, 1.45, 0.00, 1.28, 0.00, 1.23, 0.00, 1.27,
    0.00, 0.62, 0.00, 0.76, 0.00, 1.23, 0.00, 1.22, 0.00, 1.40,
    0.00, 1.36, 0.00, 1.30, 0.00, 1.29, 0.00, 1.24, 0.00, 1.28,
    0.00, 1.24, 0.00, 1.20, 0.00, 0.94, 0.00, 1.00, 0.00, 1.05,
    0.00, 0.54, 0.00, 0.60, 0.00, 0.75, 0.00, 0.75, 0.00, 0.85,
    0.00, 0.97, 0.00, 1.02, 0.00, 1.05, 0.00, 1.06, 0.00, 1.07,
    0.00, 1.06, 0.00, 1.05, 0.00, 1.02, 0.00, 0.97, 0.00, 0.91,
    0.00, 0.83, 0.00, 0.74, 0.00, 0.66, 0.00, 0.61, 0.00, 0.61,
    0.00, 0.90, 0.00, 0.52, 0.00, 0.81, 0.00, 0.68, 0.00, 0.72,
    0.00, 0.77, 0.00, 0.68, 0.00, 0.67, 0.00, 0.80, 0.00, 0.68,
    0.00, 0.64, 0.00, 0.58, 0.00, 0.55, 0.00, 0.57, 0.00, 0.57,
    0.00, 0.55, 0.00, 0.60, 0.00, 0.58, 0.00, 0.58, 0.00, 0.61,
    0.00, 0.63, 0.00, 0.65, 0.00, 0.66, 0.00, 0.65, 0.00, 0.65,
    0.00, 0.64, 0.00, 0.64, 0.00, 0.63, 0.00, 0.61, 0.00, 0.59,
    0.00, 0.55, 0.00, 0.39, 0.00, 0.36, 0.00, 0.40, 0.00, 0.40,
    0.00, 0.40, 0.00, 0.40, 0.00, 0.40, 0.00, 0.40, 0.00, 0.40
  };

  // ---------------------------------------------------------------------*
  //                                                                      *
  //     Created on 06 december 1991  by    Alfredo Ferrari & Paola Sala  *
  //                                                   Infn - Milan       *
  //                                                                      *
  //     Last change on 28-apr-92     by    Alfredo Ferrari               *
  //                                                                      *
  //                                                                      *
  // ---------------------------------------------------------------------*

  if (ja == jz) {
    *eexcon = 0.;
    *eex1st = 0.;
    *eex2nd = 0.;
    return;
  }
  if (jz == 0 || ja == jz) {
    *eexcon = 0.;
  } else {
    *eexcon = (cam4[jz - 1] + cam5[ja - jz - 1]) * .001;
  }
  //  **** Very tentative selection of Eex1st, Eex2nd, based on pairing
  //       energies according to delta = 12 MeV / A^1/2 ****
  izodd = 1 - jz % 2;
  inodd = 1 - (ja - jz) % 2;
  iodd = izodd + inodd;
  if (iodd >= 2) {
    //  Even-even nucleus
    if (ja != jaold) {
      jaold = ja;
      dja = (G4double) ja;
      sqatar = sqrt(dja);
    }
    *eex1st = .012 / sqatar;
    *eex2nd = *eex1st * 2.;
  } else if (iodd > 0) {
    //  even-odd nucleus
    if (ja != jaold) {
      jaold = ja;
      dja = (G4double) ja;
      sqatar = sqrt(dja);
    }
    *eex1st = .012 / sqatar;
    *eex2nd = *eex1st;
  } else {
    //  odd-odd nucleus
    *eex2nd = 0.;
    *eex1st = 0.;
  }
  return;

} // ExcitationEnergyLevel


G4double G4MuonMinusCaptureAtRest::CalculateIsotopicMass(G4double a,
							 G4double z)
{
  // Initialized data

  static G4bool lfirst = true;

  // System generated locals
  G4double ret_val, d__1, d__2;

  // Local variables
  static G4int n;
  static G4double a13, ec, es, ev;
  static G4int iz0;
  static G4double am13, eex, am2zoa;
  static G4double exhydr, exneut;

  static G4double cam2[130] = {
    26.17, 19.25, 24.21, 20.92, 23.15, 18.01, 19.55, 16.94, 19.73, 17.07,
    18.21, 14.99, 16.01, 12.04, 13.27, 11.09, 12.17, 10.26, 11.04,  8.41,
     9.79,  7.36,  8.15,  5.63,  5.88,  3.17,  3.32,  0.82,  1.83,  0.97,
     2.33,  1.27,  2.92,  1.61,  2.91,  1.35,  2.40,  0.89,  1.74,  0.36,
     0.95, -0.65, -0.04, -1.73, -0.96, -2.87, -2.05, -4.05, -3.40, -5.72,
    -3.75, -4.13, -2.42, -2.85, -1.01, -1.33,  0.54, -0.02,  1.74,  0.75,
     2.24,  1.00,  1.98,  0.79,  1.54,  0.39,  1.08,  0.00,  0.78, -0.35,
     0.58, -0.55,  0.59, -0.61,  0.59, -0.35,  0.32, -0.96, -0.52, -2.08,
    -2.46, -3.64, -1.55, -0.96,  0.97,  0.88,  2.37,  1.75,  2.72,  1.90,
     2.55,  1.46,  1.93,  0.86,  1.17,  0.08,  0.39, -0.76, -0.39, -1.51,
    -1.17, -2.36, -1.95, -3.06, -2.62, -3.55, -2.95, -3.75, -3.07, -3.79,
    -3.06, -3.77, -3.05, -3.78, -3.12, -3.90, -3.35, -4.24, -3.86, -4.92,
    -5.06, -6.77, -7.41, -9.18,-10.16,-11.12, -9.76, -9.23, -7.96, -7.65
  };
  static G4double cam3[200] = {
    -8.32,-15.90,-11.51,-14.31,-11.57,-15.90,-13.91,-16.03,-12.13,-13.87,
   -12.25,-14.40,-13.07,-15.80,-13.81,-14.98,-12.63,-13.76,-11.37,-12.38,
    -9.23, -9.65, -7.64, -9.17, -8.05, -9.72, -8.87,-10.76, -8.64, -8.89,
    -6.60, -7.13, -4.77, -5.33, -3.06, -3.79, -1.72, -2.79, -0.93, -2.19,
    -0.52, -1.90, -0.45, -2.20, -1.22, -3.07, -2.42, -4.37, -3.94, -6.08,
    -4.49, -4.50, -3.14, -2.93, -1.04, -1.36,  0.69,  0.21,  2.11,  1.33,
     3.29,  2.46,  4.30,  3.32,  4.79,  3.62,  4.97,  3.64,  4.63,  3.07,
     4.06,  2.49,  3.30,  1.46,  2.06,  0.51,  0.74, -1.18, -1.26, -3.54,
    -3.97, -5.26, -4.18, -3.71, -2.10, -1.70, -0.08, -0.18,  0.94,  0.27,
     1.13,  0.08,  0.91, -0.31,  0.49, -0.78,  0.08, -1.15, -0.23, -1.41,
    -0.42, -1.55, -0.55, -1.66, -0.66, -1.73, -0.75, -1.74, -0.78, -1.69,
    -0.78, -1.60, -0.75, -1.46, -0.67, -1.26, -0.51, -1.04, -0.53, -1.84,
    -2.42, -4.52, -4.76, -6.33, -6.76, -7.81, -5.80, -5.37, -3.63, -3.35,
    -1.75, -1.88, -0.61, -0.90,  0.09, -0.32,  0.55, -0.13,  0.70, -0.06,
     0.49, -0.20,  0.40, -0.22,  0.36, -0.09,  0.58,  0.12,  0.75,  0.15,
     0.70,  0.17,  1.11,  0.89,  1.85,  1.62,  2.54,  2.29,  3.20,  2.91,
     3.84,  3.53,  4.48,  4.15,  5.12,  4.78,  5.75,  5.39,  6.31,  5.91,
     6.87,  6.33,  7.13,  6.61,  7.30,  6.31,  6.27,  4.83,  4.49,  2.85,
     2.32,  0.58, -0.11, -0.98,  0.81,  1.77,  3.37,  4.13,  5.60,  6.15,
     7.29,  7.35,  7.95,  7.67,  8.16,  7.83,  8.31,  8.01,  8.53,  8.27
  };

  if (lfirst) {
    lfirst = false;
    exhydr = GetIsotopicMass(1.0, 1.0);
    exneut = GetIsotopicMass(1.0, 0.0);
  }
  iz0 = NINT(z);
  if (iz0 <= 0) {
    ret_val = a * exneut;
    return ret_val;
  }
  d__1 = a - z;
  n = NINT(d__1);
  if (n <= 0) {
    ret_val = z * exhydr;
    return ret_val;
  }
  am2zoa = (a - z - z) / a;
  am2zoa *= am2zoa;
  a13 = residualNuclearMasses[NINT(a) - 1];
  am13 = 1. / a13;
  ev = (1. - am2zoa * 1.84619) * -17.0354 * a;
  es =
    (1. - am2zoa * 1.712185) * 25.8357 *
    (1. - am13 * .62025 * am13) * (a13 * a13 - .62025);
  ec =
    z * .799 * (z - 1.) * am13 *
    (((am13 * 1.5772 + 1.2273) * am13 - 1.5849) * am13 * am13 + 1.);
  eex = am13 * -.4323 * pow(z, 1.3333333) *
    (((am13 * .49597 - .14518) * am13 - .57811) * am13 + 1.);
  ret_val =
    a * 8.367 - z * .783 + ev + es + ec + eex + cam2[iz0 - 1] + cam3[n - 1];
  ret_val = (ret_val + a * O16OLD) * O16RAT - a * ( C12NEW - ADJUST );
  d__1 = ret_val, d__2 = z * exhydr + (a - z) * exneut;
  ret_val = G4std::min(d__1,d__2);
  return ret_val;

} // CalculateIsotopicMass


void G4MuonMinusCaptureAtRest::EvaporationDeexcitation()
{
  // System generated locals
  G4double d__1, d__2;

  // Local variables
  static G4double etapcm, aiamom, egroun, deltae, hhh, umo, bre2;
  static G4int ipar, jpar;
  static G4double rndm[2];
  static G4double etax, etay, etaz, pcms, umev, p2res, cfe2e1, cfm1e1,
    echck, eegcm, gamcm, bre1m1, delta, ahelp, freje;
  static G4int jamin;
  static G4double energ;
  static G4int iamin, jamom;
  static G4double enmin;
  static G4double erncm, hhhsq, phelp, plbgx, plbgy, rnucl;
  static G4int lexpn;
  static G4double dismx, pcmsx, pcmsy, echck0;
  static G4int lmult;
  static G4double xtent, pcmsz, plbgz, energ0, eex2nd, roten0, eex1st,
    delpai, cosgam[3];
  static G4int ibhelp;
  static G4double eexdum;
  static G4int ichelp;
  static G4double ainerm, rnmass, asmall, aogmax, aogmin, tempsq, temper,
    ddlexp, xdismx;
  static G4int naiam;

  // ---------------------------------------------------------------------*
  //                                                                      *
  //  EvaporationDeexcitation :  created on 5-10-1990 by A. Ferrari       *
  //                                               & P. Sala, INFN Milan  *
  //                                                                      *
  //    Last change  on  28-jan-93   by   Alfredo Ferrari, INFN-Milan     *
  //                                                                      *
  //    This routine provides a simple model for sampling nuclear deexci- *
  //    tation gammas following the evaporation Step                      *
  //                                                                      *
  // ---------------------------------------------------------------------*

  // ---------------------------------------------------------------------*
  //  Entering the routine we have:                                       *
  //            Ammres is the atomic mass of the residual nucleus         *
  //            Ibres (Anow) its mass number                              *
  //            Icres (Znow) its atomic number                            *
  //            Eres         its total energy                             *
  //            Ptres        its momentum                                 *
  //            Pxres        x-component of the momentum                  *
  //            Pyres        y-component of the momentum                  *
  //            Pzres        z-component of the momentum                  *
  //            Tvrecl       kinetic energy                               *
  //            Tvcms        excitation energy                            *
  // ---------------------------------------------------------------------*

  if (atMassResidNucl <= 0 || excitEnResidNucl <= GAMMIN) {
    excitEnResidNucl = 0.;
    return;
  }
  echck = totEnResidNucl - massResidNucl;
  echck0 = echck;
  ibhelp = atMassResidNucl / 2;
  if (ibhelp << 1 < atMassResidNucl) {
    ExcitationEnergyLevel(atMassResidNucl, chargeResidNucl, &delta,
			  &eex2nd, &eexdum);
    ipar = 1;
    jpar = 1;
  } else {
    ichelp = chargeResidNucl / 2;
    jpar = 0;
    if (ichelp << 1 < chargeResidNucl) {
      ExcitationEnergyLevel(atMassResidNucl, chargeResidNucl, &delta,
			    &eex1st, &eexdum);
      ipar = 1;
    } else {
      ExcitationEnergyLevel(atMassResidNucl, chargeResidNucl, &eex1st,
			    &delta, &eexdum);
      ipar = 2;
    }
  }
  delpai = G4std::min(eexdum,delta);
  rnucl = residualNuclearMasses[atMassResidNucl - 1] * R0NUCL;
  ainerm = massResidNucl * .24 * rnucl * rnucl;
  roten0 = 0.5 * PLABRC * PLABRC / ainerm;
  d__1 = delta, d__2 = roten0 * 2.;
  enmin = G4std::max(d__1,d__2);
  rnmass = massResidNucl + excitEnResidNucl;
  umo = rnmass;
  gamcm = totEnResidNucl / rnmass;
  etax = momXResidNucl / rnmass;
  etay = momYResidNucl / rnmass;
  etaz = momZResidNucl / rnmass;
  if (excitEnResidNucl > enmin && atMassResidNucl > 4) {
    ahelp =
      residualNuclearMasses[atMassResidNucl - 1] * HNDFE1 *
      residualNuclearMasses[atMassResidNucl - 1];
    cfm1e1 = C0M1E1 * HNDFM1 / ahelp;
    cfe2e1 = C0E2E1 * HNDFE2 / ahelp;
    do {
      umev = excitEnResidNucl * 1e3;
      asmall =
	LevelDensity(chargeResidNucl, atMassResidNucl - chargeResidNucl,
		     &aogmax, &aogmin);
      tempsq = (excitEnResidNucl - delpai) * .001 / asmall;
      temper = sqrt(tempsq);
      hhh = excitEnResidNucl / temper;
      hhhsq = hhh * hhh;
      if (hhh > EXPMAX) {
	ahelp = 0.;
      } else {
	ahelp = exp(-hhh);
      }
      bre1m1 = 6. - ahelp * (hhhsq * hhh + hhhsq * 3. + hhh * 6.);
      bre2 = bre1m1 * 20. - ahelp * (hhhsq * hhhsq * hhh + hhhsq * 5. * hhhsq);
      bre1m1 *= cfm1e1 + 1.;
      bre2 = bre2 * tempsq * cfe2e1;
      rndm[0] = G4UniformRand();
      if (rndm[0] < bre1m1 / (bre1m1 + bre2)) {
	lmult = 1;
      } else {
	lmult = 2;
      }
      lexpn = (lmult << 1) + 1;
      ddlexp = (G4double) lexpn;
      xdismx = G4std::min(ddlexp,hhh);
      if (xdismx > EXPMAX) {
	dismx = 0.;
      } else {
	dismx = pow(xdismx, (G4double)lexpn) * exp(-xdismx);
      }
      do {
	rndm[0] = G4UniformRand();
	rndm[1] = G4UniformRand();
	xtent = rndm[0] * hhh;
	if (xtent > EXPMAX) {
	  freje = 0.;
	} else {
	  freje = pow(xtent, (G4double)lexpn) * exp(-xtent) / dismx;
	}
      } while (rndm[1] >= freje);
      energ0 = xtent * temper;
      excitEnResidNucl -= energ0;
      rnmass = massResidNucl + excitEnResidNucl;
      RanDirCos(cosgam, &cosgam[1], &cosgam[2]);
      erncm = (umo * umo + rnmass * rnmass) * .5 / umo;
      eegcm = umo - erncm;
      pcms = eegcm;
      pcmsx = pcms * cosgam[0];
      pcmsy = pcms * cosgam[1];
      pcmsz = pcms * cosgam[2];
      etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
      energ = gamcm * eegcm + etapcm;
      phelp = etapcm / (gamcm + 1.) + eegcm;
      plbgx = pcmsx + etax * phelp;
      plbgy = pcmsy + etay * phelp;
      plbgz = pcmsz + etaz * phelp;
      totEnResidNucl = gamcm * erncm - etapcm;
      kinEnResidNucl = totEnResidNucl - rnmass;
      phelp = -etapcm / (gamcm + 1.) + erncm;
      momXResidNucl = -pcmsx + etax * phelp;
      momYResidNucl = -pcmsy + etay * phelp;
      momZResidNucl = -pcmsz + etaz * phelp;
      echck -= energ;
      ++nSecPart;
      Secondaries[nSecPart - 1].SetZero();
      Secondaries[nSecPart - 1].SetMass( massGamma );
      Secondaries[nSecPart - 1].SetMomentumAndUpdate( plbgx, plbgy, plbgz );
      Secondaries[nSecPart - 1].SetParticleDef( pdefGamma );
      gamcm = totEnResidNucl / rnmass;
      etax = momXResidNucl / rnmass;
      etay = momYResidNucl / rnmass;
      etaz = momZResidNucl / rnmass;
      umo = rnmass;
    } while (excitEnResidNucl > enmin);
    if (nSecPart >= MXSECS) {
      G4cout << " **** Finuc overflow in EvaporationDeexcitation, program stopped ****" << G4endl;
      G4cout << "STOP" << G4endl;
      exit(1);
    }
  }
  if (excitEnResidNucl <= ((jpar << 1) + 6) * roten0) {
    energ0 = excitEnResidNucl;
    excitEnResidNucl = 0.;
    RanDirCos(cosgam, &cosgam[1], &cosgam[2]);
    erncm = (umo * umo + massResidNucl * massResidNucl) * .5 / umo;
    eegcm = umo - erncm;
    pcms = eegcm;
    pcmsx = pcms * cosgam[0];
    pcmsy = pcms * cosgam[1];
    pcmsz = pcms * cosgam[2];
    etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
    energ = gamcm * eegcm + etapcm;
    phelp = etapcm / (gamcm + 1.) + eegcm;
    plbgx = pcmsx + etax * phelp;
    plbgy = pcmsy + etay * phelp;
    plbgz = pcmsz + etaz * phelp;
    totEnResidNucl = gamcm * erncm - etapcm;
    kinEnResidNucl = totEnResidNucl - massResidNucl;
    phelp = -etapcm / (gamcm + 1.) + erncm;
    momXResidNucl = -pcmsx + etax * phelp;
    momYResidNucl = -pcmsy + etay * phelp;
    momZResidNucl = -pcmsz + etaz * phelp;
    echck -= energ;
    ++nSecPart;
    Secondaries[nSecPart - 1].SetZero();
    Secondaries[nSecPart - 1].SetMass( massGamma );
    Secondaries[nSecPart - 1].SetMomentumAndUpdate( plbgx, plbgy, plbgz );
    Secondaries[nSecPart - 1].SetParticleDef( pdefGamma );
  } else {
    aiamom = sqrt(excitEnResidNucl / roten0 + .25 + jpar * .75) - .25;
    naiam = G4int(aiamom);
    if (ipar == 1) {
      iamin = naiam;
      if (jpar == 1) {
	jamin = (naiam << 1) + 1;
	if (iamin % 2 > 0) {
	  egroun = roten0 * 3.75;
	} else {
	  egroun = roten0 * .75;
	}
      } else {
	jamin = iamin << 1;
	if (iamin % 2 > 0) {
	  egroun = roten0 * 2.;
	} else {
	  egroun = 0.;
	}
      }
    } else {
      iamin = ( naiam / ipar ) * ipar;
      jamin = iamin << 1;
      egroun = 0.;
    }
    deltae = excitEnResidNucl + egroun - roten0 * .25 * jamin * (jamin + 2);
    for (jamom = jamin; jamom >= 4; jamom += -4) {
      energ0 = roten0 * ((jamom << 1) - 2) + deltae;
      deltae = 0.;
      excitEnResidNucl -= energ0;
      rnmass = massResidNucl + excitEnResidNucl;
      RanDirCos(cosgam, &cosgam[1], &cosgam[2]);
      erncm = (umo * umo + rnmass * rnmass) * .5 / umo;
      eegcm = umo - erncm;
      pcms = eegcm;
      pcmsx = pcms * cosgam[0];
      pcmsy = pcms * cosgam[1];
      pcmsz = pcms * cosgam[2];
      etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
      energ = gamcm * eegcm + etapcm;
      phelp = etapcm / (gamcm + 1.) + eegcm;
      plbgx = pcmsx + etax * phelp;
      plbgy = pcmsy + etay * phelp;
      plbgz = pcmsz + etaz * phelp;
      totEnResidNucl = gamcm * erncm - etapcm;
      kinEnResidNucl = totEnResidNucl - rnmass;
      phelp = -etapcm / (gamcm + 1.) + erncm;
      momXResidNucl = -pcmsx + etax * phelp;
      momYResidNucl = -pcmsy + etay * phelp;
      momZResidNucl = -pcmsz + etaz * phelp;
      echck -= energ;
      ++nSecPart;
      Secondaries[nSecPart - 1].SetZero();
      Secondaries[nSecPart - 1].SetMass( massGamma );
      Secondaries[nSecPart - 1].SetMomentumAndUpdate( plbgx, plbgy, plbgz );
      Secondaries[nSecPart - 1].SetParticleDef( pdefGamma );
      gamcm = totEnResidNucl / rnmass;
      etax = momXResidNucl / rnmass;
      etay = momYResidNucl / rnmass;
      etaz = momZResidNucl / rnmass;
      umo = rnmass;
    }
  }
  echck -= kinEnResidNucl;
  if (abs(echck) > echck0 * 1e-7) {
    G4cout << " **** No energy conservation in EvaporationDeexcitation ****"
	 << " " << echck << G4endl;
  }
  excitEnResidNucl = 0.;
  if (nSecPart > MXSECS) {
    G4cout << " **** Finuc overflow in EvaporationDeexcitation, program stopped ****" << G4endl;
    G4cout << "STOP" << G4endl;
    exit(1);
  }
  recoilEnResidNucl = kinEnResidNucl;
  p2res =
    momXResidNucl * momXResidNucl + momYResidNucl * momYResidNucl +
    momZResidNucl * momZResidNucl;
  totMomResidNucl = sqrt(p2res);
  return;

} // EvaporationDeexcitation


void G4MuonMinusCaptureAtRest::FermiMotion(G4int pidx)
{

  // System generated locals
  G4double r__1;

  // Local variables
  static G4int iztemp;
  static G4double delctr;
  static G4double p2;
  static G4double cfe, sfe, p2sq, ferm, polc;
  static G4double pols;
  static G4double frndm[3];
  static G4double atemp, ztemp;
  static G4double tveuz, v0extr;

  ferm = nucleonMaxFermiEn[pidx - 1];
  frndm[0] = G4UniformRand();
  frndm[1] = G4UniformRand();
  frndm[2] = G4UniformRand();
  r__1 = G4std::max(frndm[0],frndm[1]);
  p2 = G4std::max(r__1,frndm[2]);
  if (atMassTarget <= 1) {
    ferm = 0.;
  }
  p2 = ferm * p2;
  p2sq = p2 * p2;
  atemp = (G4double) atMassTarget - 1.;
  if (pidx == 1) {
    iztemp = chargeTarget - 1;
  } else {
    iztemp = chargeTarget;
  }
  ztemp = (G4double) iztemp;
  delctr = ((G4double) chargeTarget - ztemp) * AMELEC;
  massResidNucl = atemp * AMUAMU + GetIsotopicMass(atemp, ztemp) * .001;
  RanPolarAng(&polc, &pols);
  RanAzimuthalAng(&sfe, &cfe);
  nucleonFermiMomX = cfe * pols * p2;
  nucleonFermiMomY = sfe * pols * p2;
  nucleonFermiMomZ = polc * p2;
  nucleonFermiEn = sqrt(nucleonMassSquared[pidx - 1] + p2sq);
  tveuz =
    nucleonPotWell[pidx - 1] - nucleonFermiEn +
    nucleonBindingEn[pidx - 1] + massTarget -
    massResidNucl - delctr;
  if (tveuz < 0.) {
    v0extr = -tveuz + 10. * TVEPSI;
    nucleonPotWell[pidx - 1] += v0extr;
  }
  return;

} // FermiMotion


void G4MuonMinusCaptureAtRest::ResidualNucleusCascade(G4int m2, G4int m3,
						      G4double t1, G4double *u,
						      G4double *erec,
						      G4bool *loppar)
{
  // Initialized data

  static G4int ievevp = 0;
  static G4bool lfirst = true;

  // System generated locals
  G4int i__1;
  G4double d__1, d__2, d__3, d__4;

  // Local variables
  static G4double a, c[3];
  static G4int i, j, k;
  static G4double q[7], r[6], s[6], v, z, e1, e2;
  static G4int n1;
  static G4double aa;
  static G4int ja;
  static G4double ar;
  static G4int jj, kk, mm, nn, js;
  static G4double um;
  static G4int jz;
  static G4double zr, zz;
  static G4int iaa, jja;
  static G4double emh, arg;
  static G4int jat;
  static G4double emn;
  static G4int jjn;
  static G4double fjs, sas, eps, ses;
  static G4int jjz;
  static G4double umo, sos[6], sum;
  static G4int jzt;
  static G4double sus;
  static G4int izz;
  static G4double eye1[6], eye0[6], umo2, ajja, fact;
  static G4double emhn, ecms;
  static G4double rndm[2];
  static G4double etax, pcms, etay, etaz;
  static G4double corr, uran, unew, umax, umin;
  static G4int ipro, ineu;
  static G4double zjjz, ccou2, umin2, p2res, smom1[6], ddjja, gamcm,
    etacm, eepcm, sigma, ccoul[6];
  static G4double deltu, emnum, erncm, ddjjz, pcmsx, qnorm, zmass[6],
    pcmsy, pcmsz, phelp, plbpx, plbpy, plbpz;
  static G4int itemp;
  static G4double strun[6];
  static G4double epsav;
  static G4double energ0;
  static G4int imass;
  static G4double eex2nd, ueu3he, rnmas0, ratio2, eex1st, z2mass[6],
    flkcou[6], thresh[6], smalla[6], bnmass[6], corrrr[6];
  static G4double rnmass, elbtot;
  static G4int jresid, jemiss;
  static G4double deudeu, protri, deupro, prprne, deuneu, prnene, etapcm;
  static G4double umxres, asmmax, asmmin;
  static G4double tmpvar, expsas, expsus;
  static G4int ncount;
  static G4double coslbp[3];
  static G4double exmass[6];

  static G4int ia[6] = { 1, 1, 2, 3, 3, 4 };
  static G4int iz[6] = { 0, 1, 1, 1, 2, 2 };
  static G4double fla[6] = { 1., 1., 2., 3., 3., 4. };
  static G4double flz[6] = { 0., 1., 1., 1., 2., 2. };

  static G4double p2[1001] = {
    0.00000, 0.00000, 0.00326851, 0.00736486, 0.0129026,
    0.0200330, 0.0286195, 0.0386427, 0.0500610, 0.0628421,
    0.0769436, 0.0923297, 0.108963, 0.126807, 0.145824,
    0.165979, 0.187236, 0.209560, 0.232915, 0.257269,
    0.282588, 0.308839, 0.335990, 0.364010, 0.392867,
    0.422532, 0.452976, 0.484169, 0.516085, 0.548696,
    0.581975, 0.615897, 0.650438, 0.685572, 0.721278,
    0.757531, 0.794311, 0.831596, 0.869366, 0.907601,
    0.946281, 0.985390, 1.02491, 1.06482, 1.10511,
    1.14576, 1.18675, 1.22808, 1.26972, 1.31166,
    1.35390, 1.39642, 1.43920, 1.48224, 1.52552,
    1.56904, 1.61278, 1.65674, 1.70090, 1.74526,
    1.78980, 1.83453, 1.87943, 1.92450, 1.96973,
    2.01510, 2.06063, 2.10629, 2.15209, 2.19802,
    2.24407, 2.29024, 2.33652, 2.38291, 2.42941,
    2.47601, 2.52270, 2.56949, 2.61636, 2.66332,
    2.71037, 2.75749, 2.80469, 2.85197, 2.89931,
    2.94672, 2.99420, 3.04174, 3.08934, 3.13700,
    3.18472, 3.23250, 3.28032, 3.32820, 3.37612,
    3.42410, 3.47211, 3.52018, 3.56828, 3.61643,
    3.66462, 3.71263, 3.76091, 3.80923, 3.85758,
    3.90596, 3.95438, 4.00283, 4.05130, 4.09981,
    4.14835, 4.19691, 4.24551, 4.29412, 4.34277,
    4.39144, 4.44013, 4.48885, 4.53759, 4.58635,
    4.63513, 4.68394, 4.73276, 4.78161, 4.83047,
    4.87936, 4.92826, 4.97719, 5.02612, 5.07508,
    5.12406, 5.17305, 5.22205, 5.27108, 5.32011,
    5.36917, 5.41823, 5.46732, 5.51641, 5.56552,
    5.61465, 5.66378, 5.71293, 5.76209, 5.81127,
    5.86045, 5.90965, 5.95886, 6.00808, 6.05731,
    6.10655, 6.15581, 6.20507, 6.25434, 6.30363,
    6.35292, 6.40222, 6.45153, 6.50085, 6.55018,
    6.59952, 6.64887, 6.69822, 6.74759, 6.79696,
    6.84634, 6.89573, 6.94512, 6.99453, 7.04394,
    7.09336, 7.14278, 7.19221, 7.24165, 7.29110,
    7.34055, 7.39001, 7.43947, 7.48894, 7.53842,
    7.58791, 7.63740, 7.68689, 7.73639, 7.78590,
    7.83541, 7.88493, 7.93445, 7.98398, 8.03352,
    8.08306, 8.13260, 8.18215, 8.23170, 8.28126,
    8.33082, 8.38039, 8.42997, 8.47954, 8.52912,
    8.57871, 8.62830, 8.67789, 8.72749, 8.77710,
    8.82670, 8.87631, 8.92593, 8.97555, 9.02517,
    9.07480, 9.12442, 9.17406, 9.22370, 9.27334,
    9.32298, 9.37263, 9.42228, 9.47193, 9.52159,
    9.57125, 9.62091, 9.67058, 9.72025, 9.76993,
    9.81960, 9.86928, 9.91896, 9.96865, 10.0183,
    10.0680, 10.1177, 10.1674, 10.2171, 10.2668,
    10.3165, 10.3662, 10.4159, 10.4656, 10.5154,
    10.5651, 10.6148, 10.6645, 10.7142, 10.7640,
    10.8137, 10.8634, 10.9132, 10.9629, 11.0126,
    11.0624, 11.1121, 11.1619, 11.2116, 11.2614,
    11.3111, 11.3609, 11.4106, 11.4604, 11.5101,
    11.5599, 11.6097, 11.6594, 11.7092, 11.7589,
    11.8087, 11.8585, 11.9083, 11.9580, 12.0078,
    12.0576, 12.1074, 12.1572, 12.2069, 12.2567,
    12.3065, 12.3563, 12.4061, 12.4559, 12.5057,
    12.5555, 12.6053, 12.6551, 12.7049, 12.7547,
    12.8045, 12.8543, 12.9041, 12.9539, 13.0037,
    13.0535, 13.1033, 13.1531, 13.2029, 13.2527,
    13.3026, 13.3524, 13.4022, 13.4520, 13.5018,
    13.5516, 13.6015, 13.6513, 13.7011, 13.7509,
    13.8008, 13.8506, 13.9004, 13.9503, 14.0001,
    14.0499, 14.0998, 14.1496, 14.1994, 14.2493,
    14.2991, 14.3490, 14.3988, 14.4486, 14.4985,
    14.5483, 14.5982, 14.6480, 14.6979, 14.7477,
    14.7976, 14.8474, 14.8973, 14.9471, 14.9970,
    15.0468, 15.0967, 15.1465, 15.1964, 15.2462,
    15.2961, 15.3460, 15.3958, 15.4457, 15.4955,
    15.5454, 15.5953, 15.6451, 15.6950, 15.7449,
    15.7947, 15.8446, 15.8945, 15.9443, 15.9942,
    16.0441, 16.0939, 16.1438, 16.1937, 16.2435,
    16.2934, 16.3433, 16.3932, 16.4430, 16.4929,
    16.5428, 16.5927, 16.6426, 16.6924, 16.7423,
    16.7922, 16.8421, 16.8920, 16.9418, 16.9917,
    17.0416, 17.0915, 17.1414, 17.1913, 17.2411,
    17.2910, 17.3409, 17.3908, 17.4407, 17.4906,
    17.5405, 17.5904, 17.6403, 17.6902, 17.7401,
    17.7899, 17.8398, 17.8897, 17.9396, 17.9895,
    18.0394, 18.0893, 18.1392, 18.1891, 18.2390,
    18.2889, 18.3388, 18.3887, 18.4386, 18.4885,
    18.5384, 18.5883, 18.6382, 18.6881, 18.7380,
    18.7879, 18.8378, 18.8877, 18.9376, 18.9875,
    19.0374, 19.0874, 19.1373, 19.1872, 19.2371,
    19.2870, 19.3369, 19.3868, 19.4367, 19.4866,
    19.5365, 19.5865, 19.6364, 19.6863, 19.7362,
    19.7861, 19.8360, 19.8859, 19.9358, 19.9857,
    20.0357, 20.0856, 20.1355, 20.1854, 20.2353,
    20.2852, 20.3352, 20.3851, 20.4350, 20.4849,
    20.5348, 20.5848, 20.6347, 20.6846, 20.7345,
    20.7844, 20.8344, 20.8843, 20.9342, 20.9841,
    21.0340, 21.0840, 21.1339, 21.1838, 21.2337,
    21.2837, 21.3336, 21.3835, 21.4334, 21.4834,
    21.5333, 21.5832, 21.6331, 21.6831, 21.7330,
    21.7829, 21.8329, 21.8828, 21.9327, 21.9826,
    22.0326, 22.0825, 22.1324, 22.1824, 22.2323,
    22.2822, 22.3321, 22.3821, 22.4320, 22.4819,
    22.5319, 22.5818, 22.6317, 22.6817, 22.7316,
    22.7815, 22.8315, 22.8814, 22.9313, 22.9813,
    23.0312, 23.0811, 23.1311, 23.1810, 23.2309,
    23.2809, 23.3308, 23.3808, 23.4307, 23.4806,
    23.5306, 23.5805, 23.6304, 23.6804, 23.7303,
    23.7803, 23.8302, 23.8801, 23.9301, 23.9800,
    24.0300, 24.0799, 24.1298, 24.1798, 24.2297,
    24.2797, 24.3296, 24.3795, 24.4295, 24.4794,
    24.5294, 24.5793, 24.6293, 24.6792, 24.7291,
    24.7791, 24.8290, 24.8790, 24.9289, 24.9789,
    25.0288, 25.0788, 25.1287, 25.1786, 25.2286,
    25.2785, 25.3285, 25.3784, 25.4284, 25.4783,
    25.5283, 25.5782, 25.6282, 25.6781, 25.7280,
    25.7780, 25.8279, 25.8779, 25.9278, 25.9778,
    26.0277, 26.0777, 26.1276, 26.1776, 26.2275,
    26.2775, 26.3274, 26.3774, 26.4273, 26.4773,
    26.5272, 26.5772, 26.6271, 26.6771, 26.7270,
    26.7770, 26.8269, 26.8769, 26.9268, 26.9768,
    27.0267, 27.0767, 27.1266, 27.1766, 27.2265,
    27.2765, 27.3265, 27.3764, 27.4264, 27.4763,
    27.5263, 27.5762, 27.6262, 27.6761, 27.7261,
    27.7760, 27.8260, 27.8760, 27.9259, 27.9759,
    28.0258, 28.0758, 28.1257, 28.1757, 28.2256,
    28.2756, 28.3256, 28.3755, 28.4255, 28.4754,
    28.5254, 28.5753, 28.6253, 28.6752, 28.7252,
    28.7752, 28.8251, 28.8751, 28.9250, 28.9750,
    29.0250, 29.0749, 29.1249, 29.1748, 29.2248,
    29.2747, 29.3247, 29.3747, 29.4246, 29.4746,
    29.5246, 29.5745, 29.6245, 29.6744, 29.7244,
    29.7744, 29.8243, 29.8743, 29.9242, 29.9742,
    30.0242, 30.0741, 30.1241, 30.1740, 30.2240,
    30.2740, 30.3239, 30.3739, 30.4238, 30.4738,
    30.5238, 30.5737, 30.6237, 30.6737, 30.7236,
    30.7736, 30.8235, 30.8735, 30.9235, 30.9734,
    31.0234, 31.0733, 31.1233, 31.1733, 31.2233,
    31.2732, 31.3232, 31.3731, 31.4231, 31.4731,
    31.5230, 31.5730, 31.6230, 31.6729, 31.7229,
    31.7729, 31.8228, 31.8728, 31.9227, 31.9727,
    32.0227, 32.0726, 32.1226, 32.1726, 32.2225,
    32.2725, 32.3225, 32.3724, 32.4224, 32.4724,
    32.5223, 32.5723, 32.6223, 32.6722, 32.7222,
    32.7722, 32.8221, 32.8721, 32.9221, 32.9720,
    33.0220, 33.0720, 33.1219, 33.1719, 33.2219,
    33.2718, 33.3218, 33.3718, 33.4218, 33.4717,
    33.5217, 33.5717, 33.6216, 33.6716, 33.7216,
    33.7715, 33.8215, 33.8715, 33.9214, 33.9714,
    34.0214, 34.0713, 34.1213, 34.1713, 34.2213,
    34.2712, 34.3212, 34.3712, 34.4211, 34.4711,
    34.5211, 34.5710, 34.6210, 34.6710, 34.7209,
    34.7709, 34.8209, 34.8708, 34.9208, 34.9708,
    35.0208, 35.0707, 35.1207, 35.1707, 35.2207,
    35.2706, 35.3206, 35.3706, 35.4205, 35.4705,
    35.5205, 35.5704, 35.6204, 35.6704, 35.7204,
    35.7703, 35.8203, 35.8703, 35.9203, 35.9702,
    36.0202, 36.0702, 36.1201, 36.1701, 36.2201,
    36.2701, 36.3200, 36.3700, 36.4200, 36.4700,
    36.5199, 36.5699, 36.6199, 36.6698, 36.7198,
    36.7698, 36.8198, 36.8697, 36.9197, 36.9697,
    37.0197, 37.0696, 37.1196, 37.1696, 37.2196,
    37.2695, 37.3195, 37.3695, 37.4194, 37.4694,
    37.5194, 37.5694, 37.6194, 37.6693, 37.7193,
    37.7693, 37.8192, 37.8692, 37.9192, 37.9692,
    38.0191, 38.0691, 38.1191, 38.1691, 38.2191,
    38.2690, 38.3190, 38.3690, 38.4189, 38.4689,
    38.5189, 38.5689, 38.6189, 38.6688, 38.7188,
    38.7688, 38.8188, 38.8687, 38.9187, 38.9687,
    39.0187, 39.0686, 39.1186, 39.1686, 39.2186,
    39.2685, 39.3185, 39.3685, 39.4185, 39.4684,
    39.5184, 39.5684, 39.6184, 39.6684, 39.7183,
    39.7683, 39.8183, 39.8683, 39.9182, 39.9682,
    40.0182, 40.0682, 40.1181, 40.1681, 40.2181,
    40.2681, 40.3181, 40.3680, 40.4180, 40.4680,
    40.5180, 40.5679, 40.6179, 40.6679, 40.7179,
    40.7679, 40.8178, 40.8678, 40.9178, 40.9678,
    41.0178, 41.0677, 41.1177, 41.1677, 41.2177,
    41.2677, 41.3176, 41.3676, 41.4176, 41.4676,
    41.5175, 41.5675, 41.6175, 41.6675, 41.7175,
    41.7674, 41.8174, 41.8674, 41.9174, 41.9674,
    42.0173, 42.0673, 42.1173, 42.1673, 42.2173,
    42.2672, 42.3172, 42.3672, 42.4172, 42.4671,
    42.5171, 42.5671, 42.6171, 42.6671, 42.7171,
    42.7670, 42.8170, 42.8670, 42.9170, 42.9669,
    43.0169, 43.0669, 43.1169, 43.1669, 43.2169,
    43.2668, 43.3168, 43.3668, 43.4168, 43.4668,
    43.5167, 43.5667, 43.6167, 43.6667, 43.7167,
    43.7666, 43.8166, 43.8666, 43.9166, 43.9666,
    44.0166, 44.0665, 44.1165, 44.1665, 44.2165,
    44.2665, 44.3164, 44.3664, 44.4164, 44.4664,
    44.5164, 44.5663, 44.6163, 44.6663, 44.7163,
    44.7663, 44.8163, 44.8662, 44.9162, 44.9662,
    45.0162, 45.0662, 45.1161, 45.1661, 45.2161,
    45.2661, 45.3161, 45.3661, 45.4160, 45.4660,
    45.5160, 45.5660, 45.6160, 45.6660, 45.7159,
    45.7659, 45.8159, 45.8659, 45.9159, 45.9659,
    46.0158, 46.0658, 46.1158, 46.1658, 46.2158,
    46.2657, 46.3157, 46.3657, 46.4157, 46.4657,
    46.5157, 46.5656, 46.6156, 46.6656, 46.7156,
    46.7656, 46.8156, 46.8655, 46.9155, 46.9655,
    47.0155, 47.0655, 47.1155, 47.1654, 47.2154,
    47.2654, 47.3154, 47.3654, 47.4154, 47.4653,
    47.5153, 47.5653, 47.6153, 47.6653, 47.7153,
    47.7652, 47.8152, 47.8652, 47.9152, 47.9652,
    48.0152, 48.0652, 48.1151, 48.1651, 48.2151,
    48.2651, 48.3151, 48.3651, 48.4150, 48.4650,
    48.5150
  };

  static G4double rho[6] = { 0.0, 0.0, 0.70588, 0.70588, 0.70588, 0.70588 };

  static G4double alph[297] = {
    2.69000, 2.29184, 2.09819, 1.97582, 1.88867,
    1.82212, 1.76892, 1.72500, 1.68785, 1.65583,
    1.62781, 1.60300, 1.58081, 1.56078, 1.54258,
    1.52592, 1.51060, 1.49643, 1.48328, 1.47102,
    1.45955, 1.44878, 1.43865, 1.42909, 1.42005,
    1.41148, 1.40333, 1.39558, 1.38819, 1.38113,
    1.37438, 1.36791, 1.36171, 1.35575, 1.35002,
    1.34451, 1.33919, 1.33407, 1.32912, 1.32434,
    1.31971, 1.31523, 1.31089, 1.30669, 1.30261,
    1.29865, 1.29480, 1.29106, 1.28742, 1.28388,
    1.28044, 1.27708, 1.27380, 1.27061, 1.26750,
    1.26446, 1.26149, 1.25859, 1.25576, 1.25299,
    1.25028, 1.24763, 1.24504, 1.24250, 1.24001,
    1.23758, 1.23519, 1.23285, 1.23055, 1.22830,
    1.22609, 1.22392, 1.22180, 1.21971, 1.21765,
    1.21564, 1.21366, 1.21171, 1.20979, 1.20791,
    1.20606, 1.20424, 1.20245, 1.20069, 1.19895,
    1.19724, 1.19556, 1.19391, 1.19228, 1.19067,
    1.18909, 1.18752, 1.18599, 1.18447, 1.18298,
    1.18150, 1.18005, 1.17861, 1.17720, 1.17581,
    1.17443, 1.17307, 1.17173, 1.17040, 1.16910,
    1.16781, 1.16653, 1.16527, 1.16403, 1.16280,
    1.16159, 1.16039, 1.15921, 1.15804, 1.15688,
    1.15574, 1.15460, 1.15349, 1.15238, 1.15129,
    1.15021, 1.14914, 1.14808, 1.14703, 1.14600,
    1.14498, 1.14396, 1.14296, 1.14197, 1.14099,
    1.14001, 1.13905, 1.13810, 1.13716, 1.13622,
    1.13530, 1.13438, 1.13348, 1.13258, 1.13169,
    1.13081, 1.12994, 1.12907, 1.12822, 1.12737,
    1.12653, 1.12569, 1.12487, 1.12405, 1.12324,
    1.12244, 1.12164, 1.12085, 1.12007, 1.11929,
    1.11852, 1.11776, 1.11700, 1.11625, 1.11551,
    1.11477, 1.11404, 1.11331, 1.11259, 1.11188,
    1.11117, 1.11047, 1.10977, 1.10908, 1.10840,
    1.10772, 1.10704, 1.10637, 1.10571, 1.10505,
    1.10439, 1.10374, 1.10310, 1.10246, 1.10182,
    1.10119, 1.10056, 1.09994, 1.09933, 1.09871,
    1.09811, 1.09750, 1.09690, 1.09631, 1.09572,
    1.09513, 1.09455, 1.09397, 1.09339, 1.09282,
    1.09225, 1.09169, 1.09113, 1.09058, 1.09002,
    1.08948, 1.08893, 1.08839, 1.08785, 1.08732,
    1.08679, 1.08626, 1.08574, 1.08522, 1.08470,
    1.08419, 1.08368, 1.08317, 1.08267, 1.08216,
    1.08167, 1.08117, 1.08068, 1.08019, 1.07970,
    1.07922, 1.07874, 1.07826, 1.07779, 1.07732,
    1.07685, 1.07638, 1.07592, 1.07546, 1.07500,
    1.07455, 1.07409, 1.07364, 1.07320, 1.07275,
    1.07231, 1.07187, 1.07143, 1.07100, 1.07057,
    1.07014, 1.06971, 1.06928, 1.06886, 1.06844,
    1.06802, 1.06760, 1.06719, 1.06678, 1.06637,
    1.06596, 1.06556, 1.06515, 1.06475, 1.06435,
    1.06396, 1.06356, 1.06317, 1.06278, 1.06239,
    1.06200, 1.06162, 1.06123, 1.06085, 1.06047,
    1.06010, 1.05972, 1.05935, 1.05898, 1.05861,
    1.05824, 1.05787, 1.05751, 1.05715, 1.05679,
    1.05643, 1.05607, 1.05572, 1.05536, 1.05501,
    1.05466, 1.05431, 1.05396, 1.05362, 1.05328,
    1.05293, 1.05259, 1.05225, 1.05192, 1.05158,
    1.05125, 1.05091, 1.05058, 1.05025, 1.04992,
    3.98462e+28, 2.17105e-56
  };

  static G4double bet[297] = {
    0.598513, 0.434469, 0.356520, 0.308110, 0.274114,
    0.248467, 0.228183, 0.211594, 0.197684, 0.185790,
    0.175462, 0.166378, 0.158303, 0.151061, 0.144516,
    0.138561, 0.133111, 0.128098, 0.123466, 0.119167,
    0.115163, 0.111422, 0.107915, 0.104620, 0.101514,
    0.0985804, 0.0958036, 0.0931700, 0.0906674, 0.0882853,
    0.0860140, 0.0838452, 0.0817713, 0.0797855, 0.0778817,
    0.0760543, 0.0742981, 0.0726089, 0.0709823, 0.0694144,
    0.0679018, 0.0664414, 0.0650300, 0.0636651, 0.0623440,
    0.0610645, 0.0598244, 0.0586217, 0.0574546, 0.0563212,
    0.0552201, 0.0541496, 0.0531083, 0.0520951, 0.0511085,
    0.0501474, 0.0492108, 0.0482976, 0.0474069, 0.0465377,
    0.0456891, 0.0448604, 0.0440508, 0.0432596, 0.0424860,
    0.0417294, 0.0409892, 0.0402648, 0.0395557, 0.0388613,
    0.0381811, 0.0375146, 0.0368614, 0.0362210, 0.0355930,
    0.0349770, 0.0343726, 0.0337795, 0.0331973, 0.0326257,
    0.0320644, 0.0315130, 0.0309713, 0.0304389, 0.0299157,
    0.0294013, 0.0288955, 0.0283981, 0.0279088, 0.0274275,
    0.0269538, 0.0264876, 0.0260288, 0.0255770, 0.0251321,
    0.0246940, 0.0242625, 0.0238374, 0.0234186, 0.0230058,
    0.0225990, 0.0221980, 0.0218027, 0.0214130, 0.0210286,
    0.0206495, 0.0202756, 0.0199068, 0.0195428, 0.0191837,
    0.0188293, 0.0184795, 0.0181343, 0.0177934, 0.0174569,
    0.0171246, 0.0167965, 0.0164724, 0.0161523, 0.0158361,
    0.0155236, 0.0152150, 0.0149100, 0.0146085, 0.0143107,
    0.0140162, 0.0137251, 0.0134374, 0.0131529, 0.0128716,
    0.0125935, 0.0123184, 0.0120463, 0.0117772, 0.0115110,
    0.0112477, 0.0109872, 0.0107294, 0.0104744, 0.0102220,
    0.0099722, 0.0097250, 0.0094804, 0.0092382, 0.0089984,
    0.0087611, 0.0085261, 0.0082934, 0.0080630, 0.0078349,
    0.0076089, 0.0073852, 0.0071636, 0.0069440, 0.0067266,
    0.0065112, 0.0062978, 0.0060864, 0.0058769, 0.0056693,
    0.0054636, 0.0052598, 0.0050578, 0.0048576, 0.0046592,
    0.0044625, 0.0042676, 0.0040743, 0.0038828, 0.0036928,
    0.0035045, 0.0033178, 0.0031327, 0.0029492, 0.0027671,
    0.0025866, 0.0024076, 0.0022301, 0.0020540, 0.0018793,
    0.0017061, 0.0015342, 0.0013637, 0.0011946, 0.0010268,
    0.0008603, 0.0006952, 0.0005313, 0.0003687, 0.0002074,
    0.0000473, -0.0001116, -0.0002693, -0.0004258, -0.0005811,
   -0.0007353, -0.0008883, -0.0010402, -0.0011909, -0.0013406,
   -0.0014892, -0.0016367, -0.0017831, -0.0019285, -0.0020728,
   -0.0022161, -0.0023584, -0.0024997, -0.0026400, -0.0027793,
   -0.0029176, -0.0030550, -0.0031915, -0.0033270, -0.0034616,
   -0.0035953, -0.0037280, -0.0038599, -0.0039909, -0.0041210,
   -0.0042503, -0.0043787, -0.0045063, -0.0046330, -0.0047589,
   -0.0048840, -0.0050083, -0.0051317, -0.0052544, -0.0053763,
   -0.0054975, -0.0056178, -0.0057374, -0.0058563, -0.0059744,
   -0.0060918, -0.0062084, -0.0063244, -0.0064396, -0.0065541,
   -0.0066679, -0.0067811, -0.0068935, -0.0070053, -0.0071164,
   -0.0072268, -0.0073366, -0.0074457, -0.0075542, -0.0076621,
   -0.0077693, -0.0078759, -0.0079819, -0.0080873, -0.0081920,
   -0.0082962, -0.0083998, -0.0085027, -0.0086051, -0.0087069,
   -0.0088082, -0.0089089, -0.0090090, -0.0091085, -0.0092075,
   -0.0093060, -0.0094039, -0.0095013, -0.0095981, -0.0096945,
   -0.0097903, -0.0098855, -0.0099803, -0.0100746, -0.0101683,
   -0.0102616, -0.0103544, -0.0104466, -0.0105384, -0.0106297,
   -0.0107206, -0.0108109, -0.0109008, -0.0109902, -0.0110792,
   -0.0111677, -0.0112558, -0.0113434, -0.0114305, -0.0115173,
   -0.0116035, -0.0116894, -0.0117748, -0.0118598, -0.0119444,
    5.07811e+17, -9.81600e-54
  };

  // ---------------------------------------------------------------------*
  //                                                                      *
  //  New version of DRES created  by A.Ferrari & P.Sala, INFN - Milan    *
  //                                                                      *
  //  Last change  on  10-apr-93  by  Alfredo Ferrari, INFN - Milan       *
  //                                                                      *
  //  Dres93: Dres91 plus the RAL fission model taken from LAHET thanks   *
  //          to R.E.Prael                                                *
  //  Dres91: new version from A. Ferrari and P. Sala, INFN - Milan       *
  //          This routine has been adapted from the original one of the  *
  //          Evap-5 module (KFA - Julich). Main modifications concern    *
  //          with kinematics which is now fully relativistic and with    *
  //          the treatment of few nucleons nuclei, which are now frag-   *
  //          mented, even though in a very rough manner. Changes have    *
  //          been made also to other routines of the Evap-5 package      *
  //                                                                      *
  // ---------------------------------------------------------------------*

  // ---------------------------------------------------------------------*
  //                                                                      *
  //  Input variables:                                                    *
  //     M2 = Mass number of the residual nucleus                         *
  //     M3 = Atomic number of the residual nucleus                       *
  //     T1 = Excitation energy of the residual nucleus before evaporation*
  //     U  = Excitation energy of the residual nucleus after evaporation *
  //     Erec = Recoil kinetic energy of the residual nucleus             *
  //            The recoil direction is given by Coslbr (i)               *
  //                                                                      *
  //  Significant variables:                                              *
  //     JA = Present mass number of the residual nucleus                 *
  //     JZ = Present atomic number of the residual nucleus               *
  //     Smom1 = Energy accumulators for the six types of evaporated      *
  //             particles                                                *
  //                                                                      *
  //    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
  //    !!!! Please note that the following variables concerning !!!!     *
  //    !!!! with the present residual nucleus must be set before!!!!     *
  //    !!!! entering DRES91: Ammres, Ptres                      !!!!     *
  //    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
  //                                                                      *
  // ---------------------------------------------------------------------*

  // ---------------------------------------------------------------------
  // SUBNAME = DRES --- EVAPORATION
  //           EVAPORATION DATA SHOUD BE READ ON INPUT STAGE
  // ---------------------------------------------------------------------
  // *****A LA EVAP III(TWA,8-68)

  ++ievevp;
  // -------------------------------------- 1.ST CALL INIT
  if (lfirst) {
    lfirst = false;
    pairCorrFlag = 0.;
    exmass[0] = 1.e3 * ( AMNEUT - AMUAMU );
    exmass[1] = GetIsotopicMass(1.0, 1.0);
    exmass[2] = GetIsotopicMass(2.0, 1.0);
    exmass[3] = GetIsotopicMass(3.0, 1.0);
    exmass[4] = GetIsotopicMass(3.0, 2.0);
    exmass[5] = GetIsotopicMass(4.0, 2.0);
    zmass[0] = exmass[0] + 1.0e3 * AMUAMU;
    zmass[1] = exmass[1] + 1.0e3 * AMUAMU;
    zmass[2] = exmass[2] + 2.0e3 * AMUAMU;
    zmass[3] = exmass[3] + 3.0e3 * AMUAMU;
    zmass[4] = exmass[4] + 3.0e3 * AMUAMU;
    zmass[5] = exmass[5] + 4.0e3 * AMUAMU;
    bnmass[0] = 0.;
    bnmass[1] = 0.;
    bnmass[2] = zmass[0] + zmass[1] - zmass[2];
    bnmass[3] = zmass[0] * 2. + zmass[1] - zmass[3];
    bnmass[4] = zmass[0] + zmass[1] * 2. - zmass[4];
    bnmass[5] = (zmass[0] + zmass[1]) * 2. - zmass[5];
    for (kk = 1; kk <= 6; ++kk) {
      z2mass[kk - 1] = zmass[kk - 1] * zmass[kk - 1];
    }
    emn = 1.0e3 * AMNEUT;
    emh = zmass[1];
    um = AMUMEV + GetIsotopicMass(16.0, 8.0) / 16.;
    emhn = emh - emn;
    emnum = emn - um;
  }
  //  End of initialization:
  //
  //     --------------------------------- START OF PROCESS
  //  Initialize nVariousFragm and smom1 if nothing has been already
  //  evaporated for this event
  for (i = 1; i <= 6; ++i) {
    nVariousFragm[i - 1] = 0;
    smom1[i - 1] = 0.;
  }
  ja = m2;
  jz = m3;
  *u = t1;
  rnmass = massResidNucl * 1e3 + *u;
  // P2res and  Ptres are the squared momentum and the momentum of the
  // residual nucleus (now in relativistic kinematics), Umo the
  // invariant mass of the system!
  umo = rnmass;
  umo2 = umo * umo;
  elbtot = rnmass + *erec;
  gamcm = elbtot / rnmass;
  etacm = totMomResidNucl * 1e3 / rnmass;
L1000:
  *loppar = false;
  // Check for starting data inconsistencies
  if (ja - jz < 0) {
    G4cout << " Dres: cascade residual nucleus has mass no. less than Z!!"
	 << G4endl;
    return;
  } else if (ja <= 6 || jz <= 2) {
    //
    //                Rough treatment for very few nucleon residual
    //                nuclei. The basic ideas are:
    //        a) as many as possible alpha particles are emitted
    //        b) particles are emitted one per time leaving a residual
    //           excitation energy proportional to number of nucleons
    //           left in the residual nucleus (so we deal only with
    //           two body kinematics)
    //       T A K E   I N T O   A C C O U N T   T H A T   T H I S
    //       T R E A T M E N T   I S   E X T R E M E L Y   R O U G H
    //       T H E   T A S K   B E I N G   O N L Y   T O   S U P P L Y
    //       S O M E T H I N G   T O   S H A R E   E N E R G Y   A N D
    //       M O M E N T U M   A M O N G   A   F E W   F R A G M E N T S
    jresid = 0;
    // First check we are not concerning with a couple of neutrons or
    // protons
    if (ja == 2 && jz != 1) {
      jemiss = jz / 2 + 1;
      jresid = jemiss;
      rnmass = zmass[jresid - 1];
      *u = 0.;
      deltu = umo - zmass[jemiss - 1] * 2.;
      if (deltu <= 0.) {
	if (deltu < umo * -2.0 * ANGLGB) {
	  G4cout << " *** Dres: insufficient Umo for a nucleon couple"
	       << " " << umo << " " << zmass[jemiss - 1] * 2. << G4endl;
	}
	umo = (umo + deltu) * ( 1.0 + ANGLGB );
      }
      goto L2500;
    }
    // Then check we are not concerning with one of the six
    // standard particles
    for (j = 6; j >= 1; --j) {
      if (jz == iz[j - 1] && ja == ia[j - 1]) {
	switch (j) {
	case 1:
	case 2:
	  //  Proton or neutron, nothing can be Done
	  return;
	case 3:
	  //  Deuteron:
	  // Split into a proton and a neutron if it is possible
	  if (*u > bnmass[2]) {
	    jemiss = 1;
	    jresid = 2;
	    rnmass = zmass[1];
	    *u = 0.;
	    goto L2500;
	  } else {
	    // Energy too low to split the deuteron, return
	    G4cout << " **Dres: energy too low to split a deuteron! M2,M3"
		 << " " << m2 << " " << m3 << G4endl;
	    return;
	  }
	case 4:
	  //  Triton:
	  d__1 = 0., d__2 = *u + bnmass[2] - bnmass[3];
	  deuneu = G4std::max(d__1,d__2);
	  d__1 = 0., d__2 = *u - bnmass[3];
	  prnene = G4std::max(d__1,d__2);
	  qnorm = deuneu + prnene;
	  //  If we cannot split then return
	  if (qnorm <= 0.) {
	    return;
	  }
	  rndm[0] = G4UniformRand();
	  v = rndm[0];
	  // Split or into a deuteron and a neutron
	  // or into two protons and one neutron: no account is
	  // made for Coulomb effects, probability is simply assumed
	  // proportional to reaction Qs
	  if (v < deuneu / qnorm) {
	    // A deuteron and a proton selected
	    jemiss = 1;
	    jresid = 3;
	    rnmass = zmass[2];
	    *u = 0.;
	    goto L2500;
	    // Split into 1 proton and 2 neutrons: part of the exci-
	    // tation energy is conserved to allow the further
	    // splitting of the deuteron
	  } else {
	    jemiss = 1;
	    jresid = 0;
	    fact = 1.;
	    // Loop to compute the residual excitation energy
	    do {
	      fact *= .6666666666666667;
	      // Erncm, Eepcm are the total energies of the residual
	      // nucleus and of the emitted particle in the CMS frame
	      *u = fact * prnene + bnmass[2];
	      rnmass = zmass[2] + *u;
	      d__1 = rnmass;
	      erncm = (umo2 + d__1 * d__1 - z2mass[jemiss - 1]) * .5 / umo;
	      eepcm = umo - erncm;
	    } while (eepcm <= zmass[jemiss - 1]);
	    goto L2600;
	  }
	case 5:
	  //  3-He:
	  d__1 = 0., d__2 = *u + bnmass[2] - bnmass[4];
	  deupro = G4std::max(d__1,d__2);
	  d__1 = 0., d__2 = *u - bnmass[4];
	  prprne = G4std::max(d__1,d__2);
	  qnorm = deupro + prprne;
	  //  If we cannot split then return
	  if (qnorm <= 0.) {
	    return;
	  }
	  rndm[0] = G4UniformRand();
	  v = rndm[0];
	  // Split or into a deuteron and a proton
	  // or into two protons and one neutron: no account is
	  // made for Coulomb effects, probability is simply assumed
	  // prportional to reaction Qs
	  if (v < deupro / qnorm) {
	    // A deuteron and a proton selected
	    jemiss = 2;
	    jresid = 3;
	    rnmass = zmass[2];
	    *u = 0.;
	    goto L2500;
	    // Split into 2 protons and 1 neutron: part of the exci-
	    // tation energy is conserved to allow the further
	    // splitting of the deuteron
	  } else {
	    jemiss = 2;
	    jresid = 0;
	    fact = 1.;
	    // Loop to compute the residual excitation energy
	    do {
	      fact *= .6666666666666667;
	      // Erncm, Eepcm are the total energies of the residual
	      // nucleus and of the emitted particle in the CMS frame
	      *u = fact * prprne + bnmass[2];
	      rnmass = zmass[2] + *u;
	      d__1 = rnmass;
	      erncm = (umo2 + d__1 * d__1 - z2mass[jemiss - 1]) * .5 / umo;
	      eepcm = umo - erncm;
	    } while (eepcm <= zmass[jemiss - 1]);
	    goto L2600;
	  }
	case 6:
	  //  Alpha:
	  //
	  d__1 = 0., d__2 = *u + bnmass[2] * 2. - bnmass[5];
	  deudeu = G4std::max(d__1,d__2);
	  d__1 = 0., d__2 = *u + bnmass[3] - bnmass[5];
	  protri = G4std::max(d__1,d__2);
	  d__1 = 0., d__2 = *u + bnmass[4] - bnmass[5];
	  ueu3he = G4std::max(d__1,d__2);
	  qnorm = deudeu + protri + ueu3he;
	  //  If we cannot split then return
	  if (qnorm <= 0.) {
	    return;
	  }
	  rndm[0] = G4UniformRand();
	  v = rndm[0];
	  // Split or into two deuterons or a triton and a proton
	  // or a 3-He and a neutron: no account is made for
	  // Coulomb effects, probability is simply assumed
	  // proportional to reaction Qs
	  if (v < deudeu / qnorm) {
	    // Two deuterons selected
	    jemiss = 3;
	    jresid = 3;
	    rnmass = zmass[2];
	    *u = 0.;
	    goto L2500;
	  } else if (v < (deudeu + protri) / qnorm) {
	    // Split into a triton and a proton
	    jemiss = 2;
	    jresid = 4;
	    rnmass = zmass[3];
	    *u = 0.;
	    goto L2500;
	  } else {
	    // Split into a 3-He and a neutron
	    jemiss = 1;
	    jresid = 5;
	    rnmass = zmass[4];
	    *u = 0.;
	    goto L2500;
	  }
	}
      }
    }
    a = (G4double) ja;
    z = (G4double) jz;
    q[0] = 0.;
    energ0 = GetIsotopicMass(a, z);
    //   Note that Q(i) are not the reaction Qs but the remaining
    //   energy after the reaction
    for (k = 1; k <= 6; ++k) {
      jja = ja - ia[k - 1];
      jjz = jz - iz[k - 1];
      jjn = jja - jjz;
      if (jjn < 0 || jjz < 0) {
	q[k] = q[k - 1];
	continue;
      }
      ddjja = (G4double) jja;
      ddjjz = (G4double) jjz;
      d__1 = *u + energ0 - GetIsotopicMass(ddjja, ddjjz) - exmass[k - 1];
      q[k] = G4std::max(d__1,0.) + q[k - 1];
    }
    //  If no emission channel is open then return
    if (q[6] <= 0.) {
      return;
    }
    rndm[0] = G4UniformRand();
    v = rndm[0];
    fact = 1.;
    for (j = 1; j <= 6; ++j) {
      if (v < q[j] / q[6]) {
	jemiss = j;
	jja = ja - ia[jemiss - 1];
	jjz = jz - iz[jemiss - 1];
	for (jj = 1; jj <= 6; ++jj) {
	  if (jja == ia[jj - 1] && jjz == iz[jj - 1]) {
	    jresid = jj;
	    rnmass = zmass[jresid - 1];
	    erncm = (umo2 + z2mass[jresid - 1] - z2mass[jemiss - 1]) *
	      .5 / umo;
	    eepcm = umo - erncm;
	    *u = 0.;
	    goto L2600;
	  }
	}
	ajja = (G4double) jja;
	zjjz = (G4double) jjz;
	rnmas0 = ajja * AMUMEV + GetIsotopicMass(ajja, zjjz);
	goto L2300;
      }
    }
    G4cout << " **** error in Dres, few nucleon treatment ****" << G4endl;
    return;
    // Loop to compute the residual excitation energy
L2300:
    fact = fact * ajja / a;
    *u = fact * (q[jemiss] - q[jemiss - 1]);
    // Erncm, Eepcm are the total energies of the residual
    // nucleus and of the emitted particle in the CMS frame
    rnmass = rnmas0 + *u;
    d__1 = rnmass;
    erncm = (umo2 + d__1 * d__1 - z2mass[jemiss - 1]) * .5 / umo;
    eepcm = umo - erncm;
    if (eepcm <= zmass[jemiss - 1]) {
      if (q[jemiss] - q[jemiss - 1] >= 1e-6) {
	goto L2300;
      }
      //  Actually there is no excitation energy available!
      *u = ANGLGB;
      rnmass = rnmas0 * ONEPLS;
      erncm = rnmass * ONEPLS;
      eepcm = zmass[jemiss - 1] * ONEPLS;
    }
    goto L2600;
    //  From here standard two bodies kinematics with Jemiss, Rnmass
    //  (new excitation energy is U)
L2500:
    //  Erncm, Eepcm are the total energies of the residual
    //  nucleus and of the emitted particle in the CMS frame
    d__1 = rnmass;
    erncm = (umo2 + d__1 * d__1 - z2mass[jemiss - 1]) * .5 / umo;
    eepcm = umo - erncm;
L2600:
    //  C(i) are the direction cosines of the emitted particle
    //  (Jemiss) in the CMS frame, of course - C(i)
    //  are the ones of the residual nucleus (Jresid if one of the
    //  standard partcles, say the proton)
    RanDirCos(c, &c[1], &c[2]);
    d__1 = eepcm;
    pcms = sqrt(d__1 * d__1 - z2mass[jemiss - 1]);
    //  Now we perform the Lorentz transformation back to the original
    //  frame (lab frame)
    //  First the emitted particle:
    etax = etacm * dirCosXResidNucl;
    etay = etacm * dirCosYResidNucl;
    etaz = etacm * dirCosZResidNucl;
    pcmsx = pcms * c[0];
    pcmsy = pcms * c[1];
    pcmsz = pcms * c[2];
    etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
    eps = gamcm * eepcm + etapcm - zmass[jemiss - 1];
    phelp = etapcm / (gamcm + 1.) + eepcm;
    plbpx = pcmsx + etax * phelp;
    plbpy = pcmsy + etay * phelp;
    plbpz = pcmsz + etaz * phelp;
    phelp = sqrt(plbpx * plbpx + plbpy * plbpy + plbpz * plbpz);
    coslbp[0] = plbpx / phelp;
    coslbp[1] = plbpy / phelp;
    coslbp[2] = plbpz / phelp;
    //  Then the residual nucleus ( for it c (i) --> - c (i) ):
    *erec = gamcm * erncm - etapcm - rnmass;
    kinEnResidNucl = *erec * .001;
    phelp = -etapcm / (gamcm + 1.) + erncm;
    momXResidNucl = (-pcmsx + etax * phelp) * .001;
    momYResidNucl = (-pcmsy + etay * phelp) * .001;
    momZResidNucl = (-pcmsz + etaz * phelp) * .001;
    p2res = momXResidNucl * momXResidNucl +
      momYResidNucl * momYResidNucl + momZResidNucl * momZResidNucl;
    totMomResidNucl = sqrt(p2res);
    dirCosXResidNucl = momXResidNucl / totMomResidNucl;
    dirCosYResidNucl = momYResidNucl / totMomResidNucl;
    dirCosZResidNucl = momZResidNucl / totMomResidNucl;
    //  Score the emitted particle
    ++nVariousFragm[jemiss - 1];
    smom1[jemiss - 1] += eps;
    itemp = nVariousFragm[jemiss - 1] - 1 + (jemiss - 1)*MXEVAP;
    Evaporates[itemp].SetZero();
    Evaporates[itemp].SetMass( massFragm[jemiss - 1] );
    Evaporates[itemp].SetMomentum( coslbp[0], coslbp[1], coslbp[2] );
    Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
    Evaporates[itemp].SetParticleDef( pdefFragm[jemiss - 1] );
    //  Check if the residual nucleus is one of the emitted particles
    if (jresid > 0) {
      j = jresid;
      if (*u > ANGLGB) {
	goto L1000;
      }
      jemiss = j;
      // *****STORE,RESIDUAL NUC IS OF EMITTED PARTICLE TYPE
      // If we are here this means that the residual nucleus is equal to
      // one of the six emitted particle (the j-th one). So give to it
      // all the energy, score it and return with 0 recoil and excitation
      // energy for the residual nucleus
      eps = *erec;
      ++nVariousFragm[jemiss - 1];
      smom1[jemiss - 1] += eps;
      itemp = nVariousFragm[jemiss - 1] - 1 + (jemiss - 1)*MXEVAP;
      Evaporates[itemp].SetZero();
      Evaporates[itemp].SetMass( massFragm[jemiss - 1] );
      Evaporates[itemp].SetMomentum( dirCosXResidNucl,
				     dirCosYResidNucl,
				     dirCosZResidNucl );
      Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
      Evaporates[itemp].SetParticleDef( pdefFragm[jemiss - 1] );
      *loppar = false;
      *erec = 0.;
      *u = 0.;
      kinEnResidNucl = 0.;
      totMomResidNucl = 0.;
      return;
      // -->-->-->-->--> go to return
    }
    ja -= ia[jemiss - 1];
    jz -= iz[jemiss - 1];
    // Umo is the invariant mass of the system!!
    umo = rnmass;
    umo2 = umo * umo;
    elbtot = rnmass + *erec;
    gamcm = elbtot / rnmass;
    etacm = totMomResidNucl * 1e3 / rnmass;
    goto L1000;
  }
  // Come here at the beginning and after the end of a "normal"
  // evaporation cycle
  do {
    a = (G4double) ja;
    z = (G4double) jz;
    if (ja - 8 == 0 && jz - 4 == 0) {
      break;
    }
    for (k = 1; k <= 6; ++k) {
      if (a - fla[k - 1] <= z - flz[k - 1]) {
	q[k] = 99999.;
	continue;
      }
      if (a < fla[k - 1] * 2. || z < flz[k - 1] * 2.) {
	continue;
      }
      q[k] =
	NuclearBindingEnergy(a - fla[k - 1],
			     z - flz[k - 1], a, z) + exmass[k - 1];
    }
    flkcou[0] = 0.;
    flkcou[1] = CoulombBarrier(1, z - flz[1]);
    flkcou[2] = flkcou[1] + .06;
    flkcou[3] = flkcou[1] + .12;
    flkcou[5] = CoulombBarrier(2, z - flz[5]);
    flkcou[4] = flkcou[5] - .06;
    ccoul[0] = 1.;
    ccou2 = CoulombBarrier(3, z - flz[1]);
    ccoul[1] = ccou2 + 1.;
    ccoul[2] = ccou2 * 1.5 + 3.;
    ccoul[3] = ccou2 + 3.;
    ccoul[5] = CoulombBarrier(4, z - flz[5]) * 2. + 2.;
    ccoul[4] = ccoul[5] * 2. - 1.;
    sigma = 0.;
    //  Initialize the flag which checks for open particle decay with
    //  zero excitation and pairing --> for particle unstable residual
    //  nuclei
    *loppar = false;
    for (j = 1; j <= 6; ++j) {
      if (a < fla[j - 1] * 2. || z < flz[j - 1] * 2. ) {
	r[j - 1] = 0.;
	s[j - 1] = 0.;
	sos[j - 1] = 0.;
	continue;
      }
      mm = ja - ia[j - 1];
      zz = z - flz[j - 1];
      aa = a - fla[j - 1];
      if (aa <= zz) {
	r[j - 1] = 0.;
	s[j - 1] = 0.;
	sos[j - 1] = 0.;
	continue;
      }
      //  Energy threshold for the emission of the jth-particle
      thresh[j - 1] =
	q[j] + flkcou[j - 1] * .88235 * flz[j - 1] * zz /
	(residualNuclearMasses[mm - 1] + rho[j - 1]);
      *loppar = *loppar || thresh[j - 1] > 0.;
      iaa = NINT(aa);
      izz = NINT(zz);
      nn = iaa - izz;
      //  The residual nucleus excitation energy ranges from 0 up
      //  to U - Q (J)
      umxres = *u - thresh[j - 1];
      //  This is the a lower case of the level density
      smalla[j - 1] = LevelDensity(izz, nn, &asmmax, &asmmin);
      ExcitationEnergyLevel(iaa, izz, &eex1st, &eex2nd, &corr);
      eex1st *= 1e3;
      eex2nd *= 1e3;
      corr = G4std::max(corr,0.) * 1e3;
      if (nn == 4 && izz == 4) {
	if (*u - thresh[j - 1] - 6.1 > 0.) {
	  corr = 6.;
	} else {
	  tmpvar = *u - thresh[j - 1] - .1;
	  corr = G4std::max(0.,tmpvar);
	}
      }
      if (NINT(pairCorrFlag) == 1) {
	corr = 0.;
      }
      corrrr[j - 1] = corr;
      // Standard calculation:
      arg = *u - thresh[j - 1] - corr;
      if (arg < 0.) {
	r[j - 1] = 0.;
	s[j - 1] = 0.;
	sos[j - 1] = 0.;
	continue;
      }
      s[j - 1] = sqrt(smalla[j - 1] * arg) * 2.;
      sos[j - 1] = s[j - 1] * 10.;
    }
    n1 = 1;
    d__1 =
      G4std::max(s[0],s[1]), d__1 = G4std::max(d__1,s[2]), d__1 =
      G4std::max(d__1,s[3]), d__1 = G4std::max(d__1,s[4]);
    ses = G4std::max(d__1,s[5]);
    for (j = 1; j <= 6; ++j) {
      js = (G4int) (sos[j - 1] + 1.);
      fjs = (G4double) js;
      strun[j - 1] = fjs - 1.;
      if (s[j - 1] > 0.) {
	mm = ja - ia[j - 1];
	d__3 = EXPMIN, d__4 = s[j - 1] - ses;
	d__1 = EXPMAX, d__2 = G4std::max(d__3,d__4);
	expsas = G4std::min(d__1,d__2);
	sas = exp(expsas);
	d__3 = EXPMIN, d__4 = -s[j - 1];
	d__1 = EXPMAX, d__2 = G4std::max(d__3,d__4);
	expsus = G4std::min(d__1,d__2);
	sus = exp(expsus);
	d__1 = s[j - 1];
	d__2 = s[j - 1];
	d__3 = smalla[j - 1];
	eye1[j - 1] =
	  (d__1 * d__1 * 2. - s[j - 1] * 6. + 6. + sus * (d__2 * d__2 - 6.)) /
	  (d__3 * d__3 * 8.);
	if (j == 1) {
	  eye0[j - 1] = (s[j - 1] - 1. + sus) / (smalla[j - 1] * 2.);
	  // Standard calculation
	  d__1 = residualNuclearMasses[mm - 1];
	  r[j - 1] =
	    d__1 * d__1 * alph[mm - 1] *
	    (eye1[j - 1] + bet[mm - 1] * eye0[j - 1]) * sas;
	} else {
	  d__1 = residualNuclearMasses[mm - 1];
	  r[j - 1] = ccoul[j - 1] * (d__1 * d__1) * eye1[j - 1] * sas;
	}
	d__1 = 0., d__2 = r[j - 1];
	r[j - 1] = G4std::max(d__1,d__2);
	sigma += r[j - 1];
      }
    }
    ncount = 0;
L6202:
    if (sigma <= 0.) {
      goto L9;
    } else {
      goto L10;
    }
L9:
    for (j = 1; j <= 6; ++j) {
      if (ja - ia[j - 1] == 0 && jz - iz[j - 1] == 0) {
	if (*u > ANGLGB) {
	  goto L1000;
	}
	jemiss = j;
	// *****STORE,RESIDUAL NUC IS OF EMITTED PARTICLE TYPE
	// If we are here this means that the residual nucleus is equal to
	// one of the six emitted particle (the j-th one). So give to it
	// all the energy, score it and return with 0 recoil and excitation
	// energy for the residual nucleus
	eps = *erec;
	++nVariousFragm[jemiss - 1];
	smom1[jemiss - 1] += eps;
	itemp = nVariousFragm[jemiss - 1] - 1 + (jemiss - 1)*MXEVAP;
	Evaporates[itemp].SetZero();
	Evaporates[itemp].SetMass( massFragm[jemiss - 1] );
	Evaporates[itemp].SetMomentum( dirCosXResidNucl,
				       dirCosYResidNucl,
				       dirCosZResidNucl );
	Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
	Evaporates[itemp].SetParticleDef( pdefFragm[jemiss - 1] );
	*loppar = false;
	*erec = 0.;
	*u = 0.;
	kinEnResidNucl = 0.;
	totMomResidNucl = 0.;
	return;
	// -->-->-->-->--> go to return
      }
    }
    if (ja - 8 != 0 || jz - 4 != 0) {
      return;
    }
    break;
    // Come here for a "normal" Step
L10:
    *loppar = false;
    rndm[0] = G4UniformRand();
    uran = rndm[0] * sigma;
    sum = 0.;
    for (j = 1; j <= 6; ++j) {
      k = j;
      sum = r[j - 1] + sum;
      if (sum - uran > 0.) {
	break;
      }
    }
    jemiss = k;
    ++nVariousFragm[jemiss - 1];
    js = (G4int) (sos[jemiss - 1] + 1.);
    if (js - 1000 >= 0) {
      d__1 = s[jemiss - 1], d__2 = d__1;
      d__3 = s[jemiss - 1];
      d__4 = s[jemiss - 1];
      ratio2 =
	(d__2 * (d__1 * d__1) - d__3 * d__3 * 6. + s[jemiss - 1] * 15. - 15.) /
	((d__4 * d__4 * 2. - s[jemiss - 1] * 6. + 6.) * smalla[jemiss - 1]);
    } else {
      ratio2 = (p2[js - 1] + (p2[js] - p2[js - 1]) *
		(sos[jemiss - 1] - strun[jemiss - 1])) / smalla[jemiss - 1];
    }
    epsav = ratio2 * 2.;
    //  Neutron channel selected:
    if (jemiss == 1) {
      mm = ja - ia[j - 1];
      epsav =
	(epsav + bet[mm - 1]) /
	(bet[mm - 1] * eye0[jemiss - 1] / eye1[jemiss - 1] + 1.);
    }
    do {
      rndm[0] = G4UniformRand();
      rndm[1] = G4UniformRand();
      e1 = log(rndm[0]) * -.5;
      e2 = log(rndm[1]) * -.5;
      // Eps should be the total kinetic energy in the CMS frame
      // Standard calculation:
      eps = (e1 + e2) * epsav + thresh[jemiss - 1] - q[jemiss];
      ar = a - ia[jemiss - 1];
      zr = z - iz[jemiss - 1];
      // The CMS energy is updated
      imass = NINT(ar);
      if (imass == 8 && NINT(zr) == 4) {
	unew = *u - eps - q[jemiss];
	umax = *u - thresh[jemiss - 1];
	if (unew > 6.) {
	  umin = 6.;
	} else if (unew > 4.47 && umax > 6.) {
	  umin = 4.47;
	  unew = 6.;
	} else if (unew > 1.47 && umax > 2.94) {
	  umin = 1.47;
	  unew = 2.94;
	} else {
	  umin = -.1;
	  unew = ANGLGB * 0.1;
	}
      } else if (imass <= 4) {
	ipro = NINT(zr);
	ineu = imass - ipro;
	if (imass == 1) {
	  //  Be sure that residual neutrons or protons are not left excited
	  umin = 0.;
	  unew = 0.;
	  eps = *u - q[jemiss];
	} else if (ipro == 0 || ineu == 0) {
	  //  Ipro protons or ineu neutrons arrived here!
	  umin = corrrr[jemiss - 1];
	  unew = *u - eps - q[jemiss];
	} else if (imass <= 2) {
	  //  Be sure that residual deuterons are not left excited!
	  umin = 0.;
	  unew = 0.;
	  eps = *u - q[jemiss];
	} else if ((i__1 = ineu - ipro, abs(i__1)) <= 1) {
	  //  For the moment also residual 3-H, 3-He and 4-He are not left
	  //  excited !
	  umin = 0.;
	  unew = 0.;
	  eps = *u - q[jemiss];
	} else {
	  umin = corrrr[jemiss - 1];
	  unew = *u - eps - q[jemiss];
	}
      } else {
	umin = corrrr[jemiss - 1];
	unew = *u - eps - q[jemiss];
      }
      // Standard calculation
      if (unew - umin >= 0.) {
	rnmass = ar * AMUMEV + GetIsotopicMass(ar, zr) + unew;
	d__1 = rnmass + zmass[jemiss - 1];
	umin2 = d__1 * d__1;
	if (umin2 >= umo2) {
	  ++ncount;
	  continue;
	}
	*u = unew;
	// C(i) are the direction cosines of the evaporated particle in the CMS
	// frame, of course - C(i) are the ones of the residual nucleus
	RanDirCos(c, &c[1], &c[2]);
	// Erncm, Eepcm are the total energies of the residual nucleus and
	// of the evaporated particle in the CMS frame
	d__1 = rnmass;
	erncm = (umo2 + d__1 * d__1 - z2mass[jemiss - 1]) * .5 / umo;
	eepcm = umo - erncm;
	d__1 = eepcm;
	pcms = sqrt(d__1 * d__1 - z2mass[jemiss - 1]);
	// Now we perform the Lorentz transformation back to the original
	// frame (lab frame)
	// First the evaporated particle:
	etax = etacm * dirCosXResidNucl;
	etay = etacm * dirCosYResidNucl;
	etaz = etacm * dirCosZResidNucl;
	pcmsx = pcms * c[0];
	pcmsy = pcms * c[1];
	pcmsz = pcms * c[2];
	etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
	eps = gamcm * eepcm + etapcm - zmass[jemiss - 1];
	phelp = etapcm / (gamcm + 1.) + eepcm;
	plbpx = pcmsx + etax * phelp;
	plbpy = pcmsy + etay * phelp;
	plbpz = pcmsz + etaz * phelp;
	phelp = sqrt(plbpx * plbpx + plbpy * plbpy + plbpz * plbpz);
	coslbp[0] = plbpx / phelp;
	coslbp[1] = plbpy / phelp;
	coslbp[2] = plbpz / phelp;
	// Then the residual nucleus ( for it c (i) --> - c (i) ):
	*erec = gamcm * erncm - etapcm - rnmass;
	kinEnResidNucl = *erec * .001;
	phelp = -etapcm / (gamcm + 1.) + erncm;
	momXResidNucl = (-pcmsx + etax * phelp) * .001;
	momYResidNucl = (-pcmsy + etay * phelp) * .001;
	momZResidNucl = (-pcmsz + etaz * phelp) * .001;
	p2res =
	  momXResidNucl * momXResidNucl + momYResidNucl * momYResidNucl +
	  momZResidNucl * momZResidNucl;
	totMomResidNucl = sqrt(p2res);
	dirCosXResidNucl = momXResidNucl / totMomResidNucl;
	dirCosYResidNucl = momYResidNucl / totMomResidNucl;
	dirCosZResidNucl = momZResidNucl / totMomResidNucl;
	// Check energy and momentum conservation !!
	if (*erec <= 0.) {
	  totMomResidNucl = 0.;
	  *erec = 0.;
	}
	// Umo is the invariant mass of the system!!
	umo = rnmass;
	umo2 = umo * umo;
	elbtot = rnmass + *erec;
	gamcm = elbtot / rnmass;
	etacm = totMomResidNucl * 1e3 / rnmass;
	goto L76;
      }
      ++ncount;
    } while (ncount <= 10);
    sigma -= r[jemiss - 1];
    // if we are here we have sampled for > 10 Times a negative energy Unew
    --nVariousFragm[jemiss - 1];
    r[jemiss - 1] = 0.;
    ncount = 0;
    goto L6202;
L76:
    jat = ja - ia[jemiss - 1];
    jzt = jz - iz[jemiss - 1];
    if (jat - jzt <= 0) {
      goto L9;
    }
    ja = jat;
    jz = jzt;
    // *****STORE,END OF NORMAL CYCLE
    smom1[jemiss - 1] += eps;
    itemp = nVariousFragm[jemiss - 1] - 1 + (jemiss - 1)*MXEVAP;
    Evaporates[itemp].SetZero();
    Evaporates[itemp].SetMass( massFragm[jemiss - 1] );
    Evaporates[itemp].SetMomentum( coslbp[0], coslbp[1], coslbp[2] );
    Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
    Evaporates[itemp].SetParticleDef( pdefFragm[jemiss - 1] );
    // The following card switch to the rough splitting treatment
    if (ja <= 2) {
      goto L1000;
    }
  } while (ja - 8 != 0 || jz - 4 != 0);
  // If we are here the residual nucleus is a 8-Be one, break it into
  // two alphas with all the available energy (U plus the Q of the breakup)
  // , score them and return with 0 recoil and excitation energy
  *loppar = false;
  if (*u < 0.) {
    eps = 0.;
  }
  // C(i) are the direction cosines of the first alpha in the CMS
  // frame, of course - C(i) are the ones of the other
  RanDirCos(c, &c[1], &c[2]);
  // Ecms is the total energy of the alphas in the CMS frame
  ecms = umo * .5;
  d__1 = ecms;
  pcms = sqrt(d__1 * d__1 - z2mass[5]);
  // Now we perform the Lorentz transformation back to the original
  // frame (lab frame)
  // First alpha:
  etax = etacm * dirCosXResidNucl;
  etay = etacm * dirCosYResidNucl;
  etaz = etacm * dirCosZResidNucl;
  pcmsx = pcms * c[0];
  pcmsy = pcms * c[1];
  pcmsz = pcms * c[2];
  etapcm = pcmsx * etax + pcmsy * etay + pcmsz * etaz;
  eps = gamcm * ecms + etapcm - zmass[5];
  phelp = etapcm / (gamcm + 1.) + ecms;
  plbpx = pcmsx + etax * phelp;
  plbpy = pcmsy + etay * phelp;
  plbpz = pcmsz + etaz * phelp;
  phelp = sqrt(plbpx * plbpx + plbpy * plbpy + plbpz * plbpz);
  // Store the first alpha!!
  smom1[5] += eps;
  ++nVariousFragm[5];
  itemp = nVariousFragm[5] - 1 + 499*MXEVAP;
  Evaporates[itemp].SetZero();
  Evaporates[itemp].SetMass( massFragm[5] );
  Evaporates[itemp].SetMomentum( plbpx/phelp, plbpy/phelp, plbpz/phelp );
  Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
  Evaporates[itemp].SetParticleDef( pdefFragm[5] );
  // Then the second alpha ( for it c (i) --> - c (i) ):
  eps = gamcm * ecms - etapcm - zmass[5];
  phelp = -etapcm / (gamcm + 1.) + ecms;
  plbpx = -pcmsx + etax * phelp;
  plbpy = -pcmsy + etay * phelp;
  plbpz = -pcmsz + etaz * phelp;
  phelp = sqrt(plbpx * plbpx + plbpy * plbpy + plbpz * plbpz);
  // Store the second alpha !!
  smom1[5] += eps;
  ++nVariousFragm[5];
  itemp = nVariousFragm[5] - 1 + 499*MXEVAP;
  Evaporates[itemp].SetZero();
  Evaporates[itemp].SetMass( massFragm[5] );
  Evaporates[itemp].SetMomentum( plbpx/phelp, plbpy/phelp, plbpz/phelp );
  Evaporates[itemp].SetKineticEnergyAndUpdate( eps / 1e3 );
  Evaporates[itemp].SetParticleDef( pdefFragm[5] );
  *loppar = false;
  *erec = 0.;
  *u = 0.;
  kinEnResidNucl = 0.;
  totMomResidNucl = 0.;
  return;

} // ResidualNucleusCascade


G4double G4MuonMinusCaptureAtRest::GetIsotopicMass(G4double a, G4double z)
{
  // System generated locals
  G4double ret_val;

  // Local variables
  static G4int n, ka0, iz0, kz0, izz;

  ka0 = NINT(a);
  kz0 = NINT(z);
  n = ka0 - kz0;
  if (n <= 0) {
    if (ka0 != 1) {
      if (n < 0) {
	G4cout << " Stopped in energy: mass number =< atomic number !!"
	     << " " << ka0 << " " << kz0 << G4endl;
      }
    } else {
      ret_val = isotopicData[2][0];
      return ret_val;
    }
  }
  izz = NINT(isotopicData[0][ka0 - 1]);
  if (kz0 < izz || kz0 > izz + 9) {
    ret_val = CalculateIsotopicMass(a, z);
  } else {
    iz0 = kz0 - izz + 2;
    ret_val = isotopicData[iz0 - 1][ka0 - 1];
    if (ret_val == 0. && (ka0 != 12 || kz0 != 6)) {
      ret_val = CalculateIsotopicMass(a, z);
    }
  }
  return ret_val;

} // GetIsotopicMass


void G4MuonMinusCaptureAtRest::Erup()
{
  // Local variables
  static G4int i, m2, m3;
  static G4double fpart[6];
  static G4bool loppar;
  static G4double fpartt;

  // ---------------------------------------------------------------------*
  //                                                                      *
  //     Created  on   15 may 1990     by     Alfredo & Paola Sala        *
  //                                              INFN - Milan            *
  //     Last change  on   10-apr-93   by     Alfredo Ferrari, INFN-Milan *
  //                                                                      *
  //     Derived from the ERUP routine of EVAP-V, HERMES, KFA-Julich      *
  //                                                                      *
  // ---------------------------------------------------------------------*

  // *****MODIFIED TO OBTAIN APR,ZPR AFTER CAS + EVAP (8-68,T.W.A.)

  //  Check the excitation energy
  if (iniExcitEnEvapNucl <= 0.) {
    //  No excitation energy:
    for (i = 1; i <= 6; ++i) {
      nVariousFragm[i - 1] = 0;
    }
    finExcitEnEvapNucl = iniExcitEnEvapNucl;
  } else {
    //  Positive excitation energy:
    //  Try evaporation
    m2 = NINT(atMassEvapNucl);
    m3 = NINT(chargeEvapNucl);
    for (;;) {
      ResidualNucleusCascade(m2, m3, iniExcitEnEvapNucl, &finExcitEnEvapNucl,
			     &recoilEnEvapNucl, &loppar);
      fpartt = 0.;
      for (i = 1; i <= 6; ++i) {
	fpart[i - 1] = (G4double) nVariousFragm[i - 1];
	fpartt += fpart[i - 1];
      }
      if (fpartt + pairCorrFlag >= ANGLGB || ! loppar) {
	break;
      }
      pairCorrFlag = 1.;
    }
    //  No more particles evaporated and pairing corrections accounted for
    pairCorrFlag = 0.;
    chargeEvapNucl =
      chargeEvapNucl - fpart[1] - fpart[2] -
      (fpart[4] + fpart[5]) * 2. - fpart[3];
    atMassEvapNucl =
      atMassEvapNucl - fpart[0] - fpart[1] -
      fpart[2] * 2. - (fpart[3] + fpart[4]) * 3. - fpart[5] * 4.;
  }
  return;

} // Erup


void G4MuonMinusCaptureAtRest::InitializeMuCapture()
{
  static G4int i, j;

  static G4double rmass[297] = {
    1.00000000000000, 1.25992104989487, 1.44224957030741,
    1.58740105196820, 1.70997594667670, 1.81712059283214,
    1.91293118277239, 2.00000000000000, 2.08008382305190,
    2.15443469003188, 2.22398009056932, 2.28942848510666,
    2.35133468772076, 2.41014226417523, 2.46621207433047,
    2.51984209978975, 2.57128159065824, 2.62074139420890,
    2.66840164872194, 2.71441761659491, 2.75892417638112,
    2.80203933065539, 2.84386697985157, 2.88449914061482,
    2.92401773821287, 2.96249606840737, 3.00000000000000,
    3.03658897187566, 3.07231682568585, 3.10723250595386,
    3.14138065239139, 3.17480210393640, 3.20753432999583,
    3.23961180127748, 3.27106631018859, 3.30192724889463,
    3.33222185164595, 3.36197540679896, 3.39121144301417,
    3.41995189335339, 3.44821724038273, 3.47602664488645,
    3.50339806038672, 3.53034833532606, 3.55689330449006,
    3.58304787101595, 3.60882608013869, 3.63424118566428,
    3.65930571002297, 3.68403149864039, 3.70842976926619,
    3.73251115681725, 3.75628575422107, 3.77976314968462,
    3.80295246076139, 3.82586236554478, 3.84850113127680,
    3.87087664062780, 3.89299641587326, 3.91486764116886,
    3.93649718310217, 3.95789160968041, 3.97905720789639,
    4.00000000000000, 4.02072575858906, 4.04124002062219,
    4.06154810044568, 4.08165510191735, 4.10156592970235,
    4.12128529980856, 4.14081774942285, 4.16016764610381,
    4.17933919638123, 4.19833645380841, 4.21716332650875,
    4.23582358425489, 4.25432086511501, 4.27265868169792,
    4.29084042702621, 4.30886938006377, 4.32674871092222,
    4.34448148576861, 4.36207067145484, 4.37951913988789,
    4.39682967215818, 4.41400496244210, 4.43104762169363,
    4.44796018113863, 4.46474509558454, 4.48140474655716,
    4.49794144527541, 4.51435743547400, 4.53065489608349,
    4.54683594377634, 4.56290263538697, 4.57885697021333,
    4.59470089220704, 4.61043629205845, 4.62606500918274,
    4.64158883361278, 4.65700950780383, 4.67232872835526,
    4.68754814765360, 4.70266937544151, 4.71769398031653,
    4.73262349116337, 4.74745939852340, 4.76220315590460,
    4.77685618103502, 4.79141985706278, 4.80589553370533,
    4.82028452835046, 4.83458812711164, 4.84880758583988,
    4.86294413109428, 4.87699896107331, 4.89097324650875,
    4.90486813152402, 4.91868473445873, 4.93242414866094,
    4.94608744324870, 4.95967566384230, 4.97318983326859,
    4.98663095223865, 5.00000000000000, 5.01329793496458,
    5.02652569531348, 5.03968419957949, 5.05277434720856,
    5.06579701910089, 5.07875307813270, 5.09164336965949,
    5.10446872200146, 5.11722994691205, 5.12992784003009,
    5.14256318131647, 5.15513673547577, 5.16764925236362,
    5.18010146738029, 5.19249410185110, 5.20482786339420,
    5.21710344627617, 5.22932153175598, 5.24148278841779,
    5.25358787249290, 5.26563742817144, 5.27763208790408,
    5.28957247269421, 5.30145919238090, 5.31329284591305,
    5.32507402161499, 5.33680329744389, 5.34848124123936,
    5.36010841096536, 5.37168535494483, 5.38321261208728,
    5.39469071210959, 5.40612017575022, 5.41750151497718,
    5.42883523318981, 5.44012182541480, 5.45136177849642,
    5.46255557128140, 5.47370367479843, 5.48480655243262,
    5.49586466009501, 5.50687844638735, 5.51784835276224,
    5.52877481367887, 5.53965825675446, 5.55049910291155,
    5.56129776652123, 5.57205465554262, 5.58277017165842,
    5.59344471040698, 5.60407866131077, 5.61467240800149,
    5.62522632834186, 5.63574079454424, 5.64621617328617,
    5.65665282582291, 5.66705110809706, 5.67741137084543,
    5.68773395970313, 5.69801921530506, 5.70826747338486,
    5.71847906487132, 5.72865431598244, 5.73879354831717,
    5.74889707894483, 5.75896522049240, 5.76899828122963,
    5.77899656515213, 5.78896037206240, 5.79888999764900,
    5.80878573356370, 5.81864786749696, 5.82847668325146,
    5.83827246081400, 5.84803547642573, 5.85776600265065,
    5.86746430844261, 5.87713065921074, 5.88676531688334,
    5.89636853997037, 5.90594058362449, 5.91548169970072,
    5.92499213681474, 5.93447214039994, 5.94392195276313,
    5.95334181313905, 5.96273195774369, 5.97209261982640,
    5.98142402972088, 5.99072641489509, 6.00000000000000,
    6.00924500691737, 6.01846165480645, 6.02765016014974,
    6.03681073679769, 6.04594359601251, 6.05504894651110,
    6.06412699450696, 6.07317794375132, 6.08220199557340,
    6.09119934891978, 6.10017020039306, 6.10911474428961,
    6.11803317263662, 6.12692567522842, 6.13579243966196,
    6.14463365137169, 6.15344949366368, 6.16224014774904,
    6.17100579277672, 6.17974660586564, 6.18846276213620,
    6.19715443474113, 6.20582179489575, 6.21446501190772,
    6.22308425320606, 6.23167968436975, 6.24025146915571,
    6.24879976952624, 6.25732474567597, 6.26582655605827,
    6.27430535741117, 6.28276130478279, 6.29119455155629,
    6.29960524947437, 6.30799354866327, 6.31635959765638,
    6.32470354341737, 6.33302553136292, 6.34132570538500,
    6.34960420787280, 6.35786117973420, 6.36609676041689,
    6.37431108792909, 6.38250429885991, 6.39067652839931,
    6.39882791035777, 6.40695857718556, 6.41506865999165,
    6.42315828856237, 6.43122759137962, 6.43927669563891,
    6.44730572726691, 6.45531481093889, 6.46330407009565,
    6.47127362696036, 6.47922360255497, 6.48715411671635,
    6.49506528811226, 6.50295723425693, 6.51083007152643,
    6.51868391517377, 6.52651887934375, 6.53433507708757,
    6.54213262037718, 6.54991162011937, 6.55767218616971,
    6.56541442734614, 6.57313845144243, 6.58084436524139,
    6.58853227452786, 6.59620228410148, 6.60385449778925,
    6.61148901845794, 6.61910594802623, 6.62670538747667,
    6.63428743686750, 6.64185219534421, 6.64939976115097,
    6.65693023164187, -5.56381e-67, 3.44267e+07
  };

  static G4double Wapstra[11][250] = {
    {0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0,
     3.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 6.0, 6.0, 7.0,
     7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 10.0, 11.0, 11.0, 11.0,
     11.0, 11.0, 12.0, 13.0, 13.0, 14.0, 14.0, 15.0, 15.0, 16.0,
     16.0, 17.0, 17.0, 18.0, 18.0, 18.0, 19.0, 19.0, 20.0, 20.0,
     21.0, 21.0, 22.0, 22.0, 23.0, 23.0, 24.0, 24.0, 25.0, 25.0,
     26.0, 26.0, 27.0, 27.0, 28.0, 28.0, 28.0, 29.0, 29.0, 29.0,
     30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 31.0, 31.0,
     32.0, 32.0, 33.0, 33.0, 34.0, 34.0, 35.0, 35.0, 36.0, 36.0,
     36.0, 36.0, 36.0, 37.0, 37.0, 37.0, 38.0, 38.0, 39.0, 40.0,
     40.0, 41.0, 41.0, 42.0, 42.0, 43.0, 43.0, 44.0, 44.0, 45.0,
     45.0, 46.0, 46.0, 47.0, 47.0, 47.0, 47.0, 48.0, 48.0, 48.0,
     49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0,
     50.0, 50.0, 51.0, 51.0, 52.0, 52.0, 53.0, 53.0, 54.0, 54.0,
     54.0, 54.0, 55.0, 55.0, 55.0, 56.0, 57.0, 58.0, 58.0, 59.0,
     59.0, 60.0, 60.0, 61.0, 61.0, 62.0, 62.0, 63.0, 63.0, 63.0,
     64.0, 64.0, 65.0, 65.0, 66.0, 66.0, 67.0, 67.0, 67.0, 67.0,
     68.0, 68.0, 68.0, 69.0, 69.0, 69.0, 70.0, 70.0, 71.0, 71.0,
     72.0, 72.0, 72.0, 72.0, 73.0, 73.0, 74.0, 74.0, 74.0, 74.0,
     75.0, 76.0, 76.0, 76.0, 76.0, 77.0, 77.0, 77.0, 78.0, 78.0,
     78.0, 79.0, 79.0, 79.0, 80.0, 80.0, 81.0, 81.0, 81.0, 81.0,
     82.0, 82.0, 82.0, 82.0, 83.0, 83.0, 84.0, 84.0, 85.0, 85.0,
     86.0, 86.0, 87.0, 87.0, 87.0, 87.0, 87.0, 88.0, 88.0, 88.0,
     89.0, 89.0, 90.0, 90.0, 90.0, 91.0, 91.0, 91.0, 92.0, 92.0,
     93.0, 93.0, 94.0, 94.0, 94.0, 94.0, 95.0, 95.0, 96.0, 96.0},
    {8.07131, 13.13573546, 14.9497831614761, 25.9197974264349,
     33.7897359197236, 17.5968624734944, 26.1107959337053, 31.6087529649761,
     24.9546049697993, 33.8297356071101, 40.9396800400558, 25.0298043820859,
     34.8997272446982, 40.9696798055956, 29.5297692130642, 37.9997030171499,
     45.2696461996415, 25.3698017248709, 34.4297309179072, 22.1998264994929,
     26.9497893766366, 9.48992583244087, 17.9498597146800, 8.64993239732492,
     12.8398996510580,-0.189998515085750, 6.66994787169448,-1.12999116866788,
     2.65997921120050, 8.37993450746622, 10.6099170792621, 16.4098717503008,
     4.12996772265340,-4.14996756634664,-0.839993435115946,-12.6699009796655,
    -7.00994521447950,-14.5598862086764,-12.2999038713406,-22.2398261868793,
    -18.0998585423793,-24.4198091494421,-23.1398191530750,-32.2707477912223,
    -29.7297676499965,-29.7297676499965,-35.6977210080584,-32.2197481898045,
    -41.2856773359488,-39.5716907314384,-43.2196622210848,-40.1396862923263,
    -46.8896335387937,-45.3296457307212,-49.0096169702768,-46.2096388532237,
    -52.7895874282986,-52.0495932116488,-55.4775664206696,-52.9495861778445,
    -59.0095388168952,-58.9295394421223,-61.8495166213348,-59.7905327131161,
    -65.1239910300100,-66.0204840235594,-63.4695039604870,-65.3894889550377,
    -65.9394846566018,-63.3895045857141,-67.3234738401737,-68.1334675097498,
    -65.0294917685595,-65.6694867667431,-62.4595118539786,-62.5495111505981,
    -58.9095395984290,-58.0795460851597,-62.8095091186102,-59.5295347529194,
    -66.3394815304665,-65.7894858289025,-69.9494533170957,-66.1594829372274,
    -72.5694328409098,-70.8594462051380,-74.2094200237552,-71.0894444076102,
    -76.7893998601827,-75.1794124428772,-71.7694390931803,-69.1494595693663,
    -64.9194926282467,-69.4594571466114,-66.5494798892455,-62.7695094312237,
    -69.0794601164399,-67.3794734025148,-71.4994412033216,-76.5994013450970,
    -73.0494290895475,-76.3594032207781,-75.4094106453494,-80.4993708652782,
    -77.1393971248144,-80.0293745384871,-79.5093786024629,-83.8193449183555,
    -80.8093684425233,-82.9293518740064,-82.5293550001417,-86.3253253331180,
    -83.6393463251163,-85.1593344458023,-84.9093363996369,-82.6193542967613,
    -82.2393572665898,-86.7063223554742,-84.2293417140668,-83.9803436600860,
    -85.8413291157417,-83.5993466377299,-83.4393478881840,-81.0993661760753,
    -80.4993708652782,-77.8993911851574,-77.1693968903542,-74.3394190077612,
    -73.1194285424738,-70.0794523011018,-77.4793944675994,-76.3894029863180,
    -78.9793827445922,-73.8694226809702,-77.5993935297588,-74.8294151782455,
    -76.7194004072564,-71.7294394057938,-75.7494079881344,-73.1794280735535,
    -68.9994607416670,-66.0494837969146,-68.3594657434834,-63.9295003654315,
    -61.7195176373288,-65.5594876264303,-67.5394721520607,-70.7094473774387,
    -67.4694726991344,-68.6794632425752,-67.4394729335945,-70.1454517852895,
    -67.3594735588216,-68.4494650401030,-67.0994755908095,-69.3674578656226,
    -66.8594774664907,-67.2394744966622,-65.9294847347552,-63.5395034134133,
    -65.5064880406432,-64.3594970048361,-64.6794945039279,-62.1095145893469,
    -63.6105028585243,-62.5825108926920,-62.3155129793873,-60.2695289695691,
    -58.7925405128236,-56.0995615595293,-57.7135489455735,-56.4905585037320,
    -53.7295800818807,-53.8495791440401,-52.2895913359676,-49.5896124373807,
    -50.9856015271686,-49.6596118903070,-49.1096161887430,-46.6796351800147,
    -47.4026295295252,-45.9896405725980,-43.2686618381332,-41.4796758197731,
    -41.3596767576137,-38.5996983279470,-39.8926882227148,-38.6566978824728,
    -35.4697227899555,-34.2197325591282,-34.3427315978416,-35.8747196247435,
    -33.3867390693049,-32.4167466501829,-29.6897679626100,-29.4397699164446,
    -28.4297778099361,-25.5198005525702,-27.4197857034277,-26.5997921120050,
    -23.7398144638721,-23.8598135260315,-22.9798204035291,-20.1998421301692,
    -22.2988257257744,-20.9548362295889,-21.0408355574698,-16.7678689524097,
    -13.6498933206341,-9.25092770030669,-10.4918180022536,-7.56194090041284,
    -3.13997545983818,-0.185298551817839, 1.70998663577175, 5.96995334243119,
     5.95995342058457, 8.35453470597581, 10.5299177044892, 14.1998890221981,
     14.3798876154373, 16.3698720629143, 18.3821563361094, 21.7098303290086,
     23.7898140731052, 27.4597853908141, 29.5797688222972, 28.9407738162983,
     32.7197442821354, 34.5597299019132, 35.9097193512067, 39.1496940295111,
     38.7319972939778, 40.6116826034867, 44.1496549528203, 45.5396440895002,
     47.6396276772901, 51.2695993076126, 50.5718047611555, 52.7115880378949,
     54.3095755489846, 57.5195504617491, 57.7520486446830, 59.8025326193320,
     63.1565064066879, 65.2894897365716, 67.1294753563494, 70.4894490968131,
     70.7474470804559, 72.9854295897291},
    {7.28897545999993, 0.0, 14.9311833068413, 2.42492104827388,
     11.3899109832984, 14.0871899029867, 14.9080834873757, 20.9467362928931,
     11.3479113115426, 12.6075014673426, 20.1758423177373, 13.3693955128365,
     16.5618705623694, 23.6568151125452, 9.87312283760328, 13.6928929845746,
     21.0598354089784, 13.2738962592013, 15.5998780807247, 3.79897030953033,
     8.11993653945414, 2.82597791385436, 3.34997381861717,-5.94895350655329,
    -2.14998319702296,-6.88794616795076,-5.62995599964616,-15.0162826417561,
    -10.7499159851148,-9.78992348783942,-3.89996952018118,-2.88997741367272,
    -9.36992677028145,-20.2498417394023,-15.0398824573141,-20.7698376754264,
    -19.0098514304216,-26.8617900643864,-22.9998202472223,-27.5397847655871,
    -27.3997858597344,-34.4197309960606,-31.9797500654857,-35.8067201561865,
    -36.6107138726546,-35.4197231807224,-42.3425690759179,-44.2156544370079,
    -46.5546361569320,-44.5386519126537,-49.7326113197873,-49.4686133830366,
    -51.8625946731170,-49.9296097801657,-55.1058693256308,-55.2645680853366,
    -57.4865507196552,-56.2095606998421,-60.6609259106458,-61.4365198490695,
    -62.8965084386758,-61.5035193254418,-65.5120879968773,-67.0973756072217,
    -67.2609743286324,-66.2561821814842,-67.3044739886652,-70.0057528770922,
    -68.4164652980092,-69.5593563658592,-70.1409518204585,-68.5904639381403,
    -69.7294550364701,-68.0194684006984,-68.5594641804158,-66.4394807489327,
    -66.4094809833928,-63.6795023192660,-69.5694562869243,-69.4294573810716,
    -72.6394322938361,-70.1894514414146,-75.4094106453494,-75.9414064875895,
    -78.6693851673470,-75.9594063469134,-80.7063692475032,-79.6883772035174,
    -81.7163613540116,-79.5693781335427,-77.9693906380837,-75.1194129117975,
    -72.9194301055414,-78.9593829008989,-75.1394127554907,-73.0694289332407,
    -76.2794038460052,-73.1894279954001,-77.8893912633108,-79.9593750855608,
    -78.9493829790523,-83.5613469347127,-80.6093700055910,-82.6993536715342,
    -82.5393549219883,-86.3323252784107,-83.7093457780427,-85.0193355399497,
    -85.1093348365692,-88.3343096321037,-86.0293276464581,-86.6193230354086,
    -87.0393197529666,-90.0188964663850,-88.0923115234155,-88.7169066419553,
    -86.4153246297376,-87.4493165486780,-87.7293143603833,-85.6993302255197,
    -89.2011028577686,-89.9452970415939,-87.8203136491875,-88.2393103745608,
    -85.9013286468214,-86.0233276933502,-83.5993466377299,-83.4393478881840,
    -80.6393697711308,-80.3793718031187,-82.0993583607371,-79.6093778209291,
    -82.9293518740064,-82.6693539059944,-83.7953451059236,-79.4293792276900,
    -82.2143574619732,-80.0293745384871,-80.6293698492842,-77.2393963432806,
    -74.9994138496380,-70.9494455017576,-74.0094215868228,-72.0294370611924,
    -67.8194699637660,-69.4594571466114,-72.2394354199714,-72.5094333098301,
    -71.3694422193156,-73.6814241502537,-70.9444455408343,-71.2894428445426,
    -70.7594469866718,-72.4534337474890,-70.1954513945226,-70.0824522776558,
    -69.4644571075348,-70.6904475259302,-68.5614641647851,-67.9424690024794,
    -67.4654727303957,-65.7594860633626,-66.3814812022223,-65.9664844455876,
    -64.8954928158148,-63.0665071100683,-63.2855053985092,-62.9845077509260,
    -60.9165239130453,-60.1035302669153,-59.2045372929043,-57.3795515558964,
    -56.2255605747967,-56.9395549946452,-54.6905725713407,-53.4895819575619,
    -52.3815906169565,-50.2996068884906,-50.4616056224058,-49.7786109602818,
    -48.4246215422496,-46.4166372354487,-45.2786461293035,-42.8206653394047,
    -43.3696610487840,-42.4976678637589,-41.2046779689911,-39.0056951549197,
    -37.9697032516101,-35.5197223991886,-36.3877156154751,-34.8257278230333,
    -34.5187302223421,-32.5137458920951,-31.6917523163031,-32.6517448135784,
    -30.4307621714445,-29.9207661572669,-29.1037725423982,-27.2997866412682,
    -26.3997936750726,-27.3557862036093,-25.2768024516973,-24.7028069377014,
    -23.8368137057843,-22.2688259602345,-22.4628244440589,-21.7588299460570,
    -17.6238622624803,-14.7378848175462,-11.8649072710127,-8.13493642222407,
    -5.24295902418203,-1.20899055125617,-0.540495775809725, 1.76898617466680,
     4.38196575318818, 8.09893670357625, 8.83063098509331, 10.5989171652308,
     13.2648963295393, 16.3378723130052, 17.2346653042099, 18.8128529700432,
     21.9871281618153, 23.6655150445517, 27.1847875400321, 28.8947741758039,
     30.7197599128117, 33.7597361541837, 33.8119357462231, 35.4469229681452,
     37.4868070256369, 40.3486846589206, 42.3196692548891, 42.4416683014178,
     45.3883452719609, 47.3066302797977, 49.3060146538106, 52.2095919611947,
     52.9525861543984, 54.7145723837726, 57.1696531963359, 59.8781320284925,
     61.8968162516693, 64.9194926282467, 65.5294878608904, 67.3884733321768,
     69.8474541142602, 72.9494298710813},
    {0.0, 0.0, 0.0, 25.1298036005521,
     11.6799087168503, 18.3748563931613, 15.7699767513357, 4.94172137847450,
     12.4160029639799, 12.0516058118891, 8.66783225743037, 0.0,
     3.12500557683379, 3.01989639830399, 0.101509206665023, 5.68155559637471,
     7.86993849328868,-0.783023880355761, 3.33137396398246,-0.01709986635771748,
    -0.04699963267910651,-8.02603727331440,-5.15505971115025,-8.41743421439104,
    -9.35692687188084,-16.2122732946116,-14.5848860132929,-16.8480683256196,
    -18.2118576670614,-15.8898758142766,-15.0998819883938,-11.2899117648322,
    -20.5698392384941,-24.5498081334482,-24.9398050854663,-30.6656603356215,
    -26.9077897048808,-29.7977671185535,-29.8027670794768,-35.0399261489878,
    -33.0674415647424,-35.0225262849747,-36.5877140524074,-41.4656759291879,
    -40.8092810591759,-43.1378628603794,-44.3301535421517,-44.4976522330826,
    -48.5583204973389,-51.4316980407462,-52.1986920463819,-51.4384979876019,
    -55.2832679391898,-56.9308550626386,-57.7095489768348,-56.9083552384838,
    -60.1785296807649,-62.1513142626658,-62.2259136796416,-61.6461182109746,
    -64.2185981060172,-66.7448783621284,-65.5779874818465,-65.4224886971316,
    -65.9090848941880,-68.8977615364869,-67.8790694979719,-67.0846757064765,
    -69.3209582290358,-68.9046614825611,-69.9052536625337,-72.5820327424365,
    -71.2929428171889,-73.4215261814601,-71.8555384202797,-73.2129278117397,
    -71.2137434361637,-71.7594391713337,-73.7194238532709,-72.0594368267322,
    -76.3904029785027,-77.5853936391736,-79.0243823929020,-77.7583922871201,
    -81.4711632703325,-83.2623492714988,-84.5950388559977,-82.6013544374374,
    -86.2023262944046,-85.9340283912599,-83.6653461219176,-82.8913521709893,
    -80.2793725846526,-82.3813561568117,-81.2323651366353,-78.4293870430282,
    -82.9535516848753,-81.2913646755303,-82.3453564381639,-86.1883264038194,
    -83.5153472942183,-84.5993388223917,-84.9093363996369,-88.0983114765235,
    -85.9373283654693,-86.3713249736125,-86.8593211597275,-89.5223003474819,
    -87.6053153294852,-87.4553165017859,-88.2253104839755,-90.5771921030817,
    -89.0496040417923,-88.5753077486072,-89.5403002068059,-88.2523102729614,
    -88.9433048725627,-91.6528836961225,-90.0659960982826,-91.1010880086261,
    -89.5876998363588,-88.3226097235431,-89.2168027350678,-87.6127152716517,
    -88.2513102807767,-86.4013247391524,-86.7033223789202,-84.7293378063978,
    -84.6293385879316,-82.3793561724424,-85.2003341253735,-85.2123340315894,
    -85.9013286468214,-83.9693437460547,-86.5053239263572,-86.4243245593996,
    -86.5593235043289,-82.7693531244606,-84.9243362824068,-83.2843490995614,
    -79.9793749292540,-77.8193918103844,-78.3093879808687,-74.9294143967117,
    -72.9194301055414,-75.7594079099811,-75.4394104108893,-77.4063950381191,
    -74.3734187420397,-73.5494251818784,-73.3854264635938,-74.7604157175039,
    -72.5564329425092,-71.7254394370552,-71.8244386633367,-72.5354331066313,
    -70.8244464786749,-69.4744570293814,-69.5354565526458,-67.8394698074593,
    -68.0554681193462,-68.1804671424290,-66.3784812256683,-64.9364924953859,
    -64.5174957700126,-64.9204926204314,-62.5365112521975,-61.3055208728788,
    -61.2685211620463,-59.7905327131161,-59.3015365348165,-59.2495369412141,
    -57.5455502585503,-55.5615657641812,-55.1585689137625,-53.3805828094337,
    -52.8785867327335,-52.4335902105589,-50.3466065211697,-48.9136177205493,
    -48.2366230115332,-48.2276230818713,-46.3466377825223,-45.6866429406455,
    -43.8016576725579,-41.9096724591777,-41.2076779455451,-41.1246785942182,
    -38.9776953737492,-38.6986975542286,-36.6977131927202,-36.2827164360856,
    -34.4577306990777,-34.7647282997689,-32.8017436412777,-31.1617564584323,
    -31.1497565522163,-29.5907687363285,-29.5567690020500,-29.5137693381096,
    -27.6717837339625,-25.9877968949919,-25.7687986065510,-24.3528096730698,
    -23.7768141747046,-23.7948140340285,-20.0578432399472,-18.8788524542309,
    -18.2678572294025,-14.8008843251799,-12.4439027459319,-10.3809188689746,
    -6.66294792640184,-4.47896499510038,-1.26199013704324, 2.23698251708854,
     3.64897148183106, 5.21195926645751, 8.61693265523108, 11.4699103580713,
     12.9568987366635, 14.3118881468803, 17.8248606915973, 20.2188419816778,
     21.6258309854970, 24.3008100794674, 25.8497979735086, 26.7577908771815,
     29.5806688152634, 30.8610588085045, 33.4228387871712, 35.9337191636386,
     36.9144114991365, 38.1423019026827, 40.9160802244977, 43.4256606111251,
     44.8689493312475, 47.4522291418845, 48.5847202910140, 50.1224082733685,
     52.9317863169575, 55.4622665402443, 57.1769531392839, 58.4491431966107,
     61.0008232542123, 62.6155106347858, 65.4834882203960, 67.9894686351585,
     69.7210551021190, 71.1692437839463},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 27.9397816394518, 22.9217208576002,
     28.9118740421616, 15.7027772765264, 10.6499167666486, 17.3378644976670,
     5.34555822232834, 2.86341762124810, 2.85537768408342,-4.73698297858683,
    -0.809893670357625, 0.872493181117456,-1.48736837562233,-7.04294495657334,
    -5.73305519388480,-5.18395948528698,-9.52952552295348,-13.9304911276502,
    -13.1906969094374,-12.2075045934779,-17.1941656207311,-21.4910320390045,
    -21.8935288933309,-24.4315090580027,-22.9485206481492,-24.0918117128731,
    -26.3366941682204,-29.9309660775505,-28.8460745564109,-29.5215692771499,
    -31.7615517707925,-34.7147286905358,-33.2407402103442,-33.5349379110718,
    -35.5594220889197,-38.5435987663875,-38.4050998488119,-37.8104044965934,
    -41.0661790514155,-41.7552736658659,-44.9306488490412,-48.4873210522279,
    -47.9565252006094,-49.2189153345265,-51.4473979180454,-55.4148669106913,
    -54.6869725994759,-55.5538658243593,-57.4781507853041,-60.6036263584647,
    -59.3419362190768,-59.8435322989032,-61.1524220694071,-64.4696961435858,
    -61.9802155998701,-62.7955092280249,-62.2106137992162,-66.0006841783031,
    -62.6533103393660,-63.7228019808619,-66.8779773219069,-66.9714765911728,
    -67.0954756220709,-70.5608485387980,-67.8924693932463,-68.2314667438467,
    -70.9484455095729,-70.8591462074826,-73.0333292153744,-72.2900350245152,
    -73.9151223238092,-72.7394315123023,-75.9200066548377,-77.7606922691448,
    -77.9753905911917,-77.4973943269233,-79.9839748933035,-82.4312557668264,
    -82.1581579011952,-82.7370533768960,-84.8682367208473,-87.9099129489332,
    -87.6946146315755,-86.4800241240852,-86.3488251494576,-84.8213370873867,
    -84.2263417375129,-87.2632180031124,-85.6627305115611,-85.4440322207756,
    -85.6109309163956,-83.5293471848035,-85.9688281192861,-86.0181277339899,
    -86.3263253253027,-89.0998036494623,-87.2607180226507,-86.9513204407164,
    -87.8543133834660,-89.9122972995001,-88.3703093507515,-87.6013153607466,
    -88.7213066075678,-90.3482938920126,-89.2533024498079,-87.9993122502420,
    -89.3713015275980,-90.5592922429763,-90.0343963452473,-91.5253846925781,
    -90.3981935020273,-87.9663125081481,-89.4823006600955,-88.4203089599846,
    -88.5073082800502,-90.3032942437029,-89.1648031414653,-90.5175925688759,
    -89.0183042864124,-90.0652961037533,-88.2843100228706,-88.9916044950819,
    -87.0063200108728,-87.3473173458425,-87.4503165408626,-85.7053301786277,
    -87.6613148918263,-88.1243112733247,-87.6643148683803,-86.3573250830272,
    -87.7323143369373,-88.2723101166547,-87.2303182602370,-84.3193410106864,
    -83.0073512644101,-80.0173746322712,-81.6093621902528,-80.4303714045365,
    -77.1193972811212,-76.8393994694158,-78.1433892782149,-76.8693992349557,
    -76.0624055419336,-77.0483978360102,-74.5734171789721,-72.8834303868936,
    -73.3624266433466,-73.7034239783163,-72.0704367407635,-70.0974521604257,
    -70.7664469319645,-70.4094497220402,-69.1704594052442,-69.6734554741291,
    -67.2024747858297,-66.0464838203606,-65.1674906900428,-65.9394846566018,
    -62.9235082276617,-61.8735164337667,-60.5825265233683,-61.5645188487062,
    -60.3605282583734,-60.7585251478688,-57.8205481093323,-56.7255566671276,
    -56.8705555339035,-55.8295636696706,-54.5475736889341,-54.5665735404427,
    -51.7205957828951,-50.5196051691162,-49.2826148366895,-49.6236121716592,
    -46.4396370556959,-45.4296449491874,-45.7906421278503,-44.1906546323914,
    -42.7866656051262,-42.9866640420585,-39.7096896529217,-38.3227004927957,
    -38.4796992657876,-36.6997131770895,-35.6977210080584,-32.7677439069992,
    -33.3597392803190,-32.2557479084523,-32.5717454388055,-31.8457511127410,
    -30.7347597955817,-30.9637580058692,-28.0797805453045,-27.0597885169494,
    -27.1847875400321,-25.9417972544975,-24.7938062265057,-25.1168037021515,
    -21.0698353308250,-20.0328434353307,-17.1498659669506,-17.4748634269657,
    -16.3728720394683,-15.9628752437570,-11.6529089278644,-8.62493259270837,
    -6.58894850473687,-3.38897351381898,-1.17899078571631, 0.244998085242151,
     4.30696633933855, 7.04994490186598, 9.37692671557408, 10.2629197911845,
     14.5178865369206, 16.6168701325258, 19.2558495078484, 19.9928437479442,
     22.3028256945130, 23.1888187701234, 25.8060983150389, 28.8697743711873,
     29.8867664229884, 32.1652486157404, 33.7797359978770, 34.5967296127457,
     38.0097029389966, 39.9506877694252, 41.0391792624296, 42.8886648079617,
     45.0866476298484, 46.1604392377383, 49.3886140082636, 51.4425979555591,
     53.6955803476022, 54.8010717077459, 58.6845413568801, 60.6455260310020,
     61.8105169261330, 64.0194996620510, 66.1494830153808, 67.2424744732161,
     71.1154442044115, 73.1694281517069},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 35.0847257988607,
     0.0, 39.4996912941427, 25.2298028190182, 32.0697493621052,
     23.1048194266118, 8.00823741242742, 17.6598619811281, 10.6919164384044,
     1.95164474711713, 5.31895843021633, 1.75088631612442, 6.84394651182564,
    -2.18578291723385,-0.394096919975231,-5.47055724541107,-0.05199959360241571,
    -8.91283034267252,-7.14304417425799,-12.3852032046923,-7.15944404608645,
    -16.9491675354889,-20.2043420950002,-24.4393089970431,-24.3045100505506,
    -26.5856922222012,-24.4381090064215,-29.0134732481232,-30.2310637321675,
    -30.9476581316962,-28.8017749026303,-33.8059357931151,-34.8465276604742,
    -35.1382253807401,-32.1204489658676,-36.1847172019887,-37.5459065637504,
    -39.0036951705504,-44.1223551661790,-42.0007717472004,-44.4724524300291,
    -45.3286457385366,-50.2576072167348,-48.2394229896503,-50.7038037295309,
    -50.9438018538497,-56.2509603762871,-54.0234777849528,-56.0362620542402,
    -56.0985615673446,-60.2238293267301,-56.3517595885010,-58.3428440273812,
    -56.5795578081669,-61.1685219435801,-56.6895569484797,-58.8355401767641,
    -56.4095591367744,-61.6205184110473,-62.4495119321319,-58.7695406925764,
    -63.1195066958554,-64.3384971689582,-63.0895069303155,-67.8934693854310,
    -68.2084669235995,-72.2121356333301,-72.1684359748604,-75.2586118239024,
    -74.6055169280997,-77.0308979727786,-76.0694054872262,-75.8904068861717,
    -77.7063926935176,-80.5903701540824,-78.9863826898848,-79.7513767111511,
    -81.0943662151520,-84.5114395093599,-83.0065512706623,-84.2973411826238,
    -84.8588367943115,-88.7639062746344,-87.8918130903908,-88.4554086856662,
    -87.1160191535302,-86.3664250119077,-86.7858217341548,-85.6073309445309,
    -87.5438158101285,-88.1147113483519,-87.3255175162168,-89.2209027030249,
    -87.9509126285043,-86.7763218084005,-88.0233120626739,-89.3993013087685,
    -88.4213089521693,-86.9283206204691,-88.4033090928453,-89.2503024732539,
    -88.5393080299594,-86.4083246844450,-88.4043090850300,-88.6573071077495,
    -88.3313096555497,-84.8693367122504,-87.0043200265034,-86.9293206126538,
    -88.6533071390108,-87.6703148214882,-87.1883185884812,-89.4033012775072,
    -86.1393267867710,-86.1593266304642,-87.9693124847021,-87.3603172442431,
    -88.8403056775426,-87.9103129458071,-88.9793045912106,-87.7333143291219,
    -88.5043083034962,-86.8963208705600,-88.4203089599846,-89.2853021997171,
    -88.0883115546769,-86.9083207767759,-87.8693132662359,-88.9053051695456,
    -87.1293190495862,-86.5233237856811,-86.9653203313016,-88.0803116171996,
    -85.4373322731383,-84.5343393303887,-83.0643508189358,-80.7493689114436,
    -79.6243777036991,-80.9223675593901,-79.0393822756719,-79.3343799701471,
    -77.1343971638911,-74.7554157565805,-74.6494165850064,-74.7024161707935,
    -73.1184285502891,-70.2394510506477,-71.2554431102641,-70.5264488076456,
    -69.4244574201483,-66.4324808036401,-67.3174738870658,-66.3874811553303,
    -65.1964904633980,-66.3344815695432,-62.7695094312237,-61.9775156209716,
    -60.1605298214410,-61.5815187158455,-57.4495510088227,-57.0995537441911,
    -57.8805476404120,-57.3185520326320,-55.2995678117998,-56.3295597620015,
    -55.2695680462599,-51.9795937587225,-52.3495908670474,-51.4695977445449,
    -49.7196114213867,-50.4296058724966,-46.5896358833952,-45.8286418308675,
    -43.4096607361705,-44.5796515922248,-43.4896601109434,-44.2326543041472,
    -40.2896851200256,-39.1556939826190,-36.8097123174024,-37.7877046740016,
    -36.5697141930835,-37.3177083472106,-33.8697352944965,-31.9697501436390,
    -31.0197575682103,-32.2057482992192,-31.0497573337502,-27.3497862505013,
    -28.3297785914699,-27.4997850782006,-25.2798024282513,-26.1597955507538,
    -25.3268020609304,-21.0398355652851,-21.5998311886958,-20.8198372846595,
    -17.5758626376165,-18.1898578389989,-13.3098959778491,-12.6399012141257,
    -12.8878992759218,-11.9759064035102,-8.76093152982238,-8.66593227227951,
    -5.70595540568046,-4.32796617521645, 0.308997585060509, 2.97497674936898,
     5.88095403799628, 6.64394807489327, 11.5599096546909, 13.7468925625463,
     16.9338676550636, 17.1968655996297, 22.3298254834989, 23.7978140105825,
     24.3198099309760, 26.0287965745630, 26.8317902988465, 29.2207716280036,
     31.2007561536341, 31.6067529806068, 35.6257215707627, 37.2897085660400,
     40.0416870582294, 40.3416847136280, 42.1596705053432, 46.0196403381379,
     46.6396354926282, 48.4166216047723, 51.0896007143734, 51.7115958532331,
     56.0995615595293, 57.7995482734544, 60.9095239677527, 61.4645196302400,
     63.3765046873135, 64.0956990665223, 68.5494642585692, 70.2194512069544,
     73.4994255726453, 74.0684211257179},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 33.6097373264845, 0.0, 24.1098115721970,
     16.4778712188578, 25.3198021156378, 12.9298989476776, 17.5678627001392,
     10.9119147190300, 18.2098576826921, 6.76794710579134, 10.7399160632682,
     3.82397011414688, 11.2599119992923,-0.589995388950486, 4.18996725373311,
    -3.15997530353142,-14.0618901007148,-19.0439511639186,-26.0148966831962,
    -21.0028358544526,-18.3790563603369,-23.0487198650523,-17.4258638099172,
    -24.7992061843029,-22.0598275936402,-27.2817867819443,-20.5268395745536,
    -28.6432761413614,-25.1218036630748,-29.3237708230238,-23.8498136041849,
    -31.8787508548348,-37.0705102791622,-34.6177294486236,-42.8176653628507,
    -37.6127060416858,-42.6253668657402,-40.2276856045765,-48.3316222690761,
    -42.6396667539809,-48.0092247887411,-45.3336456994599,-53.9015787376425,
    -47.6196278335969,-51.6612962463446,-47.5896280680570,-54.1835765337172,
    -47.7496268176029,-51.7695953999435,-47.5496283806705,-54.4295746111440,
    -47.1496315068058,-51.5195973537780,-56.6495572610933,-54.1695766431319,
    -56.2995599964616,-61.5895186533228,-56.4895585115474,-58.9295394421223,
    -63.6695023974194,-65.2944896974949,-69.1584594990282,-70.3024505582814,
    -73.2409275929102,-73.4574259008895,-74.4384182340427,-77.8963912086034,
    -75.4444103718126,-76.2124043696328,-76.7364002743956,-80.6403697633155,
    -77.8343916931544,-79.2383807204196,-79.4293792276900,-83.6203464736078,
    -80.6203699196222,-82.6533540310398,-86.6362229033294,-86.4474243788653,
    -87.2083184321745,-88.4116090279780,-87.7114145002778,-88.7942060378297,
    -87.2233183149444,-86.4333244890615,-87.6191152216335,-85.5913310695763,
    -87.4093168612915,-87.9243128363923,-87.4773163298485,-85.1493345239557,
    -87.0743194794298,-87.1303190417708,-86.9863201671795,-84.0993427300608,
    -86.5233237856811,-85.8333291782644,-85.9403283420232,-81.7393611742588,
    -84.4423400493998,-82.1893576573567,-82.4193558598289,-85.3693328045813,
    -85.1633344145410,-81.3693640659340,-83.8193449183555,-83.9993435115946,
    -82.3493564069026,-85.1593344458023,-85.2893334298084,-87.4493165486780,
    -87.1093192058929,-89.1613031688190,-88.3153097805951,-89.8604977043346,
    -88.6968067990436,-89.8803975488094,-88.0653117344297,-87.1743186978959,
    -87.5683156186527,-88.9673046849946,-86.6693226446417,-86.0393275683048,
    -85.9093285842987,-87.5643156499141,-84.8533368372958,-84.6923380955653,
    -86.0173277402422,-83.7893451528156,-83.9993435115946,-83.7453454966905,
    -81.4293635970137,-79.4413791339059,-79.2643805172208,-76.2344041976954,
    -76.4384026033664,-75.7644078709043,-74.1674203519994,-70.8524462598454,
    -71.3284425397444,-70.3914498627163,-69.1564595146589,-65.4294886424242,
    -66.8894772320305,-65.0294917685595,-64.3894967703759,-66.0514837812839,
    -61.6795179499423,-61.5395190440897,-59.1695375664411,-60.8795242022129,
    -56.1595610906090,-56.0995615595293,-53.1495846147768,-55.0995693748674,
    -54.5295738296102,-56.1195614032225,-51.5995967285510,-51.4095982134652,
    -52.3695907107406,-50.0796086078650,-49.4496135315280,-50.5696047783493,
    -46.1196395566041,-45.7696422919725,-42.8896648001463,-44.2196544057466,
    -39.3396925445968,-38.9796953581185,-40.0896866830932,-39.5096912159893,
    -36.4897148183106,-37.8297043457574,-32.8697431098347,-32.4897460796632,
    -33.4097388895521,-32.8757430629427,-30.4797617884929,-25.5898000054965,
    -27.0197888295629,-26.8097904707839,-27.8497823428323,-25.1498034442453,
    -24.6298075082211,-25.8997975827417,-20.6098389258805,-20.4598400981813,
    -21.4098326736100,-17.7798610432875,-17.3598643257296,-18.2498573700786,
    -12.9598987132175,-12.7299005107452,-8.68993208471139,-9.55992528536720,
    -8.99392970884859,-9.60792491023096,-4.21996701927297,-3.68997116140219,
    -3.55597220865751,-0.964992458198676, 2.53098021937912, 3.28497432661415,
     8.70093199874267, 10.8369153051804, 14.4698869120568, 14.6628854036966,
     0.0, 21.9588283829894, 0.0, 0.0,
     0.0, 27.1857875322168, 28.8797742930340, 0.0,
     33.7577361698144, 35.2317246500060, 0.0, 38.3617001879975,
     0.0, 44.4596525300654, 44.6496510451512, 47.8896257234555,
     49.1696157198227, 49.3976139379256, 54.2795757834447, 55.7095646075111,
     59.1895374101343, 59.3315363003563, 64.7994935660873, 65.9694844221416,
     66.3794812178530, 67.9294691040788, 71.5394408907081, 71.8904381475244,
     77.2593961869738, 78.5993857144207},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 4.83996217376331,-7.06994474555921,-13.3288958293577,
    -9.38492665305137,-1.47998843329952,-11.1689127104881,-6.64994802800124,
    -13.1638971188885,-4.45996514359181,-14.0798899600387,-9.03992934934304,
    -15.7798766739638,-8.01993732098796,-18.0198591676064,-13.4998944929348,
    -19.4598479135194,-29.4607697523225,-22.9198208724494,-29.1697720265859,
    -24.4698087586752,-34.4297309179072,-27.2297871883419,-34.2297324809748,
    -29.4097701509047,-39.2096935605908,-31.5297535823878,-38.4996991094809,
    -32.6297449855159,-42.2596697238094, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0,-47.0796320538795, 0.0,
    -46.6196356489350,-51.1396003236065,-46.4996365867756,-53.8695789877333,
    -56.9795546820317,-62.0195152927274,-64.1594985679037,-69.0994599601332,
    -70.2354510819090,-74.1494204926755,-70.8594462051380,-72.1894358107383,
    -71.4594415159351,-75.9984060421152,-72.4394338569037,-73.6914240721003,
    -73.1294284643204,-77.9393908725438,-74.4294183043808,-76.4194027518579,
    -75.2194121302636,-80.1663734677858,-82.1983575870186,-86.8063215739404,
    -86.8023216052018,-84.1553422924019,-86.0123277793189,-85.8203292798638,
    -86.0693273338446,-88.2253104839755,-85.5163316557266,-85.2293338987287,
    -85.4273323512917,-82.3293565632093,-84.7993372593241,-83.8493446838953,
    -84.3353408856410,-80.5853701931591,-83.4993474192637,-81.8993599238047,
    -82.6193542967613,-77.4293948583663,-80.8393682080632,-77.5493939205257,
    -78.5393861833410,-73.0694289332407,-76.6194011887902,-77.6093934516054,
    -80.8493681299098,-78.0693898565499,-78.8293839168929,-82.0493587515040,
    -77.1493970466610,-78.0093903254702,-80.8893678172963,-81.5293628154799,
    -84.0393431989811,-84.3293409325330,-86.2053262709586,-85.9343283889153,
    -87.5623156655448,-86.8623211362815,-86.7253222069828,-88.4523087098938,
    -85.5693312415137,-85.2673336017458,-84.5493392131586,-86.4993239732492,
    -83.2093496857117,-83.1273503265695,-82.0493587515040,-84.2193417922202,
    -84.2023419250810,-85.9483282795005,-82.9583516473616,-81.4153637064284,
    -81.2693648474678,-80.9833670826545,-77.5343940377558,-76.2674039397893,
    -75.1304128258287,-71.0974443450876,-71.6074403592651,-70.1154520197496,
    -69.1544595302896,-64.6344948556181,-66.0544837578379,-63.9295003654315,
    -63.0895069303155,-58.4295433497913,-60.1895295947962,-60.1295300637165,
    -57.3995513995896,-59.3395362378336,-54.3695750800643,-54.5795734388433,
    -51.2595993857659,-53.4795820357152,-47.9496252545353,-48.3996217376331,
    -50.0296089986319,-50.1196082952514,-46.8996334606403,-48.8096185333444,
    -48.4696211905594,-43.5796594075630,-45.1496471374821,-44.9696485442430,
    -41.6196747256258,-43.3496612050908,-37.8897038768371,-37.9297035642236,
    -34.0597338095823,-35.9797188041330,-35.6297215395014,-37.2097091912671,
    -31.7297520193202,-31.6897523319337,-28.0597807016112,-29.8797664776958,
    -29.2097717139724,-30.9597580371306,-25.6697993802694,-22.2898257961124,
    -22.0698275154868,-23.8098139167984,-23.5498159487863,-17.7598611995943,
    -19.4098483042863,-19.2998491639735,-15.0498823791607,-16.7398691712392,
    -16.4098717503008,-10.5199177826426,-11.9699064504022,-11.9699064504022,
    -7.59994060342999,-8.96992989641671,-2.64997928935388,-2.76997835151330,
    -3.75997061432852,-3.63997155216910, 0.779993904036236,-0.109999140312802,
     0.289997733551934, 0.08999929661956566, 5.94995349873795, 7.97993763360149,
     12.1409051139794, 12.3619033867897, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 54.2795757834447, 0.0, 58.0295464759266,
     0.0, 0.0, 0.0, 0.0,
     70.0194527700221, 70.1304519025195, 0.0, 76.9993982189617,
     0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0,-2.20998272810267,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0,-12.4699025427332, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0,-57.5095505399024,-60.6095263123542,
    -65.1094911433324,-67.0894756689629,-65.4594884079641,-70.3894498783470,
     0.0,-67.9094692603856, 0.0,-71.4394416722419,
     0.0,-69.3394580844520, 0.0,-72.9194301055414,
     0.0,-70.9694453454508,-75.9794061906067,-78.9353830884670,
    -83.6093465595765,-82.5703546797128,-83.4513477943999,-86.0743272947679,
    -82.5593547656816,-83.1673500139560,-82.1113582669530,-77.9293909506972,
    -81.3293643785475,-79.4293792276900,-80.5993700837443,-75.8494072066006,
    -79.3393799310704,-76.9893982971151,-78.3993872774883,-72.3994341695173,
    -76.1194050964593,-71.7594391713337,-73.4694258071054,-67.4394729335945,
    -71.4394416722419,-67.0894756689629,-68.6994630862685,-73.2694273701730,
    -74.4794179136139,-68.6694633207286,-72.5294331535233,-73.6394244784979,
    -70.5694484715861,-74.2594196329883,-75.3894108016562,-78.7493845421199,
    -79.4593789932298,-82.5593547656816,-82.7593532026139,-85.4813319292635,
    -85.1153347896772,-87.3023176975327,-83.7693453091224,-83.7393455435825,
    -82.1693578136634,-84.7693374937842,-80.9893670357625,-81.3993638314738,
    -79.4093793839968,-82.0293589078108,-77.4993943112926,-78.1793889968627,
    -80.4693710997383,-81.0593664886888,-79.5103785946476,-81.9633594236231,
    -80.6553696460854,-77.1103973514592,-75.2064122318630,-70.6394479245124,
    -71.4334417191339,-69.1394596475197,-68.6004638599869,-63.7095020848059,
    -64.9534923625252,-62.4395120102853,-62.0565150035598,-56.9395549946452,
    -58.4895428808711,-55.5295660142720,-55.2895678899532,-57.5495502272889,
    -52.0795929771887,-52.3395909452007,-48.7696188459580,-51.2795992294592,
    -45.3596454962611,-46.0996397129109,-41.9496721465642,-44.4996522174519,
    -44.8896491694700,-46.9196333043336,-41.0996787896016,-41.5096755853130,
    -43.3696610487840,-39.6196903563021,-39.7096896529217,-41.8096732407116,
    -35.8197200545871,-36.2697165376850,-32.0097498310255,-34.1197333406620,
    -27.6397839840533,-28.1797797637707,-30.0097654617018,-30.2197638204808,
    -26.1397957070605,-28.3497784351632,-21.9298286096342,-22.2898257961124,
    -24.0198122755774,-24.1598111814301,-20.2298418957090,-13.6698931643274,
    -15.5598783933382,-15.9798751108962,-17.6798618248213,-13.2098967593829,
    -13.2298966030762,-15.0698822228539,-8.46993380408579,-8.66993224101816,
    -10.5199177826426,-5.87995404581162,-5.99995310797104,-7.76993927482250,
    -1.03999187204831,-1.17999077790097, 3.69997108324881, 1.92998491639735,
     1.96998460378383, 0.609995232643723, 7.39994216649762, 7.17994388587201,
     6.16995177936356, 6.13995201382370, 10.8699150472742, 10.3899187986365,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
    -57.9595470230003,-63.8495009906585, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
    -77.3093957962069, 0.0,-78.3393877464086,-79.6323776411763,
    -77.7593922793047,-81.2693648474678,-76.5094020484774, 0.0,
    -75.5294097075088,-70.1294519103349,-74.0994208834424,-71.1494439386900,
    -73.0894287769339,-66.1894827027672,-70.3994498001936,-65.3194895021114,
    -67.4694726991344, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0,-62.6295105253711,
    -66.8494775446440, 0.0,-64.5294956762286,-69.0494603509001,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0,-77.7593922793047,-78.6793850891936,
    -81.1193660197685,-81.5993622684062,-79.4693789150765,-82.3393564850560,
    -77.9693906380837,-78.4693867304146,-76.2894037678518,-79.1893811033712,
    -74.2094200237552,-75.0294136151779,-72.2994349510511,-75.4794100982757,
    -75.9094067376803,-78.9773827602228,-74.4094184606875,-75.6354088790830,
    -77.9353909038052,-75.9094067376803,-70.5094489405064,-67.7694703545329,
    -67.5294722302141,-62.0395151364206,-63.4395041949472,-60.4095278754218,
    -60.3095286569556,-54.5295738296102,-56.4495588241609,-53.0595853181573,
    -53.2695836769363,-47.2296308815787,-49.4896132189145,-49.9296097801657,
    -46.1296394784507,-48.7596189241114,-42.3696688641222,-43.0796633152321,
    -38.6696977808734,-41.4796758197731,-34.6497291985328,-35.6997209924277,
    -38.1297020011560,-38.9196958270388,-34.1597330280485,-36.8397120829422,
    -37.4097076281995,-30.8897585842043,-33.1597408433866,-33.8397355289567,
    -29.3497706198250,-31.6297528008540,-24.7498065703806,-25.6297996928830,
    -20.7898375191197,-23.2098186060013,-23.6898148546390,-26.0397964885943,
    -19.1098506488878,-19.8598447873842,-14.9398832388479,-17.4998632315822,
    -17.8598604180605,-20.2198419738624,-13.0498980098370, 0.0,
    -8.30993505453990,-10.8099155161945,-11.0599135623600,-4.04996834788045,
    -6.02995287351090,-6.66994787169448, 0.0,-3.73997077063528,
    -3.94996912941427, 3.15997530353142, 1.22999038713406, 0.869993200655801,
     5.97995326427781, 3.95996905126089, 0.0, 0.0,
     9.11992872411599, 8.85993075610391, 0.0, 0.0,
     12.2399043402609, 10.8699150472742, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
    -69.8794538641694,-71.4794413596284,-68.5094645711827,-71.9394377645728,
    -72.9394299492347,-67.8094700419194,-64.2094981771368, 0.0,
    -61.5295191222431,-57.9395471793070,-58.1995451473191,-51.8095950873300,
    -53.8695789877333,-50.0496088423251,-50.4496057161899,-43.8096576100352,
    -46.4696368212357,-42.2196700364229,-42.7996655035268,-45.7496424482792,
    -38.8396964522659,-39.7096896529217,-34.8497276354651,-38.0397027045364,
     0.0, 0.0, 0.0, 0.0,
    -30.5497612414192,-33.5297379517115,-26.1797953944470,-27.3197864849615,
    -29.9097662432357,-24.9298051636197,-25.6397996147296,-28.5397769502489,
    -21.1898343929844,-22.4098248582718,-16.7998687023189,-19.8598447873842,
     0.0, 0.0,-15.8298762831969,-16.8998679207851,
    -11.7399082479300,-14.3298880062042, 0.0, 0.0,
    -9.86992286261236,-10.8499152035810, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 6.27995091967636,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0}
  };

  excitEnResidNucl = 0.;
  recoilEnResidNucl = 0.;

  // A. Ferrari: first of all read isotopic data
  for (i = 0; i < 297; ++i) {
    residualNuclearMasses[i] = rmass[i];
  }
  for (i = 0; i < 11; ++i) {
    for (j = 0; j < 250; ++j) {
      isotopicData[i][j] = Wapstra[i][j];
    }
  }

  nucleonMass[0] = AMPROT;
  nucleonMass[1] = AMNEUT;
  nucleonMassSquared[0] = AMPROT * AMPROT;
  nucleonMassSquared[1] = AMNEUT * AMNEUT;
  nucleonBindingEn[0] = EBNDAV;
  nucleonBindingEn[1] = EBNDAV;

  return;

} // InitializeMuCapture


G4double G4MuonMinusCaptureAtRest::LevelDensity(G4int jz, G4int jn,
						G4double *aogmax,
						G4double *aogmin)
{
  // System generated locals
  G4double ret_val, d__1;

  // Local variables
  static G4double aa;
  static G4int ja;
  static G4double zz, temp;
  static G4bool lasmll;

  // ---------------------------------------------------------------------*
  //                                                                      *
  //     Created on 18 january 1993   by    Alfredo Ferrari & Paola Sala  *
  //                                                   Infn - Milan       *
  //                                                                      *
  //     Last change on 28-jan-93     by    Alfredo Ferrari               *
  //                                                                      *
  // ---------------------------------------------------------------------*

  lasmll = true;
  if (jz <= 0 || jn <= 0) {
    ret_val = (jz + jn) / B0;
    *aogmax = ret_val;
    *aogmin = ret_val;
    return ret_val;
  }
  ja = jn + jz;
  aa = (G4double) ja;
  zz = (G4double) jz;
  //  Standard EVAP parametrization for the level density
  d__1 = (aa - zz * 2.) / aa;
  temp = aa * (Y0 * (d__1 * d__1) + 1.) / B0;
  *aogmax = temp;
  *aogmin = temp;
  if (lasmll) {
    ret_val = temp;
  } else {
    ret_val = temp * ASMTOG;
    *aogmax *= ASMTOG;
    *aogmin *= ASMTOG;
  }
  return ret_val;

} // LevelDensity


void G4MuonMinusCaptureAtRest::GetCaptureIsotope(const G4Track& track)
{

  // Ask selector to choose the element

  G4Material * aMaterial = track.GetMaterial();
  G4Element* theElement  = pSelector->GetElement(aMaterial);
  targetCharge           = theElement->GetZ();
  targetAtomicMass       = theElement->GetN();

  // Calculate total capture velosity

  G4double lambdac  = pSelector->GetMuonCaptureRate(targetCharge,targetAtomicMass);
           lambdac += pSelector->GetMuonDecayRate(targetCharge,targetAtomicMass);

  // ===  Throw for capture  time.

  tDelay = -(G4double)log(G4UniformRand()) / lambdac;

  return;

} // GetCaptureIsotope


void G4MuonMinusCaptureAtRest::NuclearExcitation(G4double txi, G4double tyi,
						 G4double tzi, G4double *eke,
						 G4double *enu, G4double delm,
						 G4double emubnd)
{
  // System generated locals
  G4double d__1, d__2, d__3;

  // Local variables
  static G4double pnu, delm1, delm2, ekfer, pfcos, etest;

  // =======
  // === this routine interfaces to FLUKA routines to evaporate
  // === nucleons following muon capture
  // =======

  // =======
  // Calculates neutron and neutrino energies from muon capture with
  // Fermi motion
  // ==============
  *eke = 0.;
  delm1 =
    nucleonPotWell[1] + AMMUON + AMPROT - AMNEUT -
    nucleonPotWell[0] - emubnd;
  ekfer = 1e10;
  etest = 0.;
  //  Call FermiMotion for setting up anything for the Fermi motion
  //  Loop until total energy is conserved and neutron not fermi-blocked
  while(etest < ekfer || *eke <= 0.) {
    //  need proton fermi motion
    FermiMotion(1);
    //  energy and momentum conservation for mu + p > n + nu in potential wells.
    d__1 = nucleonFermiMomX;
    d__2 = nucleonFermiMomY;
    d__3 = nucleonFermiMomZ;
    delm2 =
      delm1 + (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) / (2.0 * AMPROT );
    pfcos =
      txi * nucleonFermiMomX + tyi * nucleonFermiMomY + tzi * nucleonFermiMomZ;
    pnu = pfcos - AMNEUT;
    d__1 = pnu;
    pnu += sqrt(d__1 * d__1 + delm1 * 2.0 * AMNEUT);
    *enu = abs(pnu);
    *eke = delm - *enu;
    etest = delm2 - *enu;
    //  now the neutron fermi energy
    FermiMotion(2);
    ekfer = nucleonFermiEn - AMNEUT;
  }
  return;

} // NuclearExcitation


void G4MuonMinusCaptureAtRest::DoMuCapture()
{

  // System generated locals
  G4int i__1, i__2;

  // Local variables
  static G4double emubnd;
  static G4double rndiso;
  static G4int ibtarm;
  static G4int nstak1;
  static G4int nstak2;
  static G4int k;
  static G4int is;
  static G4double cfe;
  static G4double eke;
  static G4double sfe, enu;
  static G4double txi, tyi, tzi;
  static G4double rndm[3];
  static G4double delm, polc, pols;
  static G4double bbtar, zztar, amntar1, amntar2;

  static G4double abuiso[304] = {
    9.998500000000000e-01, 1.500000000000000e-04, 9.999987000000000e-01
    , 1.300000000000000e-06, 9.258000000000000e-01, 7.420000000000000e-02
    , 1.000000000000000e+00, 8.022000000000000e-01, 1.978000000000000e-01
    , 9.889000000000000e-01, 1.110000000000000e-02, 9.963000000000000e-01
    , 3.700000000000000e-03, 9.975900000000000e-01, 2.040000000000000e-03
    , 3.700000000000000e-04, 1.000000000000000e+00, 9.113352585106910e-01
    , 8.840713792415629e-02, 2.576035651531538e-04, 1.000000000000000e+00
    , 7.870000000000000e-01, 1.117000000000000e-01, 1.013000000000000e-01
    , 1.000000000000000e+00, 9.220298268441600e-01, 4.704233810429386e-02
    , 3.092783505154639e-02, 1.000000000000000e+00, 9.500570034202050e-01
    , 4.220253215192912e-02, 7.600456027361642e-03, 1.400084005040303e-04
    , 7.553000000000000e-01, 2.447000000000000e-01, 9.960000000000000e-01
    , 3.370000000000000e-03, 6.300000000000000e-04, 9.310763482605570e-01
    , 6.880564206264915e-02, 1.180096767934971e-04, 9.677570699861240e-01
    , 2.079319604300388e-02, 6.448586002775087e-03, 1.856713616898090e-03
    , 1.796819629256216e-03, 1.347614721942162e-03, 1.000000000000000e+00
    , 7.393999999999999e-01, 7.929999999999998e-02, 7.279999999999999e-02
    , 5.509999999999999e-02, 5.339999999999999e-02, 9.976000000000000e-01
    , 2.400000000000000e-03, 8.376000000000000e-01, 9.550000000000000e-02
    , 4.310000000000000e-02, 2.380000000000000e-02, 1.000000000000000e+00
    , 9.166000000000000e-01, 5.820000000000000e-02, 2.190000000000000e-02
    , 3.300000000000000e-03, 1.000000000000000e+00, 6.827400000000000e-01
    , 2.609500000000000e-01, 3.593000000000000e-02, 1.134000000000000e-02
    , 9.039999999999999e-03, 6.909000000000000e-01, 3.091000000000000e-01
    , 4.868473895582330e-01, 2.792168674698800e-01, 1.864457831325300e-01
    , 4.126506024096385e-02, 6.224899598393573e-03, 6.040000000000000e-01
    , 3.960000000000000e-01, 3.653634636536350e-01, 2.742725727427260e-01
    , 2.051794820517950e-01, 7.759224077592240e-02, 7.759224077592240e-02
    , 1.000000000000000e+00, 4.982000000000000e-01, 2.352000000000000e-01
    , 9.190000000000000e-02, 9.020000000000000e-02, 7.580000000000001e-02
    , 8.699999999999999e-03, 5.054000000000000e-01, 4.946000000000000e-01
    , 5.690000000000000e-01, 1.737000000000000e-01, 1.156000000000000e-01
    , 1.155000000000000e-01, 2.270000000000000e-02, 3.499999999999999e-03
    , 7.215000000000000e-01, 2.785000000000000e-01, 8.256000000000000e-01
    , 9.859999999999999e-02, 7.020000000000000e-02, 5.600000000000001e-03
    , 1.000000000000000e+00, 5.145999999999999e-01, 1.740000000000000e-01
    , 1.711000000000000e-01, 1.123000000000000e-01, 2.799999999999999e-02
    , 1.000000000000000e+00, 2.378000000000000e-01, 1.653000000000000e-01
    , 1.584000000000000e-01, 1.572000000000000e-01, 9.630000000000001e-02
    , 9.460000000000000e-02, 9.039999999999999e-02, 1.000000000000000e+00
    , 3.161632326465290e-01, 1.858371674334870e-01, 1.707341468293660e-01
    , 1.272254450890180e-01, 1.262252450490100e-01, 5.511102220444088e-02
    , 1.870374074814963e-02, 1.000000000000000e+00, 2.732726727327270e-01
    , 2.670732926707330e-01, 2.222777722227780e-01, 1.180881911808820e-01
    , 1.096890310968900e-01, 9.599040095990399e-03, 5.182000000000000e-01
    , 4.818000000000000e-01, 2.885711428857110e-01, 2.406759324067590e-01
    , 1.274872512748730e-01, 1.238876112388760e-01, 1.225877412258770e-01
    , 7.579242075792421e-02, 1.219878012198780e-02, 8.799120087991200e-03
    , 9.572000000000001e-01, 4.280000000000000e-02, 3.285000000000000e-01
    , 2.403000000000000e-01, 1.430000000000000e-01, 8.580000000000000e-02
    , 7.610000000000000e-02, 5.940000000000000e-02, 4.720000000000000e-02
    , 9.599999999999999e-03, 6.600000000000000e-03, 3.500000000000000e-03
    , 5.725000000000000e-01, 4.275000000000000e-01, 3.448034480344800e-01
    , 3.179031790317900e-01, 1.871018710187100e-01, 6.990069900699007e-02
    , 4.610046100461004e-02, 2.460024600246002e-02, 8.700087000870008e-03
    , 8.900089000890008e-04, 1.000000000000000e+00, 2.688838669679820e-01
    , 2.643841369517830e-01, 2.117872927624340e-01, 1.043937363758170e-01
    , 8.869467831930082e-02, 4.079755214687118e-02, 1.919884806911585e-02
    , 9.599424034557925e-04, 8.999460032398054e-04, 1.000000000000000e+00
    , 7.166143322866459e-01, 1.132022640452810e-01, 7.810156203124062e-02
    , 6.590131802636053e-02, 2.420048400968019e-02, 1.010020200404008e-03
    , 9.700194003880077e-04, 9.991100000000001e-01, 8.899999999999999e-04
    , 8.848619403358230e-01, 1.107077495424680e-01, 2.500175012250857e-03
    , 1.930135109457662e-03, 1.000000000000000e+00, 2.700199203187250e-01
    , 2.375498007968130e-01, 1.754980079681270e-01, 1.212151394422310e-01
    , 8.266932270916334e-02, 5.707171314741035e-02, 5.597609561752987e-02
    , 1.000000000000000e+00, 2.672000000000000e-01, 2.271000000000000e-01
    , 1.497000000000000e-01, 1.383000000000000e-01, 1.124000000000000e-01
    , 7.439999999999999e-02, 3.089999999999999e-02, 5.218000000000000e-01
    , 4.782000000000000e-01, 2.487000000000000e-01, 2.190000000000000e-01
    , 2.047000000000000e-01, 1.568000000000000e-01, 1.473000000000000e-01
    , 2.149999999999999e-02, 2.000000000000000e-03, 1.000000000000000e+00
    , 2.795579452788640e-01, 2.532687843495170e-01, 2.477133390210510e-01
    , 1.872978710739870e-01, 2.271780321819011e-02, 8.928394277891312e-03
    , 5.158627805003868e-04, 1.000000000000000e+00, 3.341133645345810e-01
    , 2.707108284331370e-01, 2.294091763670550e-01, 1.488059522380900e-01
    , 1.560062402496100e-02, 1.360054402176087e-03, 1.000000000000000e+00
    , 3.184159207960400e-01, 2.182109105455270e-01, 1.613080654032700e-01
    , 1.431071553577680e-01, 1.273063653182660e-01, 3.030151507575378e-02
    , 1.350067503375169e-03, 9.741000000000000e-01, 2.590000000000000e-02
    , 3.523647635236480e-01, 2.713728627137290e-01, 1.849815018498150e-01
    , 1.374862513748620e-01, 5.199480051994800e-02, 1.799820017998200e-03
    , 9.998770003689990e-01, 1.229996310011070e-04, 3.064000000000000e-01
    , 2.841000000000000e-01, 2.641000000000000e-01, 1.440000000000000e-01
    , 1.400000000000000e-03, 6.250000000000000e-01, 3.750000000000000e-01
    , 4.035555555555560e-01, 2.666666666666670e-01, 1.626262626262630e-01
    , 1.343434343434340e-01, 1.656565656565656e-02, 1.606060606060606e-02
    , 1.818181818181818e-04, 6.260000000000000e-01, 3.740000000000000e-01
    , 3.379908742463950e-01, 3.289911172398350e-01, 2.529931691844320e-01
    , 7.209805335255949e-02, 7.799789405686048e-03, 1.269965710925805e-04
    , 1.000000000000000e+00, 2.979821210727360e-01, 2.312861228326300e-01
    , 1.683898966062040e-01, 1.321920684758910e-01, 1.001939883606980e-01
    , 6.849589024658520e-02, 1.459912405255685e-03, 7.050000000000000e-01
    , 2.950000000000000e-01, 5.231046209241850e-01, 2.360472094418880e-01
    , 2.260452090418080e-01, 1.480296059211842e-02, 1.000000000000000e+00
    , 1.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00
    , 1.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00
    , 1.000000000000000e+00, 1.000000000000000e+00, 9.922542601569140e-01
    , 7.196042176802758e-03, 5.496976662835440e-04, 1.000000000000000e+00
    , 1.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00
    , 1.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00
    , 1.000000000000000e+00
  };

  static G4int isomnm[304] = {
    1, 2, 4, 3, 7, 6, 9, 11, 10, 12, 13, 14, 15, 16, 18,
    17, 19, 20, 22, 21, 23, 24, 26, 25, 27, 28, 29, 30, 31, 32,
    34, 33, 36, 35, 37, 40, 36, 38, 39, 41, 40, 40, 44, 42, 46,
    48, 43, 45, 48, 46, 47, 49, 50, 51, 50, 52, 53, 50, 54, 55,
    56, 54, 57, 58, 59, 58, 60, 62, 61, 64, 63, 65, 64, 66, 68,
    67, 70, 69, 71, 74, 72, 70, 73, 76, 75, 80, 78, 82, 76, 77,
    74, 79, 81, 84, 86, 82, 83, 80, 78, 85, 87, 88, 86, 87, 84,
    89, 90, 94, 92, 91, 96, 93, 98, 96, 92, 95, 100, 97, 94, 97,
    102, 104, 101, 99, 100, 96, 98, 103, 106, 108, 105, 110, 104, 102, 107,
    109, 114, 112, 111, 110, 113, 116, 106, 108, 115, 113, 120, 118, 116, 119,
    117, 124, 112, 112, 114, 115, 121, 123, 130, 128, 126, 125, 124, 122, 123,
    120, 127, 132, 129, 131, 134, 136, 130, 128, 124, 126, 133, 138, 137, 136,
    135, 134, 130, 132, 139, 138, 140, 142, 138, 136, 141, 142, 144, 146, 143,
    145, 148, 150, 145, 152, 154, 147, 149, 148, 150, 144, 153, 151, 158, 160,
    156, 157, 155, 154, 152, 159, 164, 162, 163, 161, 160, 158, 156, 165, 166,
    168, 167, 170, 164, 162, 169, 174, 172, 173, 171, 176, 170, 168, 175, 176,
    180, 178, 177, 179, 176, 174, 181, 180, 184, 186, 182, 183, 180, 187, 185,
    192, 190, 189, 188, 187, 186, 184, 193, 191, 195, 194, 196, 198, 192, 190,
    197, 202, 200, 199, 201, 198, 204, 196, 205, 203, 208, 206, 207, 204, 209,
    209, 210, 222, 223, 226, 227, 232, 231, 238, 235, 234, 237, 244, 243, 247,
    247, 251, 254, 257
  };

  static G4int isondx[200] = {
    1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14,
    16, 17, 17, 18, 20, 21, 21, 22, 24, 25, 25, 26, 28, 29, 29,
    30, 33, 34, 35, 36, 38, 39, 41, 42, 47, 48, 48, 49, 53, 54,
    55, 56, 59, 60, 60, 61, 64, 65, 65, 66, 70, 71, 72, 73, 77,
    78, 79, 80, 84, 85, 85, 86, 91, 92, 93, 94, 99, 100, 101, 102,
    105, 106, 106, 107, 111, 112, 112, 113, 119, 120, 120, 121, 127, 128, 128,
    129, 134, 135, 136, 137, 144, 145, 146, 147, 156, 157, 158, 159, 166, 167,
    167, 168, 176, 177, 177, 178, 184, 185, 186, 187, 190, 191, 191, 192, 198,
    199, 199, 200, 206, 207, 208, 209, 215, 216, 216, 217, 223, 224, 224, 225,
    230, 231, 231, 232, 238, 239, 240, 241, 246, 247, 248, 249, 253, 254, 255,
    256, 262, 263, 264, 265, 270, 271, 271, 272, 278, 279, 280, 281, 284, 285,
    285, 286, 286, 287, 287, 288, 288, 289, 289, 290, 290, 291, 291, 292, 292,
    293, 293, 294, 296, 297, 297, 298, 298, 299, 299, 300, 300, 301, 301, 302,
    302, 303, 303, 304, 304
  };

  // =======
  // === this routine interfaces to FLUKA routines to evaporate
  // === nucleons following muon capture
  // =======

  nGkine = 0;            // number of generated secondary particles
  nSecPart = 0;
  nallFragm = 0;
  // ==      choose random direction
  RanPolarAng(&polc, &pols);
  RanAzimuthalAng(&sfe, &cfe);
  txi = cfe * pols;
  tyi = sfe * pols;
  tzi = polc;
  chargeTarget = G4int(targetCharge);
  //  Choice of the mass number of the target nucleus: use the input
  rndm[0] = G4UniformRand();
  rndm[1] = G4UniformRand();
  rndiso = rndm[0];
  //  Loop on the stable isotopes
  i__1 = isondx[(chargeTarget << 1) - 1] - 1;
  for (is = isondx[(chargeTarget << 1) - 2]; is <= i__1; ++is) {
    rndiso -= abuiso[is - 1];
    if (rndiso < 0.) {
      break;
    }
  }
  if (rndiso >= 0.) {
    is = isondx[(chargeTarget << 1) - 1];
  }
  atMassTarget = isomnm[is - 1];
  bbtar = (G4double) atMassTarget;
  zztar = (G4double) chargeTarget;
  if (atMassTarget != 1) {
    //  The following should be Done with the proper mass of the nuclide
    //        AMNTAR = BBTAR * AMUC12
    massTarget = bbtar * AMUAMU + GetIsotopicMass(bbtar, zztar) * .001;
    amntar1 = massTarget - zztar * AMELEC;
    CascadeCorrection(zztar, bbtar);
  }
  enu = 0.;
  eke = 0.;
  //  ---------------------
  //   mass of Z-1,A
  //  ---------------------
  --chargeTarget;
  ibtarm = atMassTarget;
  bbtar = (G4double) atMassTarget;
  zztar = (G4double) chargeTarget;
  if (atMassTarget == 1) {
    switch (chargeTarget) {
    case 0:
      amntar2 = massNeutron;
      break;
    case 1:
      amntar2 = massProton;
      break;
    default:
      G4cout << "? This cannot happen!" << G4endl;
      exit(1);
    }
  } else {
    //  The following should be Done with the proper mass of the nuclide
    //        AMNTAR = BBTAR * AMUC12
    massTarget = bbtar * AMUAMU + GetIsotopicMass(bbtar, zztar) * .001;
    amntar2 = massTarget - zztar * AMELEC;
  }
  //  ---------------------
  //   Calculate the nuclear excitation energy gkin and neutrino energy ENU
  //  ---------------------
  emubnd = zeff2 * AMMUON * FSCTO2 * .5;
  delm = amntar1 - amntar2 + AMMUON - emubnd;
  // --------------------------------
  //   Takes account of fermi motion of the proton and available phase space
  //  -------------------------------
  //   (for nuclei with A > 1 only!)
  //  -------------------------------
  if (atMassTarget > 1) {
    NuclearExcitation(txi, tyi, tzi, &eke, &enu, delm, emubnd);
  }
  recoilEnResidNucl = 0.;
  excitEnResidNucl = eke;
  if (eke > 1e-4) {
    momXResidNucl = txi;
    momYResidNucl = tyi;
    momZResidNucl = tzi;
    totMomResidNucl = enu;
    atMassCurrent = bbtar;
    chargeCurrent = zztar;
    massResidNucl = massTarget;
    totEnResidNucl = massResidNucl + eke;
    MuEvaporation();
  } else {
    if (atMassTarget == 1) {
      // mu- + p --> n + neutrino for mu- and p at rest
      nSecPart = 1;
      enu = 0.5 * (AMPROT+AMMUON) * (1.0-AMMRED*AMMRED);
      Secondaries[0].SetZero();
      Secondaries[0].SetMass( massNeutron );
      Secondaries[0].SetMomentumAndUpdate( txi*enu, tyi*enu, tzi*enu );
      Secondaries[0].SetParticleDef( pdefNeutron );
    } else {
      //  no neutrons or gammas below 100 KeV
      nstak1 = 0;
      nSecPart = 0;
    }
  }
  i__1 = nSecPart, i__2 = MXGKIN - nGkine;
  nstak1 = G4std::min(i__1,i__2);
  if (nSecPart > nstak1) {
    G4cout << " **** FLUFIN: Stack overflow, "
	 << nSecPart - nstak1 << " particles lost" << G4endl;
  }
  // ==     Put neutrino on the stack IF DESIRED.
  ++nGkine;
  Gkin[nGkine-1].SetZero();
  Gkin[nGkine-1].SetMass( massNeutrinoE  );
  Gkin[nGkine-1].SetMomentumAndUpdate( -txi * enu, -tyi * enu, -tzi * enu  );
  Gkin[nGkine-1].SetParticleDef( pdefNeutrinoE );
  Gkin[nGkine-1].SetTOF( tDelay );
  // ==
  i__1 = nstak1;
  for (k = 1; k <= i__1; ++k) {
    ++nGkine;
    Gkin[nGkine - 1] = Secondaries[k - 1];
    Gkin[nGkine - 1].SetTOF( tDelay );
  }
  i__1 = nallFragm, i__2 = MXGKIN - nGkine;
  nstak2 = G4std::min(i__1,i__2);
  if (nallFragm > nstak2) {
    G4cout << " **** FLUFIN: Stack overflow, "
	 << nallFragm - nstak2 << " heavy particles lost" << G4endl;
  }
  i__1 = nstak2;
  for (k = 1; k <= i__1; ++k) {
    ++nGkine;
    Gkin[nGkine - 1] = Fragments[k - 1];
    Gkin[nGkine - 1].SetTOF( tDelay );
  }
  return;

} // DoMuCapture


void G4MuonMinusCaptureAtRest::MuEvaporation()
{
  // System generated locals
  G4int i__1;
  G4double d__1;

  // Local variables
  static G4double eotest, etevap;
  static G4int ip, jp, itemp;

  // ---------------------------------------------------------------------*
  //                                                                      *
  //  EVent EVAPoration: this routine is used to steer both the evapora-  *
  //  tion, the high energy fission, possibly a future fragmentation      *
  //  and the gamma deexcitation routines                                 *
  //                                                                      *
  //  Created  on  15  may  1991   by   Alfredo Ferrari & Paola Sala      *
  //                                             INFN - Milan             *
  //                                                                      *
  //  Last change  on 19-apr-93    By   Alfredo Ferrari, INFN - Milan     *
  //                                                                      *
  // ---------------------------------------------------------------------*

  //  The initial excitation energy, mass and charge of the nucleus are
  //  put into Ex, Apr, Zpr (common Hetc5)
  d__1 = excitEnResidNucl * 1000;
  iniExcitEnEvapNucl = G4std::max(d__1,ANGLGB);
  atMassEvapNucl = atMassCurrent;
  chargeEvapNucl = chargeCurrent;
  //  Ammres is the atomic mass of the residual nucleus
  //  Reset accumulators for the energy conservation check (they are only
  //  local)
  eotest = massResidNucl + excitEnResidNucl + recoilEnResidNucl;
  etevap = 0.;
  if (totMomResidNucl > 0.) {
    //  Set the variables recording the recoil direction of the residual
    //  nucleus:
    dirCosXResidNucl = momXResidNucl / totMomResidNucl;
    dirCosYResidNucl = momYResidNucl / totMomResidNucl;
    dirCosZResidNucl = momZResidNucl / totMomResidNucl;
  } else {
    //  It can happen for pion capture for example that ptres=0
    //  ( it is always 0 if no "direct" particle is emitted )
    dirCosXResidNucl = 0.;
    dirCosYResidNucl = 0.;
    dirCosZResidNucl = 1.;
  }
  recoilEnEvapNucl = recoilEnResidNucl * 1e3;
  Erup();
  //  Add to the secondary stack the evaporated neutrons
  i__1 = nVariousFragm[0];
  for (ip = 1; ip <= i__1; ++ip) {
    ++nSecPart;
    Secondaries[nSecPart - 1] = Evaporates[ip - 1];
    etevap += Evaporates[ip - 1].GetEnergy();
  }
  //  Add to the secondary stack the evaporated protons
  i__1 = nVariousFragm[1];
  for (ip = 1; ip <= i__1; ++ip) {
    ++nSecPart;
    Secondaries[nSecPart - 1] = Evaporates[ip - 1 + MXEVAP];
    etevap += Evaporates[ip - 1 + MXEVAP].GetEnergy();
  }
  //  Loop over the particle types:
  for (jp = 3; jp <= 6; ++jp) {
    i__1 = nVariousFragm[jp - 1];
    for (ip = 1; ip <= i__1; ++ip) {
      ++nallFragm;
      itemp = ip - 1 + (jp - 1)*MXEVAP;
      Fragments[nallFragm - 1] = Evaporates[itemp];
      etevap += Evaporates[itemp].GetEnergy();
    }
  }
  //  Normal evaporation:
  atMassCurrent = atMassEvapNucl;
  chargeCurrent = chargeEvapNucl;
  atMassResidNucl = NINT(atMassCurrent);
  chargeResidNucl = NINT(chargeCurrent);
  //  Ammres is the atomic mass of the residual nucleus
  //  Check the residual nucleus:
  if (atMassResidNucl == 0) {
    massResidNucl = 0.;
    excitEnResidNucl = 0.;
    recoilEnResidNucl = 0.;
    totMomResidNucl = 0.;
    momXResidNucl = 0.;
    momYResidNucl = 0.;
    momZResidNucl = 0.;
    totEnResidNucl = 0.;
  } else {
    massResidNucl =
      atMassCurrent * AMUAMU +
      GetIsotopicMass(atMassCurrent, chargeCurrent) * .001;
    excitEnResidNucl = finExcitEnEvapNucl * .001;
    d__1 = recoilEnEvapNucl * .001;
    recoilEnResidNucl = G4std::max(d__1,0.);
    totMomResidNucl =
      sqrt(recoilEnResidNucl *
	   (recoilEnResidNucl + (massResidNucl + excitEnResidNucl) * 2.));
    momXResidNucl = totMomResidNucl * dirCosXResidNucl;
    momYResidNucl = totMomResidNucl * dirCosYResidNucl;
    momZResidNucl = totMomResidNucl * dirCosZResidNucl;
    totEnResidNucl = massResidNucl + excitEnResidNucl + recoilEnResidNucl;
    kinEnResidNucl = recoilEnResidNucl;
  }
  etevap += totEnResidNucl;
  if ((d__1 = etevap - eotest, abs(d__1)) / eotest > .05) {
    G4cout << " Evevap: failure in energy conservation!!"
	 << " " << etevap << " " << eotest << " " << etevap - eotest
	 << " " << totEnResidNucl << G4endl;
  }
  //  Check if the deexcitation module have to be called
  EvaporationDeexcitation();
  return;

} // MuEvaporation


void G4MuonMinusCaptureAtRest::RanPolarAng(G4double *cs, G4double *si)
{
  // Local variables
  static G4double rndm[1];

  // --------------------------------------------------
  // *** RANDOM CHOICE OF ANGLE TETA (CS = COS(TETA),SI = SIN(TETA)
  // --------------------------------------------------
  rndm[0] = G4UniformRand();
  *cs = rndm[0] * 2. - 1.;
  *si = sqrt((1. - *cs) * (*cs + 1.));
  return;

} // RanPolarAng


G4double G4MuonMinusCaptureAtRest::NuclearBindingEnergy(G4double a1,
							G4double z1,
							G4double a2,
							G4double z2)
{
  // System generated locals
  G4double ret_val;

  // Local variables
  static G4int n1, n2, ka1, ka2, iz1, iz2, kz1, kz2, izz1, izz2;
  static G4double enrg1, enrg2;

  ka1 = NINT(a1);
  kz1 = NINT(z1);
  ka2 = NINT(a2);
  kz2 = NINT(z2);
  n1 = ka1 - kz1;
  n2 = ka2 - kz2;
  if (n1 <= 0) {
    G4cout << "STOP - 105" << G4endl;
    exit(1);
  }
  if (n2 <= 0) {
    G4cout << "STOP - 106" << G4endl;
    exit(1);
  }
  izz1 = NINT(isotopicData[0][ka1 - 1]);
  izz2 = NINT(isotopicData[0][ka2 - 1]);
  if (kz1 < izz1 || kz1 > izz1 + 9) {
    enrg1 = CalculateIsotopicMass(a1, z1);
  } else {
    iz1 = kz1 - izz1 + 2;
    enrg1 = isotopicData[iz1 - 1][ka1 - 1];
    if (enrg1 == 0. && (ka1 != 12 || kz1 != 6)) {
      enrg1 = CalculateIsotopicMass(a1, z1);
    }
  }
  if (kz2 < izz2 || kz2 > izz2 + 9) {
    enrg2 = CalculateIsotopicMass(a2, z2);
  } else {
    iz2 = kz2 - izz2 + 2;
    enrg2 = isotopicData[iz2 - 1][ka2 - 1];
    if (enrg2 == 0. && (ka2 != 12 || kz2 != 6)) {
      enrg2 = CalculateIsotopicMass(a2, z2);
    }
  }
  ret_val = enrg1 - enrg2;
  return ret_val;

} // NuclearBindingEnergy


void G4MuonMinusCaptureAtRest::RanDirCos(G4double *wx, G4double *wy,
					 G4double *wz)
{
  // Local variables
  static G4double x, y, z, x2, y2, z2, cfe, sfe;
  static G4double rndm[2];

  // ********************************************************************
  //     VERSION JUNE 81 BY             PERTTI AARNIO
  //     LAST CHANGE 22. JUNE 81 BY     PERTTI AARNIO
  //                                    HELSINKI UNIVERSITY OF
  //                                    TECHNOLOGY, FINLAND
  //
  //     SUBROUTINE OF FLUKA TO GIVE THE DIRECTION COSINES OF RANDOM
  //     UNIFORM (ISOTROPIC) DIRECTION IN THREE DIMENSIONAL SPACE.
  // ********************************************************************

  do {
    rndm[0] = G4UniformRand();
    rndm[1] = G4UniformRand();
    x = rndm[0] * 2. - 1.;
    y = rndm[1];
    x2 = x * x;
    y2 = y * y;
  } while (x2 + y2 > 1.);
  cfe = (x2 - y2) / (x2 + y2);
  sfe = x * 2. * y / (x2 + y2);
  rndm[0] = G4UniformRand();
  z = rndm[0];
  z2 = z * z;
  *wz = sqrt(z - z2);
  *wx = *wz * 2. * cfe;
  *wy = *wz * 2. * sfe;
  *wz = z * 2. - 1.;
  return;

} // RanDirCos


void G4MuonMinusCaptureAtRest::RanAzimuthalAng(G4double *sfe, G4double *cfe)
{
  static G4double x, y, x2, y2;
  static G4double rndm[2];

  // ********************************************************************
  //     VERSION JUNE 81 BY             PERTTI AARNIO
  //     LAST CHANGE 11. DECEMBER 85 BY PERTTI AARNIO
  //                                    HELSINKI UNIVERSITY OF
  //                                    TECHNOLOGY, FINLAND
  //
  //     SUBROUTINE OF FLUKA TO GIVE SINE AND COSINE OF AN
  //     RANDOM ANGLE UNIFORMLY DISTRIBUTED BETWEEN 0 AND 2*PI
  // ********************************************************************

  do {
    rndm[0] = G4UniformRand();
    rndm[1] = G4UniformRand();
    x = rndm[0] * 2. - 1.;
    y = rndm[1];
    x2 = x * x;
    y2 = y * y;
  } while (x2 + y2 > 1.);
  *cfe = (x2 - y2) / (x2 + y2);
  *sfe = x * 2. * y / (x2 + y2);
  return;

} // RanAzimuthalAng

