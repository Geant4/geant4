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
// Geant4 Header : G4HadronElastic
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  

#include "G4HadronElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsModelCatalog.hh"


G4HadronElastic::G4HadronElastic(const G4String& name) 
  : G4HadronicInteraction(name), secID(-1)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  lowestEnergyLimit= 1.e-6*eV;
  pLocalTmax  = 0.0;
  nwarn = 0;

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();

  secID = G4PhysicsModelCatalog::GetModelID( "model_" + name );
}

G4HadronElastic::~G4HadronElastic()
{}


void G4HadronElastic::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4HadronElastic is the base class for all hadron-nucleus\n" 
          << "elastic scattering models except HP.\n" 
          << "By default it uses the Gheisha two-exponential momentum\n"
	  << "transfer parameterization.  The model is fully relativistic\n"
	  << "as opposed to the original Gheisha model which was not.\n"
	  << "This model may be used for all long-lived hadrons at all\n"
	  << "incident energies but fit the data only for relativistic scattering.\n";
}

G4HadFinalState* G4HadronElastic::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  // no scattering below the limit
  if(ekin <= lowestEnergyLimit) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(0.,0.,1.);
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = std::sqrt(ekin*(ekin + 2.0*m1));

  if (verboseLevel>1) {
    G4cout << "G4HadronElastic: " 
	   << aParticle->GetDefinition()->GetParticleName() 
	   << " Plab(GeV/c)= " << plab/GeV  
	   << " Ekin(MeV) = " << ekin/MeV 
	   << " scattered off Z= " << Z 
	   << " A= " << A 
	   << G4endl;
  }

  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double e1 = m1 + ekin;
  G4LorentzVector lv(0.0,0.0,plab,e1+mass2);
  G4ThreeVector bst = lv.boostVector();
  G4double momentumCMS = plab*mass2/std::sqrt(m1*m1 + mass2*mass2 + 2.*mass2*e1);

  pLocalTmax = 4.0*momentumCMS*momentumCMS;

  // Sampling in CM system
  G4double t = SampleInvariantT(theParticle, plab, Z, A);

  if(t < 0.0 || t > pLocalTmax) {
    // For the very rare cases where cos(theta) is greater than 1 or smaller than -1,
    // print some debugging information via a "JustWarning" exception, and resample
    // using the default algorithm
#ifdef G4VERBOSE
    if(nwarn < 2) {
      G4ExceptionDescription ed;
      ed << GetModelName() << " wrong sampling t= " << t << " tmax= " << pLocalTmax
	 << " for " << aParticle->GetDefinition()->GetParticleName() 
	 << " ekin=" << ekin << " MeV" 
	 << " off (Z,A)=(" << Z << "," << A << ") - will be resampled" << G4endl;
      G4Exception( "G4HadronElastic::ApplyYourself", "hadEla001", JustWarning, ed);
      ++nwarn;
    }
#endif
    t = G4HadronElastic::SampleInvariantT(theParticle, plab, Z, A);
  }

  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/pLocalTmax;

  if (cost > 1.0) { cost = 1.0; }
  else if(cost < -1.0) { cost = -1.0; } 

  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  if (verboseLevel>1) {
    G4cout << " t= " << t << " tmax(GeV^2)= " << pLocalTmax/(GeV*GeV) 
	   << " Pcms(GeV)= " << momentumCMS/GeV << " cos(t)=" << cost 
	   << " sin(t)=" << sint << G4endl;
  }
  G4LorentzVector nlv1(momentumCMS*sint*std::cos(phi),
		       momentumCMS*sint*std::sin(phi),
                       momentumCMS*cost,
		       std::sqrt(momentumCMS*momentumCMS + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) {
    G4cout <<"G4HadronElastic: m= " << m1 << " Efin(MeV)= " << eFinal 
	   << " 4-M Final: " << nlv1 
	   << G4endl;
  }

  if(eFinal <= 0.0) { 
    theParticleChange.SetMomentumChange(0.0,0.0,1.0);
    theParticleChange.SetEnergyChange(0.0);
  } else {
    theParticleChange.SetMomentumChange(nlv1.vect().unit());
    theParticleChange.SetEnergyChange(eFinal);
  }
  lv -= nlv1;
  G4double erec =  std::max(lv.e() - mass2, 0.0);
  if (verboseLevel > 1) {
    G4cout << "Recoil: " <<" m= " << mass2 << " Erec(MeV)= " << erec
	   << " 4-mom: " << lv 
	   << G4endl;
  }
 
  // the recoil is created if kinetic energy above the threshold
  if(erec > GetRecoilEnergyThreshold()) {
    G4ParticleDefinition * theDef = nullptr;
    if(Z == 1 && A == 1)       { theDef = theProton; }
    else if (Z == 1 && A == 2) { theDef = theDeuteron; }
    else if (Z == 1 && A == 3) { theDef = G4Triton::Triton(); }
    else if (Z == 2 && A == 3) { theDef = G4He3::He3(); }
    else if (Z == 2 && A == 4) { theDef = theAlpha; }
    else {
      theDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
    }
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, lv.vect().unit(), erec);
    theParticleChange.AddSecondary(aSec, secID);
  } else {
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}

// sample momentum transfer in the CMS system 
G4double 
G4HadronElastic::SampleInvariantT(const G4ParticleDefinition* part,
				  G4double mom, G4int, G4int A)
{
  const G4double plabLowLimit = 400.0*CLHEP::MeV;
  const G4double GeV2 = GeV*GeV;
  const G4double z07in13 = std::pow(0.7, 0.3333333333);
  const G4double numLimit = 18.;

  G4int pdg = std::abs(part->GetPDGEncoding());
  G4double tmax = pLocalTmax/GeV2;

  G4double aa, bb, cc, dd;
  G4Pow* g4pow = G4Pow::GetInstance();
  if (A <= 62) {
    if (pdg == 211){ //Pions
      if(mom >= plabLowLimit){     //High energy
	bb = 14.5*g4pow->Z23(A);/*14.5*/
	dd = 10.;
	cc = 0.075*g4pow->Z13(A)/dd;//1.4
	//aa = g4pow->powZ(A, 1.93)/bb;//1.63
	aa = (A*A)/bb;//1.63
      } else {                       //Low energy
	bb = 29.*z07in13*z07in13*g4pow->Z23(A);
	dd = 15.;
	cc = 0.04*g4pow->Z13(A)/dd;//1.4
	aa = g4pow->powZ(A, 1.63)/bb;//1.63
      }
    } else { //Other particles
      bb = 14.5*g4pow->Z23(A);
      dd = 20.;
      aa = (A*A)/bb;//1.63
      cc = 1.4*g4pow->Z13(A)/dd;          
    }
      //===========================
  } else { //(A>62)
    if (pdg == 211) {
      if(mom >= plabLowLimit){ //high
	bb = 60.*z07in13*g4pow->Z13(A);//60
	dd = 30.;
	aa = 0.5*(A*A)/bb;//1.33
	cc = 4.*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
      } else { //low
	bb = 120.*z07in13*g4pow->Z13(A);//60
	dd = 30.;
	aa = 2.*g4pow->powZ(A,1.33)/bb;
	cc = 4.*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
      }
    } else {
      bb = 60.*g4pow->Z13(A);
      dd = 25.;
      aa = g4pow->powZ(A,1.33)/bb;//1.33
      cc = 0.2*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
    }
  }
  G4double q1 = 1.0 - G4Exp(-std::min(bb*tmax, numLimit));
  G4double q2 = 1.0 - G4Exp(-std::min(dd*tmax, numLimit));
  G4double s1 = q1*aa;
  G4double s2 = q2*cc;
  if((s1 + s2)*G4UniformRand() < s2) {
    q1 = q2;
    bb = dd;
  }
  return -GeV2*G4Log(1.0 - G4UniformRand()*q1)/bb;
}

//////////////////////////////////////////////
//
// Cofs for s-,c-,b-particles ds/dt slopes

G4double G4HadronElastic::GetSlopeCof(const G4int pdg )
{
  // The input parameter "pdg" should be the absolute value of the PDG code
  // (i.e. the same value for a particle and its antiparticle).

  G4double coeff = 1.0;

  // heavy barions

  static const G4double  lBarCof1S  = 0.88;
  static const G4double  lBarCof2S  = 0.76;
  static const G4double  lBarCof3S  = 0.64;
  static const G4double  lBarCof1C  = 0.784378;
  static const G4double  lBarCofSC  = 0.664378;
  static const G4double  lBarCof2SC = 0.544378;
  static const G4double  lBarCof1B  = 0.740659;
  static const G4double  lBarCofSB  = 0.620659;
  static const G4double  lBarCof2SB = 0.500659;
  
  if( pdg == 3122 || pdg == 3222 ||  pdg == 3112 || pdg == 3212  )
  {
    coeff = lBarCof1S; // Lambda, Sigma+, Sigma-, Sigma0

  } else if( pdg == 3322 || pdg == 3312   )
  {
    coeff = lBarCof2S; // Xi-, Xi0
  }
  else if( pdg == 3324)
  {
    coeff = lBarCof3S; // Omega
  }
  else if( pdg == 4122 ||  pdg == 4212 ||   pdg == 4222 ||   pdg == 4112   )
  {
    coeff = lBarCof1C; // LambdaC+, SigmaC+, SigmaC++, SigmaC0
  }
  else if( pdg == 4332 )
  {
    coeff = lBarCof2SC; // OmegaC
  }
  else if( pdg == 4232 || pdg == 4132 )
  {
    coeff = lBarCofSC; // XiC+, XiC0
  }
  else if( pdg == 5122 || pdg == 5222 || pdg == 5112 || pdg == 5212    )
  {
    coeff = lBarCof1B; // LambdaB, SigmaB+, SigmaB-, SigmaB0
  }
  else if( pdg == 5332 )
  {
    coeff = lBarCof2SB; // OmegaB-
  }
  else if( pdg == 5132 || pdg == 5232 ) // XiB-, XiB0
  {
    coeff = lBarCofSB;
  }
  // heavy mesons Kaons?
  static const G4double lMesCof1S = 0.82; // Kp/piP kaons?
  static const G4double llMesCof1C = 0.676568;
  static const G4double llMesCof1B = 0.610989;
  static const G4double llMesCof2C = 0.353135;
  static const G4double llMesCof2B = 0.221978;
  static const G4double llMesCofSC = 0.496568;
  static const G4double llMesCofSB = 0.430989;
  static const G4double llMesCofCB = 0.287557;
  static const G4double llMesCofEtaP = 0.88;
  static const G4double llMesCofEta = 0.76;

  if( pdg == 321 || pdg == 311 || pdg == 310 )
  {
    coeff = lMesCof1S; //K+-0
  }
  else if( pdg == 511 ||  pdg == 521  )
  {
    coeff = llMesCof1B; // BMeson0, BMeson+
  }
  else if(pdg == 421 ||  pdg == 411 )
  {
    coeff = llMesCof1C; // DMeson+, DMeson0
  }
  else if( pdg == 531  )
  {
    coeff = llMesCofSB; // BSMeson0
  }
  else if( pdg == 541 )
  {
    coeff = llMesCofCB; // BCMeson+-
  }
  else if(pdg == 431 ) 
  {
    coeff = llMesCofSC; // DSMeson+-
  }
  else if(pdg == 441 || pdg == 443 )
  {
    coeff = llMesCof2C; // Etac, JPsi
  }
  else if(pdg == 553 )
  {
    coeff = llMesCof2B; // Upsilon
  }
  else if(pdg == 221 )
  {
    coeff = llMesCofEta; // Eta
  }
  else if(pdg == 331 )
  {
    coeff = llMesCofEtaP; // Eta'
  } 
  return coeff;
}


