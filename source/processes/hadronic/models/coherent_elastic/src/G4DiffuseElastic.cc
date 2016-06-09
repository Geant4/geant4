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
// $Id: G4DiffuseElastic.cc,v 1.7 2007/06/12 14:46:26 grichine Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// Physics model class G4DiffuseElastic 
//
//
// G4 Model: optical diffuse elastic scattering with 4-momentum balance
//                         
// 24-May-07 V. Grichine
//

#include "G4DiffuseElastic.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4QElasticCrossSection.hh"
#include "G4VQCrossSection.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "Randomize.hh"
#include "G4Integrator.hh"
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

G4DiffuseElastic::G4DiffuseElastic() 
  : G4HadronicInteraction(), fParticle(0)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  verboseLevel = 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ  = 0.0*GeV;  
  lowEnergyLimitHE = 0.0*GeV;  
  lowestEnergyLimit= 0.0*keV;  
  plabLowLimit     = 20.0*MeV;

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();
  thePionPlus = G4PionPlus::PionPlus();
  thePionMinus= G4PionMinus::PionMinus();
}


G4DiffuseElastic::~G4DiffuseElastic()
{
}



G4HadFinalState* 
G4DiffuseElastic::ApplyYourself( const G4HadProjectile& aTrack, 
                                       G4Nucleus& targetNucleus )
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;

  G4double ekin = aParticle->GetKineticEnergy();

  if(ekin <= lowestEnergyLimit) 
  {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4double aTarget = targetNucleus.GetN();
  G4double zTarget = targetNucleus.GetZ();

  G4double plab = aParticle->GetTotalMomentum();

  if (verboseLevel >1)
  { 
    G4cout << "G4DiffuseElastic::DoIt: Incident particle plab=" 
	   << plab/GeV << " GeV/c " 
	   << " ekin(MeV) = " << ekin/MeV << "  " 
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;
  }
  // Scattered particle referred to axis of incident particle

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();

  G4int Z = static_cast<G4int>(zTarget+0.5);
  G4int A = static_cast<G4int>(aTarget+0.5);
  G4int N = A - Z;

  G4int projPDG = theParticle->GetPDGEncoding();

  if (verboseLevel>1) 
  {
    G4cout << "G4DiffuseElastic for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z 
	   << " A= " << A << " N= " << N 
	   << G4endl;
  }
  G4ParticleDefinition * theDef = 0;

  if(Z == 1 && A == 1)       theDef = theProton;
  else if (Z == 1 && A == 2) theDef = theDeuteron;
  else if (Z == 1 && A == 3) theDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theDef = G4He3::He3();
  else if (Z == 2 && A == 4) theDef = theAlpha;
  else theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);
 
  G4double m2 = theDef->GetPDGMass();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,m2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();
  G4double tmax = 4.0*ptot*ptot;
  G4double t = 0.0;


  //
  // Sample t
  //
  
  t = SampleT( theParticle, ptot, A);

  // NaN finder
  if(!(t < 0.0 || t >= 0.0)) 
  {
    if (verboseLevel > 0) 
    {
      G4cout << "G4DiffuseElastic:WARNING: Z= " << Z << " N= " 
	     << N << " pdg= " <<  projPDG
	     << " mom(GeV)= " << plab/GeV 
              << " S-wave will be sampled" 
	     << G4endl; 
    }
    t = G4UniformRand()*tmax; 
  }
  if(verboseLevel>1)
  {
    G4cout <<" t= " << t << " tmax= " << tmax 
	   << " ptot= " << ptot << G4endl;
  }
  // Sampling of angles in CM system

  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  if( cost >= 1.0 ) 
  {
    cost = 1.0;
    sint = 0.0;
  }
  else if( cost <= -1.0) 
  {
    cost = -1.0;
    sint =  0.0;
  }
  else  
  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }    
  if (verboseLevel>1) 
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;

  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;

  if (verboseLevel > 1)
  { 
    G4cout << "Scattered: "
	   << nlv1<<" m= " << m1 << " ekin(MeV)= " << eFinal 
	   << " Proj: 4-mom " << lv1 
	   <<G4endl;
  }
  if(eFinal < 0.0) 
  {
    G4cout << "G4DiffuseElastic WARNING ekin= " << eFinal
	   << " after scattering of " 
	   << aParticle->GetDefinition()->GetParticleName()
	   << " p(GeV/c)= " << plab
	   << " on " << theDef->GetParticleName()
	   << G4endl;
    eFinal = 0.0;
    nlv1.setE(m1);
  }

  theParticleChange.SetMomentumChange(nlv1.vect().unit());
  theParticleChange.SetEnergyChange(eFinal);
  
  G4LorentzVector nlv0 = lv - nlv1;
  G4double erec =  nlv0.e() - m2;

  if (verboseLevel > 1) 
  {
    G4cout << "Recoil: "
	   << nlv0<<" m= " << m2 << " ekin(MeV)= " << erec 
	   <<G4endl;
  }
  if(erec > lowEnergyRecoilLimit) 
  {
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, nlv0);
    theParticleChange.AddSecondary(aSec);
  } else {
    if(erec < 0.0) erec = 0.0;
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}


////////////////////////////////////////////////////////////////////////////
//
// return differential elastic cross section d(sigma)/d(omega) 

G4double 
G4DiffuseElastic::GetDiffuseElasticXsc( const G4ParticleDefinition* particle, 
                                        G4double theta, 
			                G4double momentum, 
                                        G4double A         )
{
  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;
  
  G4double r0;
  if(A > 10.) r0  = 1.16*( 1 - std::pow(A, -2./3.) )*fermi;   // 1.08*fermi;
  else        r0  = 1.1*fermi;
  fNuclearRadius = r0*std::pow(A, 1./3.);

  G4double sigma = fNuclearRadius*fNuclearRadius*GetDiffElasticProb(theta);

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability d(probability)/d(omega) 

G4double 
G4DiffuseElastic::GetDiffElasticProb( // G4ParticleDefinition* particle, 
                                        G4double theta 
					//  G4double momentum, 
					// G4double A         
                                     )
{
  G4double sigma, bzero, bzero2, bonebyarg, bonebyarg2, damp, damp2;
  G4double delta, diffuse, gamma;
  G4double e1, e2, bone, bone2;

  // G4double wavek = momentum/hbarc;  // wave vector
  // G4double r0    = 1.08*fermi;
  // G4double rad   = r0*std::pow(A, 1./3.);
  G4double kr    = fWaveVector*fNuclearRadius; // wavek*rad;
  G4double kr2   = kr*kr;
  G4double krt   = kr*theta;

  bzero      = BesselJzero(krt);
  bzero2     = bzero*bzero;    
  bone       = BesselJone(krt);
  bone2      = bone*bone;
  bonebyarg  = BesselOneByArg(krt);
  bonebyarg2 = bonebyarg*bonebyarg;  

  if (fParticle == theProton)
  {
    diffuse = 0.63*fermi;
    gamma   = 0.3*fermi;
    delta   = 0.1*fermi*fermi;
    e1      = 0.3*fermi;
    e2      = 0.35*fermi;
  }
  else // as proton, if were not defined 
  {
    diffuse = 0.63*fermi;
    gamma   = 0.3*fermi;
    delta   = 0.1*fermi*fermi;
    e1      = 0.3*fermi;
    e2      = 0.35*fermi;
  }
  G4double kg    = fWaveVector*gamma;   // wavek*delta;
  G4double kg2   = kg*kg;
  G4double dk2t  = delta*fWaveVector*fWaveVector*theta; // delta*wavek*wavek*theta;
  G4double dk2t2 = dk2t*dk2t;
  G4double pikdt = pi*fWaveVector*diffuse*theta;// pi*wavek*diffuse*theta;

  G4double mode2k2 = (e1*e1+e2*e2)*fWaveVector*fWaveVector;  
  G4double e2dk3t  = -2.*e2*delta*fWaveVector*fWaveVector*fWaveVector*theta;


  damp           = DampFactor(pikdt);
  damp2          = damp*damp;

  sigma  = kg2 + dk2t2;
  sigma *= bzero2;
  sigma += mode2k2*bone2 + e2dk3t*bzero*bone;
  sigma += kr2*bonebyarg2;
  sigma *= damp2;          // *rad*rad;

  return sigma;
}


////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability 2*pi*sin(theta)*d(probability)/d(omega) 

G4double 
G4DiffuseElastic::GetIntegrandFunction( G4double theta )
{
  G4double result;

  result  = 2*pi*std::sin(theta);
  result *= GetDiffElasticProb(theta);
  return result;
}

////////////////////////////////////////////////////////////////////////////
//
// return integral elastic cross section d(sigma)/d(omega) integrated 0 - theta 

G4double 
G4DiffuseElastic::IntegralElasticProb(  const G4ParticleDefinition* particle, 
                                        G4double theta, 
			                G4double momentum, 
                                        G4double A         )
{
  G4double result;
  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;
  G4double r0;
  if(A > 10.) r0  = 1.16*( 1 - std::pow(A, -2./3.) )*fermi;   // 1.08*fermi;
  else        r0  = 1.1*fermi;
  fNuclearRadius = r0*std::pow(A, 1./3.);


  G4Integrator<G4DiffuseElastic,G4double(G4DiffuseElastic::*)(G4double)> integral;

  // result = integral.Legendre10(this,&G4DiffuseElastic::GetIntegrandFunction, 0., theta ); 
  result = integral.Legendre96(this,&G4DiffuseElastic::GetIntegrandFunction, 0., theta ); 

  return result;
}

////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0

G4double G4DiffuseElastic::SampleT( const G4ParticleDefinition* aParticle, G4double p, G4double A)
{
  G4double theta = SampleThetaCMS( aParticle,  p, A); // sample theta in cms
  G4double t     = 2*p*p*( 1 - std::cos(theta) ); // -t !!!
  return t;
}

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle sampled in cms


G4double 
G4DiffuseElastic::SampleThetaCMS(const G4ParticleDefinition* particle, 
                                       G4double momentum, G4double A)
{
  G4int i, iMax = 100;  
  G4double r0, norm, result, theta1, theta2, thetaMax, sum = 0.;

  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;
  
  if(A > 10.) r0  = 1.16*( 1 - std::pow(A, -2./3.) )*fermi;   // 1.08*fermi;
  else        r0  = 1.1*fermi;

  fNuclearRadius  = r0*std::pow(A, 1./3.);
  
  thetaMax = 10.174/fWaveVector/fNuclearRadius;
  if (thetaMax > pi) thetaMax = pi;

  G4Integrator<G4DiffuseElastic,G4double(G4DiffuseElastic::*)(G4double)> integral;

  // result = integral.Legendre10(this,&G4DiffuseElastic::GetIntegrandFunction, 0., theta ); 
  norm = integral.Legendre96(this,&G4DiffuseElastic::GetIntegrandFunction, 0., thetaMax );

  norm *= G4UniformRand();

  for(i = 1; i <= iMax; i++)
  {
    theta1 = (i-1)*thetaMax/iMax; 
    theta2 = i*thetaMax/iMax;
    sum   += integral.Legendre10(this,&G4DiffuseElastic::GetIntegrandFunction, theta1, theta2);

    if ( sum >= norm ) 
    {
      result = 0.5*(theta1 + theta2);
      break;
    }
  }
  if (i > iMax ) result = 0.5*(theta1 + theta2);
  return result;
}


////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle sampled in lab system (target at rest)



G4double 
G4DiffuseElastic::SampleThetaLab( const G4HadProjectile* aParticle, 
                                        G4double tmass, G4double A)
{
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = aParticle->GetTotalMomentum();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,tmass);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();
  G4double tmax = 4.0*ptot*ptot;
  G4double t = 0.0;


  //
  // Sample t
  //
  
  t = SampleT( theParticle, ptot, A);

  // NaN finder
  if(!(t < 0.0 || t >= 0.0)) 
  {
    if (verboseLevel > 0) 
    {
      G4cout << "G4DiffuseElastic:WARNING: A = " << A 
	     << " mom(GeV)= " << plab/GeV 
             << " S-wave will be sampled" 
	     << G4endl; 
    }
    t = G4UniformRand()*tmax; 
  }
  if(verboseLevel>1)
  {
    G4cout <<" t= " << t << " tmax= " << tmax 
	   << " ptot= " << ptot << G4endl;
  }
  // Sampling of angles in CM system

  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  if( cost >= 1.0 ) 
  {
    cost = 1.0;
    sint = 0.0;
  }
  else if( cost <= -1.0) 
  {
    cost = -1.0;
    sint =  0.0;
  }
  else  
  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }    
  if (verboseLevel>1) 
  {
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4ThreeVector np1 = nlv1.vect();

    // G4double theta = std::acos( np1.z()/np1.mag() );  // degree;

  G4double theta = np1.theta();

  return theta;
}
