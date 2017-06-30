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
// $Id: G4NuclNuclDiffuseElastic.cc 104887 2017-06-26 07:12:43Z gcosmo $
//
//
// Physics model class G4NuclNuclDiffuseElastic 
//
//
// G4 Model: optical diffuse elastic scattering with 4-momentum balance
//                         
// 24-May-07 V. Grichine
//

#include "G4NuclNuclDiffuseElastic.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"

#include "Randomize.hh"
#include "G4Integrator.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"

/////////////////////////////////////////////////////////////////////////
//
// Test Constructor. Just to check xsc


G4NuclNuclDiffuseElastic::G4NuclNuclDiffuseElastic() 
  : G4HadronElastic("NNDiffuseElastic"), fParticle(0)
{
  SetMinEnergy( 50*MeV );
  SetMaxEnergy( 1.*TeV );
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

  fEnergyBin = 200;
  fAngleBin  = 200;

  fEnergyVector =  new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );
  fAngleTable = 0;

  fParticle = 0;
  fWaveVector = 0.;
  fAtomicWeight = 0.;
  fAtomicNumber = 0.;
  fNuclearRadius = 0.;
  fBeta = 0.;
  fZommerfeld = 0.;
  fAm = 0.;
  fAddCoulomb = false;
  // Ranges of angle table relative to current Rutherford (Coulomb grazing) angle

  // Empirical parameters

  fCofAlphaMax     = 1.5;
  fCofAlphaCoulomb = 0.5;

  fProfileDelta  = 1.;
  fProfileAlpha  = 0.5;

  fCofLambda = 1.0;
  fCofDelta  = 0.04;
  fCofAlpha  = 0.095;

  fNuclearRadius1 = fNuclearRadius2 = fNuclearRadiusSquare 
    = fRutherfordRatio = fCoulombPhase0 = fHalfRutThetaTg = fHalfRutThetaTg2 
    = fRutherfordTheta = fProfileLambda = fCofPhase = fCofFar
    = fSumSigma = fEtaRatio = fReZ = 0.0;
  fMaxL = 0;

  fNuclearRadiusCof = 1.0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4NuclNuclDiffuseElastic::~G4NuclNuclDiffuseElastic()
{
  if ( fEnergyVector ) {
    delete fEnergyVector;
    fEnergyVector = 0;
  }

  for ( std::vector<G4PhysicsTable*>::iterator it = fAngleBank.begin();
        it != fAngleBank.end(); ++it ) {
    if ( (*it) ) (*it)->clearAndDestroy();
    delete *it;
    *it = 0;
  }
  fAngleTable = 0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Initialisation for given particle using element table of application

void G4NuclNuclDiffuseElastic::Initialise() 
{

  // fEnergyVector = new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );

  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  size_t jEl, numOfEl = G4Element::GetNumberOfElements();

  // projectile radius

  G4double A1 = G4double( fParticle->GetBaryonNumber() );
  G4double R1 = CalculateNuclearRad(A1);

  for(jEl = 0 ; jEl < numOfEl; ++jEl) // application element loop
  {
    fAtomicNumber = (*theElementTable)[jEl]->GetZ();     // atomic number
    fAtomicWeight = G4NistManager::Instance()->GetAtomicMassAmu( static_cast< G4int >( fAtomicNumber ) );

    fNuclearRadius = CalculateNuclearRad(fAtomicWeight);
    fNuclearRadius += R1;

    if(verboseLevel > 0) 
    {   
      G4cout<<"G4NuclNuclDiffuseElastic::Initialise() the element: "
	    <<(*theElementTable)[jEl]->GetName()<<G4endl;
    }
    fElementNumberVector.push_back(fAtomicNumber);
    fElementNameVector.push_back((*theElementTable)[jEl]->GetName());

    BuildAngleTable();
    fAngleBank.push_back(fAngleTable);
  }  
}


////////////////////////////////////////////////////////////////////////////
//
// return differential elastic cross section d(sigma)/d(omega) 

G4double 
G4NuclNuclDiffuseElastic::GetDiffuseElasticXsc( const G4ParticleDefinition* particle, 
                                        G4double theta, 
			                G4double momentum, 
                                        G4double A         )
{
  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;
  fAddCoulomb    = false;
  fNuclearRadius = CalculateNuclearRad(A);

  G4double sigma = fNuclearRadius*fNuclearRadius*GetDiffElasticProb(theta);

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return invariant differential elastic cross section d(sigma)/d(tMand) 

G4double 
G4NuclNuclDiffuseElastic::GetInvElasticXsc( const G4ParticleDefinition* particle, 
                                        G4double tMand, 
			                G4double plab, 
                                        G4double A, G4double Z         )
{
  G4double m1 = particle->GetPDGMass();
  G4LorentzVector lv1(0.,0.,plab,std::sqrt(plab*plab+m1*m1));

  G4int iZ = static_cast<G4int>(Z+0.5);
  G4int iA = static_cast<G4int>(A+0.5);
  G4ParticleDefinition * theDef = 0;

  if      (iZ == 1 && iA == 1) theDef = theProton;
  else if (iZ == 1 && iA == 2) theDef = theDeuteron;
  else if (iZ == 1 && iA == 3) theDef = G4Triton::Triton();
  else if (iZ == 2 && iA == 3) theDef = G4He3::He3();
  else if (iZ == 2 && iA == 4) theDef = theAlpha;
  else theDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(iZ,iA,0);
 
  G4double tmass = theDef->GetPDGMass();

  G4LorentzVector lv(0.0,0.0,0.0,tmass);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot    = p1.mag();
  G4double ptot2 = ptot*ptot;
  G4double cost = 1 - 0.5*std::fabs(tMand)/ptot2;

  if( cost >= 1.0 )      cost = 1.0;  
  else if( cost <= -1.0) cost = -1.0;
  
  G4double thetaCMS = std::acos(cost);

  G4double sigma = GetDiffuseElasticXsc( particle, thetaCMS, ptot, A);

  sigma *= pi/ptot2;

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return differential elastic cross section d(sigma)/d(omega) with Coulomb
// correction

G4double 
G4NuclNuclDiffuseElastic::GetDiffuseElasticSumXsc( const G4ParticleDefinition* particle, 
                                        G4double theta, 
			                G4double momentum, 
                                        G4double A, G4double Z         )
{
  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;
  fAtomicNumber  = Z;
  fNuclearRadius = CalculateNuclearRad(A);
  fAddCoulomb    = false;

  G4double z     = particle->GetPDGCharge();

  G4double kRt   = fWaveVector*fNuclearRadius*theta;
  G4double kRtC  = 1.9;

  if( z && (kRt > kRtC) )
  {
    fAddCoulomb = true;
    fBeta       = CalculateParticleBeta( particle, momentum);
    fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
    fAm         = CalculateAm( momentum, fZommerfeld, fAtomicNumber);
  }
  G4double sigma = fNuclearRadius*fNuclearRadius*GetDiffElasticSumProb(theta);

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return invariant differential elastic cross section d(sigma)/d(tMand) with Coulomb
// correction

G4double 
G4NuclNuclDiffuseElastic::GetInvElasticSumXsc( const G4ParticleDefinition* particle, 
                                        G4double tMand, 
			                G4double plab, 
                                        G4double A, G4double Z         )
{
  G4double m1 = particle->GetPDGMass();

  G4LorentzVector lv1(0.,0.,plab,std::sqrt(plab*plab+m1*m1));

  G4int iZ = static_cast<G4int>(Z+0.5);
  G4int iA = static_cast<G4int>(A+0.5);

  G4ParticleDefinition* theDef = 0;

  if      (iZ == 1 && iA == 1) theDef = theProton;
  else if (iZ == 1 && iA == 2) theDef = theDeuteron;
  else if (iZ == 1 && iA == 3) theDef = G4Triton::Triton();
  else if (iZ == 2 && iA == 3) theDef = G4He3::He3();
  else if (iZ == 2 && iA == 4) theDef = theAlpha;
  else theDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(iZ,iA,0);
 
  G4double tmass = theDef->GetPDGMass();

  G4LorentzVector lv(0.0,0.0,0.0,tmass);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot    = p1.mag();
  G4double ptot2   = ptot*ptot;
  G4double cost    = 1 - 0.5*std::fabs(tMand)/ptot2;

  if( cost >= 1.0 )      cost = 1.0;  
  else if( cost <= -1.0) cost = -1.0;
  
  G4double thetaCMS = std::acos(cost);

  G4double sigma = GetDiffuseElasticSumXsc( particle, thetaCMS, ptot, A, Z );

  sigma *= pi/ptot2;

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return invariant differential elastic cross section d(sigma)/d(tMand) with Coulomb
// correction

G4double 
G4NuclNuclDiffuseElastic::GetInvCoulombElasticXsc( const G4ParticleDefinition* particle, 
                                        G4double tMand, 
			                G4double plab, 
                                        G4double A, G4double Z         )
{
  G4double m1 = particle->GetPDGMass();
  G4LorentzVector lv1(0.,0.,plab,std::sqrt(plab*plab+m1*m1));

  G4int iZ = static_cast<G4int>(Z+0.5);
  G4int iA = static_cast<G4int>(A+0.5);
  G4ParticleDefinition * theDef = 0;

  if      (iZ == 1 && iA == 1) theDef = theProton;
  else if (iZ == 1 && iA == 2) theDef = theDeuteron;
  else if (iZ == 1 && iA == 3) theDef = G4Triton::Triton();
  else if (iZ == 2 && iA == 3) theDef = G4He3::He3();
  else if (iZ == 2 && iA == 4) theDef = theAlpha;
  else theDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(iZ,iA,0);
 
  G4double tmass = theDef->GetPDGMass();

  G4LorentzVector lv(0.0,0.0,0.0,tmass);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot    = p1.mag();
  G4double ptot2 = ptot*ptot;
  G4double cost = 1 - 0.5*std::fabs(tMand)/ptot2;

  if( cost >= 1.0 )      cost = 1.0;  
  else if( cost <= -1.0) cost = -1.0;
  
  G4double thetaCMS = std::acos(cost);

  G4double sigma = GetCoulombElasticXsc( particle, thetaCMS, ptot, Z );

  sigma *= pi/ptot2;

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability d(probability)/d(omega) 

G4double 
G4NuclNuclDiffuseElastic::GetDiffElasticProb( // G4ParticleDefinition* particle, 
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
  // G4double rad   = r0*G4Pow::GetInstance()->A13(A);

  G4double kr    = fWaveVector*fNuclearRadius; // wavek*rad;
  G4double kr2   = kr*kr;
  G4double krt   = kr*theta;

  bzero      = BesselJzero(krt);
  bzero2     = bzero*bzero;    
  bone       = BesselJone(krt);
  bone2      = bone*bone;
  bonebyarg  = BesselOneByArg(krt);
  bonebyarg2 = bonebyarg*bonebyarg;  

  // VI - Coverity complains
  /*
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
  */
    diffuse = 0.63*fermi;
    gamma   = 0.3*fermi;
    delta   = 0.1*fermi*fermi;
    e1      = 0.3*fermi;
    e2      = 0.35*fermi;
    //}
  G4double lambda = 15.; // 15 ok

  //  G4double kgamma    = fWaveVector*gamma;   // wavek*delta;

  G4double kgamma    = lambda*(1.-G4Exp(-fWaveVector*gamma/lambda));   // wavek*delta;
  G4double kgamma2   = kgamma*kgamma;

  // G4double dk2t  = delta*fWaveVector*fWaveVector*theta; // delta*wavek*wavek*theta;
  // G4double dk2t2 = dk2t*dk2t;
  // G4double pikdt = pi*fWaveVector*diffuse*theta;// pi*wavek*diffuse*theta;

  G4double pikdt    = lambda*(1.-G4Exp(-pi*fWaveVector*diffuse*theta/lambda));   // wavek*delta;

  damp           = DampFactor(pikdt);
  damp2          = damp*damp;

  G4double mode2k2 = (e1*e1+e2*e2)*fWaveVector*fWaveVector;  
  G4double e2dk3t  = -2.*e2*delta*fWaveVector*fWaveVector*fWaveVector*theta;


  sigma  = kgamma2;
  // sigma  += dk2t2;
  sigma *= bzero2;
  sigma += mode2k2*bone2 + e2dk3t*bzero*bone;
  sigma += kr2*bonebyarg2;
  sigma *= damp2;          // *rad*rad;

  return sigma;
}

////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability d(probability)/d(omega) with 
// Coulomb correction

G4double 
G4NuclNuclDiffuseElastic::GetDiffElasticSumProb( // G4ParticleDefinition* particle, 
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
  // G4double rad   = r0*G4Pow::GetInstance()->A13(A);

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
    // diffuse = 0.6*fermi;
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
  G4double lambda = 15.; // 15 ok
  // G4double kgamma    = fWaveVector*gamma;   // wavek*delta;
  G4double kgamma    = lambda*(1.-G4Exp(-fWaveVector*gamma/lambda));   // wavek*delta;

  // G4cout<<"kgamma = "<<kgamma<<G4endl;

  if(fAddCoulomb)  // add Coulomb correction
  {
    G4double sinHalfTheta  = std::sin(0.5*theta);
    G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta;

    kgamma += 0.5*fZommerfeld/kr/(sinHalfTheta2+fAm); // correction at J0()
  // kgamma += 0.65*fZommerfeld/kr/(sinHalfTheta2+fAm); // correction at J0()
  }

  G4double kgamma2   = kgamma*kgamma;

  // G4double dk2t  = delta*fWaveVector*fWaveVector*theta; // delta*wavek*wavek*theta;
  //   G4cout<<"dk2t = "<<dk2t<<G4endl;
  // G4double dk2t2 = dk2t*dk2t;
  // G4double pikdt = pi*fWaveVector*diffuse*theta;// pi*wavek*diffuse*theta;

  G4double pikdt    = lambda*(1.-G4Exp(-pi*fWaveVector*diffuse*theta/lambda));   // wavek*delta;

  // G4cout<<"pikdt = "<<pikdt<<G4endl;

  damp           = DampFactor(pikdt);
  damp2          = damp*damp;

  G4double mode2k2 = (e1*e1+e2*e2)*fWaveVector*fWaveVector;  
  G4double e2dk3t  = -2.*e2*delta*fWaveVector*fWaveVector*fWaveVector*theta;

  sigma  = kgamma2;
  // sigma += dk2t2;
  sigma *= bzero2;
  sigma += mode2k2*bone2; 
  sigma += e2dk3t*bzero*bone;

  // sigma += kr2*(1 + 8.*fZommerfeld*fZommerfeld/kr2)*bonebyarg2;  // correction at J1()/()
  sigma += kr2*bonebyarg2;  // correction at J1()/()

  sigma *= damp2;          // *rad*rad;

  return sigma;
}


////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability d(probability)/d(t) with 
// Coulomb correction

G4double 
G4NuclNuclDiffuseElastic::GetDiffElasticSumProbA( G4double alpha )
{
  G4double theta; 

  theta = std::sqrt(alpha);

  // theta = std::acos( 1 - alpha/2. );

  G4double sigma, bzero, bzero2, bonebyarg, bonebyarg2, damp, damp2;
  G4double delta, diffuse, gamma;
  G4double e1, e2, bone, bone2;

  // G4double wavek = momentum/hbarc;  // wave vector
  // G4double r0    = 1.08*fermi;
  // G4double rad   = r0*G4Pow::GetInstance()->A13(A);

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
    // diffuse = 0.6*fermi;
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
  G4double lambda = 15.; // 15 ok
  // G4double kgamma    = fWaveVector*gamma;   // wavek*delta;
  G4double kgamma    = lambda*(1.-G4Exp(-fWaveVector*gamma/lambda));   // wavek*delta;

  // G4cout<<"kgamma = "<<kgamma<<G4endl;

  if(fAddCoulomb)  // add Coulomb correction
  {
    G4double sinHalfTheta  = theta*0.5; // std::sin(0.5*theta);
    G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta;

    kgamma += 0.5*fZommerfeld/kr/(sinHalfTheta2+fAm); // correction at J0()
  // kgamma += 0.65*fZommerfeld/kr/(sinHalfTheta2+fAm); // correction at J0()
  }

  G4double kgamma2   = kgamma*kgamma;

  // G4double dk2t  = delta*fWaveVector*fWaveVector*theta; // delta*wavek*wavek*theta;
  //   G4cout<<"dk2t = "<<dk2t<<G4endl;
  // G4double dk2t2 = dk2t*dk2t;
  // G4double pikdt = pi*fWaveVector*diffuse*theta;// pi*wavek*diffuse*theta;

  G4double pikdt    = lambda*(1.-G4Exp(-pi*fWaveVector*diffuse*theta/lambda));   // wavek*delta;

  // G4cout<<"pikdt = "<<pikdt<<G4endl;

  damp           = DampFactor(pikdt);
  damp2          = damp*damp;

  G4double mode2k2 = (e1*e1+e2*e2)*fWaveVector*fWaveVector;  
  G4double e2dk3t  = -2.*e2*delta*fWaveVector*fWaveVector*fWaveVector*theta;

  sigma  = kgamma2;
  // sigma += dk2t2;
  sigma *= bzero2;
  sigma += mode2k2*bone2; 
  sigma += e2dk3t*bzero*bone;

  // sigma += kr2*(1 + 8.*fZommerfeld*fZommerfeld/kr2)*bonebyarg2;  // correction at J1()/()
  sigma += kr2*bonebyarg2;  // correction at J1()/()

  sigma *= damp2;          // *rad*rad;

  return sigma;
}


////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability 2*pi*sin(theta)*d(probability)/d(omega) 

G4double 
G4NuclNuclDiffuseElastic::GetIntegrandFunction( G4double alpha )
{
  G4double result;

  result  = GetDiffElasticSumProbA(alpha);

  // result *= 2*pi*std::sin(theta);

  return result;
}

////////////////////////////////////////////////////////////////////////////
//
// return integral elastic cross section d(sigma)/d(omega) integrated 0 - theta 

G4double 
G4NuclNuclDiffuseElastic::IntegralElasticProb(  const G4ParticleDefinition* particle, 
                                        G4double theta, 
			                G4double momentum, 
                                        G4double A         )
{
  G4double result;
  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;

  fNuclearRadius = CalculateNuclearRad(A);


  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  // result = integral.Legendre10(this,&G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., theta ); 
  result = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., theta ); 

  return result;
}

////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0

G4double G4NuclNuclDiffuseElastic::SampleT( const G4ParticleDefinition* aParticle, 
					    G4double p, G4double A)
{
  G4double theta = SampleThetaCMS( aParticle,  p, A); // sample theta in cms
  G4double t     = 2*p*p*( 1 - std::cos(theta) ); // -t !!!
  return t;
}

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle sampled in cms


G4double 
G4NuclNuclDiffuseElastic::SampleThetaCMS(const G4ParticleDefinition* particle, 
                                       G4double momentum, G4double A)
{
  G4int i, iMax = 100;  
  G4double norm, result, theta1, theta2, thetaMax, sum = 0.;

  fParticle      = particle;
  fWaveVector    = momentum/hbarc;
  fAtomicWeight  = A;

  fNuclearRadius = CalculateNuclearRad(A);
  
  thetaMax = 10.174/fWaveVector/fNuclearRadius;

  if (thetaMax > pi) thetaMax = pi;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  // result = integral.Legendre10(this,&G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., theta ); 
  norm = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., thetaMax );

  norm *= G4UniformRand();

  for(i = 1; i <= iMax; i++)
  {
    theta1 = (i-1)*thetaMax/iMax; 
    theta2 = i*thetaMax/iMax;
    sum   += integral.Legendre10(this,&G4NuclNuclDiffuseElastic::GetIntegrandFunction, theta1, theta2);

    if ( sum >= norm ) 
    {
      result = 0.5*(theta1 + theta2);
      break;
    }
  }
  if (i > iMax ) result = 0.5*(theta1 + theta2);

  G4double sigma = pi*thetaMax/iMax;

  result += G4RandGauss::shoot(0.,sigma);

  if(result < 0.) result = 0.;
  if(result > thetaMax) result = thetaMax;

  return result;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////  Table preparation and reading ////////////////////////

////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4NuclNuclDiffuseElastic::SampleInvariantT( const G4ParticleDefinition* aParticle, G4double p, 
                                               G4int Z, G4int A)
{
  fParticle = aParticle;
  G4double m1 = fParticle->GetPDGMass();
  G4double totElab = std::sqrt(m1*m1+p*p);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv1(p,0.0,0.0,totElab);
  G4LorentzVector  lv(0.0,0.0,0.0,mass2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double momentumCMS = p1.mag();

  G4double t = SampleTableT( aParticle,  momentumCMS, G4double(Z), G4double(A) ); // sample theta2 in cms
  return t;
}

////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4NuclNuclDiffuseElastic::SampleTableT( const G4ParticleDefinition* aParticle, G4double p, 
                                               G4double Z, G4double A)
{
  G4double alpha = SampleTableThetaCMS( aParticle,  p, Z, A); // sample theta2 in cms
  // G4double t     = 2*p*p*( 1 - std::cos(std::sqrt(alpha)) );             // -t !!!
  G4double t     = p*p*alpha;             // -t !!!
  return t;
}

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle2 sampled in cms according to precalculated table.


G4double 
G4NuclNuclDiffuseElastic::SampleTableThetaCMS(const G4ParticleDefinition* particle, 
                                       G4double momentum, G4double Z, G4double A)
{
  size_t iElement;
  G4int iMomentum, iAngle;  
  G4double randAngle, position, theta1, theta2, E1, E2, W1, W2, W;  
  G4double m1 = particle->GetPDGMass();

  for(iElement = 0; iElement < fElementNumberVector.size(); iElement++)
  {
    if( std::fabs(Z - fElementNumberVector[iElement]) < 0.5) break;
  }
  if ( iElement == fElementNumberVector.size() ) 
  {
    InitialiseOnFly(Z,A); // table preparation, if needed

    // iElement--;

    // G4cout << "G4NuclNuclDiffuseElastic: Element with atomic number " << Z
    // << " is not found, return zero angle" << G4endl;
    // return 0.; // no table for this element
  }
  // G4cout<<"iElement = "<<iElement<<G4endl;

  fAngleTable = fAngleBank[iElement];

  G4double kinE = std::sqrt(momentum*momentum + m1*m1) - m1;

  for( iMomentum = 0; iMomentum < fEnergyBin; iMomentum++)
  {
    // G4cout<<iMomentum<<", kinE = "<<kinE/MeV<<", vectorE = "<<fEnergyVector->GetLowEdgeEnergy(iMomentum)/MeV<<G4endl;
    
    if( kinE < fEnergyVector->GetLowEdgeEnergy(iMomentum) ) break;
  }
  // G4cout<<"iMomentum = "<<iMomentum<<", fEnergyBin -1 = "<<fEnergyBin -1<<G4endl;


  if ( iMomentum >= fEnergyBin ) iMomentum = fEnergyBin-1;   // kinE is more then theMaxEnergy
  if ( iMomentum < 0 )           iMomentum = 0; // against negative index, kinE < theMinEnergy


  if (iMomentum == fEnergyBin -1 || iMomentum == 0 )   // the table edges
  {
    position = (*(*fAngleTable)(iMomentum))(fAngleBin-2)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iAngle = 0; iAngle < fAngleBin-1; iAngle++)
    {
      if( position < (*(*fAngleTable)(iMomentum))(iAngle) ) break;
    }
    if (iAngle >= fAngleBin-1) iAngle = fAngleBin-2;

    // G4cout<<"iAngle = "<<iAngle<<G4endl;

    randAngle = GetScatteringAngle(iMomentum, iAngle, position);

    // G4cout<<"randAngle = "<<randAngle<<G4endl;
  }
  else  // kinE inside between energy table edges
  {
    // position = (*(*fAngleTable)(iMomentum))(fAngleBin-2)*G4UniformRand();
    position = (*(*fAngleTable)(iMomentum))(0)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iAngle = 0; iAngle < fAngleBin-1; iAngle++)
    {
      // if( position < (*(*fAngleTable)(iMomentum))(iAngle) ) break;
      if( position > (*(*fAngleTable)(iMomentum))(iAngle) ) break;
    }
    if (iAngle >= fAngleBin-1) iAngle = fAngleBin-2;

    // G4cout<<"iAngle = "<<iAngle<<G4endl;

    theta2  = GetScatteringAngle(iMomentum, iAngle, position);

    // G4cout<<"theta2 = "<<theta2<<G4endl;

    E2 = fEnergyVector->GetLowEdgeEnergy(iMomentum);

    // G4cout<<"E2 = "<<E2<<G4endl;
    
    iMomentum--;
    
    // position = (*(*fAngleTable)(iMomentum))(fAngleBin-2)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iAngle = 0; iAngle < fAngleBin-1; iAngle++)
    {
      // if( position < (*(*fAngleTable)(iMomentum))(iAngle) ) break;
      if( position > (*(*fAngleTable)(iMomentum))(iAngle) ) break;
    }
    if (iAngle >= fAngleBin-1) iAngle = fAngleBin-2;
    
    theta1  = GetScatteringAngle(iMomentum, iAngle, position);

    // G4cout<<"theta1 = "<<theta1<<G4endl;

    E1 = fEnergyVector->GetLowEdgeEnergy(iMomentum);

    // G4cout<<"E1 = "<<E1<<G4endl;

    W  = 1.0/(E2 - E1);
    W1 = (E2 - kinE)*W;
    W2 = (kinE - E1)*W;

    randAngle = W1*theta1 + W2*theta2;
    
    // randAngle = theta2;
  }
  // G4double angle = randAngle;
  // if (randAngle > 0.) randAngle /= 2*pi*std::sin(angle);
  // G4cout<<"randAngle = "<<randAngle/degree<<G4endl;

  return randAngle;
}

//////////////////////////////////////////////////////////////////////////////
//
// Initialisation for given particle on fly using new element number

void G4NuclNuclDiffuseElastic::InitialiseOnFly(G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // atomic number
  fAtomicWeight  = G4NistManager::Instance()->GetAtomicMassAmu( static_cast< G4int >( Z ) );

  G4double A1 = G4double( fParticle->GetBaryonNumber() );
  G4double R1 = CalculateNuclearRad(A1);

  fNuclearRadius = CalculateNuclearRad(fAtomicWeight) + R1;
  
  if( verboseLevel > 0 )    
  {
    G4cout<<"G4NuclNuclDiffuseElastic::Initialise() the element with Z = "
	  <<Z<<"; and A = "<<A<<G4endl;
  }
  fElementNumberVector.push_back(fAtomicNumber);

  BuildAngleTable();

  fAngleBank.push_back(fAngleTable);

  return;
}

///////////////////////////////////////////////////////////////////////////////
//
// Build for given particle and element table of momentum, angle probability.
// For the moment in lab system. 

void G4NuclNuclDiffuseElastic::BuildAngleTable() 
{
  G4int i, j;
  G4double partMom, kinE, m1 = fParticle->GetPDGMass();
  G4double alpha1, alpha2, alphaMax, alphaCoulomb, delta = 0., sum = 0.;

  // G4cout<<"particle z = "<<z<<"; particle m1 = "<<m1/GeV<<" GeV"<<G4endl;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;
  
  fAngleTable = new G4PhysicsTable(fEnergyBin);

  for( i = 0; i < fEnergyBin; i++)
  {
    kinE        = fEnergyVector->GetLowEdgeEnergy(i);

    // G4cout<<G4endl;
    // G4cout<<"kinE = "<<kinE/MeV<<" MeV"<<G4endl;

    partMom     = std::sqrt( kinE*(kinE + 2*m1) );

    InitDynParameters(fParticle, partMom);

    alphaMax = fRutherfordTheta*fCofAlphaMax;

    if(alphaMax > pi) alphaMax = pi;

    // VI: Coverity complain
    //alphaMax = pi2;

    alphaCoulomb = fRutherfordTheta*fCofAlphaCoulomb;

    // G4cout<<"alphaCoulomb = "<<alphaCoulomb/degree<<"; alphaMax = "<<alphaMax/degree<<G4endl;

    G4PhysicsFreeVector* angleVector = new G4PhysicsFreeVector(fAngleBin-1);

    // G4PhysicsLogVector*  angleBins = new G4PhysicsLogVector( 0.001*alphaMax, alphaMax, fAngleBin );

    G4double delth = (alphaMax-alphaCoulomb)/fAngleBin;
        
    sum = 0.;

    // fAddCoulomb = false;
    fAddCoulomb = true;

    // for(j = 1; j < fAngleBin; j++)
    for(j = fAngleBin-1; j >= 1; j--)
    {
      // alpha1 = angleBins->GetLowEdgeEnergy(j-1);
      // alpha2 = angleBins->GetLowEdgeEnergy(j);

      // alpha1 = alphaMax*(j-1)/fAngleBin;
      // alpha2 = alphaMax*( j )/fAngleBin;

      alpha1 = alphaCoulomb + delth*(j-1);
      // if(alpha1 < kRlim2) alpha1 = kRlim2;
      alpha2 = alpha1 + delth;

      delta = integral.Legendre10(this, &G4NuclNuclDiffuseElastic::GetFresnelIntegrandXsc, alpha1, alpha2);
      // delta = integral.Legendre96(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, alpha1, alpha2);

      sum += delta;
      
      angleVector->PutValue( j-1 , alpha1, sum ); // alpha2
      // G4cout<<"j-1 = "<<j-1<<"; alpha2 = "<<alpha2/degree<<"; sum = "<<sum<<G4endl;
    }
    fAngleTable->insertAt(i,angleVector);

    // delete[] angleVector; 
    // delete[] angleBins; 
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4NuclNuclDiffuseElastic:: GetScatteringAngle( G4int iMomentum, G4int iAngle, G4double position )
{
 G4double x1, x2, y1, y2, randAngle;

  if( iAngle == 0 )
  {
    randAngle = (*fAngleTable)(iMomentum)->GetLowEdgeEnergy(iAngle);
    // iAngle++;
  }
  else
  {
    if ( iAngle >= G4int((*fAngleTable)(iMomentum)->GetVectorLength()) )
    {
      iAngle = (*fAngleTable)(iMomentum)->GetVectorLength() - 1;
    }
    y1 = (*(*fAngleTable)(iMomentum))(iAngle-1);
    y2 = (*(*fAngleTable)(iMomentum))(iAngle);

    x1 = (*fAngleTable)(iMomentum)->GetLowEdgeEnergy(iAngle-1);
    x2 = (*fAngleTable)(iMomentum)->GetLowEdgeEnergy(iAngle);

    if ( x1 == x2 )   randAngle = x2;
    else
    {
      if ( y1 == y2 ) randAngle = x1 + ( x2 - x1 )*G4UniformRand();
      else
      {
        randAngle = x1 + ( position - y1 )*( x2 - x1 )/( y2 - y1 );
      }
    }
  }
  return randAngle;
}



////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle sampled in lab system (target at rest)



G4double 
G4NuclNuclDiffuseElastic::SampleThetaLab( const G4HadProjectile* aParticle, 
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
  G4double ptot    = p1.mag();
  G4double tmax    = 4.0*ptot*ptot;
  G4double t       = 0.0;


  //
  // Sample t
  //
  
  t = SampleT( theParticle, ptot, A);

  // NaN finder
  if(!(t < 0.0 || t >= 0.0)) 
  {
    if (verboseLevel > 0) 
    {
      G4cout << "G4NuclNuclDiffuseElastic:WARNING: A = " << A 
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

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle in lab system (target at rest) knowing theta in CMS



G4double 
G4NuclNuclDiffuseElastic::ThetaCMStoThetaLab( const G4DynamicParticle* aParticle, 
                                        G4double tmass, G4double thetaCMS)
{
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  // G4double plab = aParticle->GetTotalMomentum();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,tmass);   

  lv += lv1;

  G4ThreeVector bst = lv.boostVector();

  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot    = p1.mag();

  G4double phi  = G4UniformRand()*twopi;
  G4double cost = std::cos(thetaCMS);
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
    G4cout << "cos(tcms)=" << cost << " std::sin(tcms)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4ThreeVector np1 = nlv1.vect();


  G4double thetaLab = np1.theta();

  return thetaLab;
}

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle in CMS system (target at rest) knowing theta in Lab



G4double 
G4NuclNuclDiffuseElastic::ThetaLabToThetaCMS( const G4DynamicParticle* aParticle, 
                                        G4double tmass, G4double thetaLab)
{
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = aParticle->GetTotalMomentum();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,tmass);   

  lv += lv1;

  G4ThreeVector bst = lv.boostVector();

  // lv1.boost(-bst);

  // G4ThreeVector p1 = lv1.vect();
  // G4double ptot    = p1.mag();

  G4double phi  = G4UniformRand()*twopi;
  G4double cost = std::cos(thetaLab);
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
    G4cout << "cos(tlab)=" << cost << " std::sin(tlab)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= plab;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(plab*plab + m1*m1));

  nlv1.boost(-bst); 

  G4ThreeVector np1 = nlv1.vect();


  G4double thetaCMS = np1.theta();

  return thetaCMS;
}

///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the moment in lab system. 

void G4NuclNuclDiffuseElastic::TestAngleTable(const G4ParticleDefinition* theParticle, G4double partMom,
                                            G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // atomic number
  fAtomicWeight  = A;     // number of nucleons
  fNuclearRadius = CalculateNuclearRad(fAtomicWeight);
  
     
  
  G4cout<<"G4NuclNuclDiffuseElastic::TestAngleTable() init the element with Z = "
	  <<Z<<"; and A = "<<A<<G4endl;
 
  fElementNumberVector.push_back(fAtomicNumber);

 


  G4int i=0, j;
  G4double a = 0., z = theParticle->GetPDGCharge(),  m1 = fParticle->GetPDGMass();
  G4double alpha1=0., alpha2=0., alphaMax=0., alphaCoulomb=0.;
  G4double deltaL10 = 0., deltaL96 = 0., deltaAG = 0.;
  G4double sumL10 = 0.,sumL96 = 0.,sumAG = 0.;
  G4double epsilon = 0.001;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;
  
  fAngleTable = new G4PhysicsTable(fEnergyBin);

  fWaveVector = partMom/hbarc;

  G4double kR     = fWaveVector*fNuclearRadius;
  G4double kR2    = kR*kR;
  G4double kRmax  = 10.6; // 10.6, 18, 10.174; ~ 3 maxima of J1 or 15., 25.
  G4double kRcoul = 1.2; // 1.4, 2.5; // on the first slope of J1

  alphaMax = kRmax*kRmax/kR2;

  if (alphaMax > 4.) alphaMax = 4.;  // vmg05-02-09: was pi2 

  alphaCoulomb = kRcoul*kRcoul/kR2;

  if( z )
  {
      a           = partMom/m1; // beta*gamma for m1
      fBeta       = a/std::sqrt(1+a*a);
      fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
      fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }
  G4PhysicsFreeVector* angleVector = new G4PhysicsFreeVector(fAngleBin-1);

  // G4PhysicsLogVector*  angleBins = new G4PhysicsLogVector( 0.001*alphaMax, alphaMax, fAngleBin );
    
  
  fAddCoulomb = false;

  for(j = 1; j < fAngleBin; j++)
  {
      // alpha1 = angleBins->GetLowEdgeEnergy(j-1);
      // alpha2 = angleBins->GetLowEdgeEnergy(j);

    alpha1 = alphaMax*(j-1)/fAngleBin;
    alpha2 = alphaMax*( j )/fAngleBin;

    if( ( alpha2 > alphaCoulomb ) && z ) fAddCoulomb = true;

    deltaL10 = integral.Legendre10(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, alpha1, alpha2);
    deltaL96 = integral.Legendre96(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, alpha1, alpha2);
    deltaAG  = integral.AdaptiveGauss(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, 
                                       alpha1, alpha2,epsilon);

      // G4cout<<alpha1<<"\t"<<std::sqrt(alpha1)/degree<<"\t"
      //     <<deltaL10<<"\t"<<deltaL96<<"\t"<<deltaAG<<G4endl;

    sumL10 += deltaL10;
    sumL96 += deltaL96;
    sumAG  += deltaAG;

    G4cout<<alpha1<<"\t"<<std::sqrt(alpha1)/degree<<"\t"
            <<sumL10<<"\t"<<sumL96<<"\t"<<sumAG<<G4endl;
      
    angleVector->PutValue( j-1 , alpha1, sumL10 ); // alpha2
  }
  fAngleTable->insertAt(i,angleVector);
  fAngleBank.push_back(fAngleTable);

  /*
  // Integral over all angle range - Bad accuracy !!!

  sumL10 = integral.Legendre10(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., alpha2);
  sumL96 = integral.Legendre96(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, 0., alpha2);
  sumAG  = integral.AdaptiveGauss(this, &G4NuclNuclDiffuseElastic::GetIntegrandFunction, 
                                       0., alpha2,epsilon);
  G4cout<<G4endl;
  G4cout<<alpha2<<"\t"<<std::sqrt(alpha2)/degree<<"\t"
            <<sumL10<<"\t"<<sumL96<<"\t"<<sumAG<<G4endl;
  */
  return;
}

/////////////////////////////////////////////////////////////////
//
//

 G4double G4NuclNuclDiffuseElastic::GetLegendrePol(G4int n, G4double theta)
{
  G4double legPol, epsilon = 1.e-6;
  G4double x = std::cos(theta);

  if     ( n  < 0 ) legPol = 0.;
  else if( n == 0 ) legPol = 1.;
  else if( n == 1 ) legPol = x;
  else if( n == 2 ) legPol = (3.*x*x-1.)/2.;
  else if( n == 3 ) legPol = (5.*x*x*x-3.*x)/2.;
  else if( n == 4 ) legPol = (35.*x*x*x*x-30.*x*x+3.)/8.;
  else if( n == 5 ) legPol = (63.*x*x*x*x*x-70.*x*x*x+15.*x)/8.;
  else if( n == 6 ) legPol = (231.*x*x*x*x*x*x-315.*x*x*x*x+105.*x*x-5.)/16.;
  else           
  {
    // legPol = ( (2*n-1)*x*GetLegendrePol(n-1,x) - (n-1)*GetLegendrePol(n-2,x) )/n;

    legPol = std::sqrt( 2./(n*CLHEP::pi*std::sin(theta+epsilon)) )*std::sin( (n+0.5)*theta+0.25*CLHEP::pi );
  }
  return legPol; 
}

/////////////////////////////////////////////////////////////////
//
//

G4complex G4NuclNuclDiffuseElastic::GetErfComp(G4complex z, G4int nMax)
{
  G4int n;
  G4double n2, cofn, shny, chny, fn, gn;

  G4double x = z.real();
  G4double y = z.imag();

  G4double outRe = 0., outIm = 0.;

  G4double twox  = 2.*x;
  G4double twoxy = twox*y;
  G4double twox2 = twox*twox;

  G4double cof1 = G4Exp(-x*x)/CLHEP::pi;

  G4double cos2xy = std::cos(twoxy);
  G4double sin2xy = std::sin(twoxy);

  G4double twoxcos2xy = twox*cos2xy;
  G4double twoxsin2xy = twox*sin2xy;

  for( n = 1; n <= nMax; n++)
  {
    n2   = n*n;

    cofn = G4Exp(-0.5*n2)/(n2+twox2);  // /(n2+0.5*twox2);

    chny = std::cosh(n*y);
    shny = std::sinh(n*y);

    fn  = twox - twoxcos2xy*chny + n*sin2xy*shny;
    gn  =        twoxsin2xy*chny + n*cos2xy*shny;

    fn *= cofn;
    gn *= cofn;

    outRe += fn;
    outIm += gn;
  }
  outRe *= 2*cof1;
  outIm *= 2*cof1;

  if(std::abs(x) < 0.0001)
  {
    outRe += GetErf(x);
    outIm += cof1*y;
  }
  else
  {
    outRe += GetErf(x) + cof1*(1-cos2xy)/twox;
    outIm += cof1*sin2xy/twox;
  }
  return G4complex(outRe, outIm);
}


/////////////////////////////////////////////////////////////////
//
//

G4complex G4NuclNuclDiffuseElastic::GetErfInt(G4complex z) // , G4int nMax)
{
  G4double outRe, outIm;

  G4double x = z.real();
  G4double y = z.imag();
  fReZ       = x;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  outRe = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetExpSin, 0., y );
  outIm = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetExpCos, 0., y );

  outRe *= 2./std::sqrt(CLHEP::pi);
  outIm *= 2./std::sqrt(CLHEP::pi);

  outRe += GetErf(x);

  return G4complex(outRe, outIm);
}

/////////////////////////////////////////////////////////////////
//
//


G4complex G4NuclNuclDiffuseElastic::GammaLess(G4double theta)
{
  G4double sinThetaR      = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double cosHalfThetaR2 = 1./(1. + fHalfRutThetaTg2);

  G4double u              = std::sqrt(0.5*fProfileLambda/sinThetaR);
  G4double kappa          = u/std::sqrt(CLHEP::pi);
  G4double dTheta         = theta - fRutherfordTheta;
  u                      *= dTheta;
  G4double u2             = u*u;
  G4double u2m2p3         = u2*2./3.;

  G4complex im            = G4complex(0.,1.);
  G4complex order         = G4complex(u,u);
  order                  /= std::sqrt(2.);

  G4complex gamma         = CLHEP::pi*kappa*GetErfcInt(-order)*std::exp(im*(u*u+0.25*CLHEP::pi));
  G4complex a0            = 0.5*(1. + 4.*(1.+im*u2)*cosHalfThetaR2/3.)/sinThetaR;
  G4complex a1            = 0.5*(1. + 2.*(1.+im*u2m2p3)*cosHalfThetaR2)/sinThetaR;
  G4complex out           = gamma*(1. - a1*dTheta) - a0;

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

G4complex G4NuclNuclDiffuseElastic::GammaMore(G4double theta)
{
  G4double sinThetaR      = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double cosHalfThetaR2 = 1./(1. + fHalfRutThetaTg2);

  G4double u              = std::sqrt(0.5*fProfileLambda/sinThetaR);
  G4double kappa          = u/std::sqrt(CLHEP::pi);
  G4double dTheta         = theta - fRutherfordTheta;
  u                      *= dTheta;
  G4double u2             = u*u;
  G4double u2m2p3         = u2*2./3.;

  G4complex im            = G4complex(0.,1.);
  G4complex order         = G4complex(u,u);
  order                  /= std::sqrt(2.);
  G4complex gamma         = CLHEP::pi*kappa*GetErfcInt(order)*std::exp(im*(u*u+0.25*CLHEP::pi));
  G4complex a0            = 0.5*(1. + 4.*(1.+im*u2)*cosHalfThetaR2/3.)/sinThetaR;
  G4complex a1            = 0.5*(1. + 2.*(1.+im*u2m2p3)*cosHalfThetaR2)/sinThetaR;
  G4complex out           = -gamma*(1. - a1*dTheta) - a0;

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

 G4complex G4NuclNuclDiffuseElastic::AmplitudeNear(G4double theta)
{
  G4double kappa = std::sqrt(0.5*fProfileLambda/std::sin(theta)/CLHEP::pi);
  G4complex out = G4complex(kappa/fWaveVector,0.);

  out *= PhaseNear(theta);

  if( theta <= fRutherfordTheta )
  {
    out *= GammaLess(theta) + ProfileNear(theta);
    // out *= GammaMore(theta) + ProfileNear(theta);
    out += CoulombAmplitude(theta);
  }
  else
  {
    out *= GammaMore(theta) + ProfileNear(theta);
    // out *= GammaLess(theta) + ProfileNear(theta);
  }
  return out;
}

/////////////////////////////////////////////////////////////////
//
//

G4complex G4NuclNuclDiffuseElastic::AmplitudeSim(G4double theta)
{
  G4double sinThetaR  = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double dTheta     = 0.5*(theta - fRutherfordTheta);
  G4double sindTheta  = std::sin(dTheta);
  G4double persqrt2   = std::sqrt(0.5);

  G4complex order     = G4complex(persqrt2,persqrt2);
  order              *= std::sqrt(0.5*fProfileLambda/sinThetaR)*2.*sindTheta;
  // order              *= std::sqrt(0.5*fProfileLambda/sinThetaR)*2.*dTheta;

  G4complex out;

  if ( theta <= fRutherfordTheta )
  {
    out = 1. - 0.5*GetErfcInt(-order)*ProfileNear(theta);
  }
  else
  {
    out = 0.5*GetErfcInt(order)*ProfileNear(theta);
  }

  out *= CoulombAmplitude(theta);

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

 G4complex G4NuclNuclDiffuseElastic::AmplitudeGla(G4double theta)
{
  G4int n;
  G4double T12b, b, b2; // cosTheta = std::cos(theta);
  G4complex out = G4complex(0.,0.), shiftC, shiftN; 
  G4complex im  = G4complex(0.,1.);

  for( n = 0; n < fMaxL; n++)
  {
    shiftC = std::exp( im*2.*CalculateCoulombPhase(n) );
    // b = ( fZommerfeld + std::sqrt( fZommerfeld*fZommerfeld + n*(n+1) ) )/fWaveVector;
    b = ( std::sqrt( G4double(n*(n+1)) ) )/fWaveVector;
    b2 = b*b;
    T12b = fSumSigma*G4Exp(-b2/fNuclearRadiusSquare)/CLHEP::pi/fNuclearRadiusSquare;         
    shiftN = std::exp( -0.5*(1.-im*fEtaRatio)*T12b ) - 1.;
    out +=  (2.*n+1.)*shiftC*shiftN*GetLegendrePol(n, theta);   
  }
  out /= 2.*im*fWaveVector;
  out += CoulombAmplitude(theta);
  return out;
}


/////////////////////////////////////////////////////////////////
//
//

 G4complex G4NuclNuclDiffuseElastic::AmplitudeGG(G4double theta)
{
  G4int n;
  G4double T12b, a, aTemp, b2, sinThetaH = std::sin(0.5*theta);
  G4double sinThetaH2 = sinThetaH*sinThetaH;
  G4complex out = G4complex(0.,0.); 
  G4complex im  = G4complex(0.,1.);

  a  = -fSumSigma/CLHEP::twopi/fNuclearRadiusSquare;
  b2 = fWaveVector*fWaveVector*fNuclearRadiusSquare*sinThetaH2;

  aTemp = a;

  for( n = 1; n < fMaxL; n++)
  {    
    T12b   = aTemp*G4Exp(-b2/n)/n;         
    aTemp *= a;  
    out   += T12b; 
    G4cout<<"out = "<<out<<G4endl;  
  }
  out *= -4.*im*fWaveVector/CLHEP::pi;
  out += CoulombAmplitude(theta);
  return out;
}


///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

void G4NuclNuclDiffuseElastic::InitParameters(const G4ParticleDefinition* theParticle,  
                                          G4double partMom, G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // atomic number
  fAtomicWeight  = A;     // number of nucleons

  fNuclearRadius2 = CalculateNuclearRad(fAtomicWeight);
  G4double A1     = G4double( theParticle->GetBaryonNumber() );   
  fNuclearRadius1 = CalculateNuclearRad(A1);
  // fNuclearRadius = std::sqrt(fNuclearRadius1*fNuclearRadius1+fNuclearRadius2*fNuclearRadius2);
  fNuclearRadius = fNuclearRadius1 + fNuclearRadius2;

  G4double a = 0.;
  G4double z = theParticle->GetPDGCharge();
  G4double m1 = theParticle->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double lambda = fCofLambda*fWaveVector*fNuclearRadius;
  G4cout<<"kR = "<<lambda<<G4endl;

  if( z )
  {
    a           = partMom/m1; // beta*gamma for m1
    fBeta       = a/std::sqrt(1+a*a);
    fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
    fRutherfordRatio = fZommerfeld/fWaveVector; 
    fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }
  G4cout<<"fZommerfeld = "<<fZommerfeld<<G4endl;
  fProfileLambda = lambda; // *std::sqrt(1.-2*fZommerfeld/lambda);
  G4cout<<"fProfileLambda = "<<fProfileLambda<<G4endl;
  fProfileDelta  = fCofDelta*fProfileLambda;
  fProfileAlpha  = fCofAlpha*fProfileLambda;

  CalculateCoulombPhaseZero();
  CalculateRutherfordAnglePar();

  return;
}
///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

void G4NuclNuclDiffuseElastic::InitDynParameters(const G4ParticleDefinition* theParticle,  
                                          G4double partMom) 
{
  G4double a = 0.;
  G4double z = theParticle->GetPDGCharge();
  G4double m1 = theParticle->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double lambda = fCofLambda*fWaveVector*fNuclearRadius;

  if( z )
  {
    a           = partMom/m1; // beta*gamma for m1
    fBeta       = a/std::sqrt(1+a*a);
    fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
    fRutherfordRatio = fZommerfeld/fWaveVector; 
    fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }
  fProfileLambda = lambda; // *std::sqrt(1.-2*fZommerfeld/lambda);
  fProfileDelta  = fCofDelta*fProfileLambda;
  fProfileAlpha  = fCofAlpha*fProfileLambda;

  CalculateCoulombPhaseZero();
  CalculateRutherfordAnglePar();

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

void G4NuclNuclDiffuseElastic::InitParametersGla(const G4DynamicParticle* aParticle,  
                                          G4double partMom, G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // target atomic number
  fAtomicWeight  = A;     // target number of nucleons

  fNuclearRadius2 = CalculateNuclearRad(fAtomicWeight); // target nucleus radius
  G4double A1     = G4double( aParticle->GetDefinition()->GetBaryonNumber() );   
  fNuclearRadius1 = CalculateNuclearRad(A1); // projectile nucleus radius
  fNuclearRadiusSquare = fNuclearRadius1*fNuclearRadius1+fNuclearRadius2*fNuclearRadius2;
 

  G4double a  = 0., kR12;
  G4double z  = aParticle->GetDefinition()->GetPDGCharge();
  G4double m1 = aParticle->GetDefinition()->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double pN = A1 - z;
  if( pN < 0. ) pN = 0.;

  G4double tN = A - Z;
  if( tN < 0. ) tN = 0.;

  G4double pTkin = aParticle->GetKineticEnergy();  
  pTkin /= A1;


  fSumSigma = (Z*z+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
              (z*tN+pN*Z)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  G4cout<<"fSumSigma = "<<fSumSigma/CLHEP::millibarn<<" mb"<<G4endl;
  G4cout<<"pi*R2 = "<<CLHEP::pi*fNuclearRadiusSquare/CLHEP::millibarn<<" mb"<<G4endl;
  kR12 = fWaveVector*std::sqrt(fNuclearRadiusSquare);
  G4cout<<"k*sqrt(R2) = "<<kR12<<" "<<G4endl;
  fMaxL = (G4int(kR12)+1)*4;
  G4cout<<"fMaxL = "<<fMaxL<<" "<<G4endl;

  if( z )
  {
      a           = partMom/m1; // beta*gamma for m1
      fBeta       = a/std::sqrt(1+a*a);
      fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
      fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }

  CalculateCoulombPhaseZero();
 

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns nucleon-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
// projectile nucleon is pParticle with pTkin shooting target nucleon tParticle

G4double 
G4NuclNuclDiffuseElastic::GetHadronNucleonXscNS( G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
                                                 G4ParticleDefinition* tParticle)
{
  G4double xsection(0), /*Delta,*/ A0, B0;
  G4double hpXsc(0);
  G4double hnXsc(0);


  G4double targ_mass     = tParticle->GetPDGMass();
  G4double proj_mass     = pParticle->GetPDGMass(); 

  G4double proj_energy   = proj_mass + pTkin; 
  G4double proj_momentum = std::sqrt(pTkin*(pTkin+2*proj_mass));

  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= CLHEP::GeV*CLHEP::GeV;  // in GeV for parametrisation
  proj_momentum /= CLHEP::GeV;
  proj_energy   /= CLHEP::GeV;
  proj_mass     /= CLHEP::GeV;
  G4double logS = G4Log(sMand);

  // General PDG fit constants


  // fEtaRatio=Re[f(0)]/Im[f(0)]

  if( proj_momentum >= 1.2 )
  {
    fEtaRatio  = 0.13*(logS - 5.8579332)*G4Pow::GetInstance()->powA(sMand,-0.18);
  }       
  else if( proj_momentum >= 0.6 )
  { 
    fEtaRatio = -75.5*(G4Pow::GetInstance()->powA(proj_momentum,0.25)-0.95)/
	  (G4Pow::GetInstance()->powA(3*proj_momentum,2.2)+1);     
  }
  else 
  {
    fEtaRatio = 15.5*proj_momentum/(27*proj_momentum*proj_momentum*proj_momentum+2);
  }
  G4cout<<"fEtaRatio = "<<fEtaRatio<<G4endl;

  // xsc
  
  if( proj_momentum >= 10. ) // high energy: pp = nn = np
    // if( proj_momentum >= 2.)
  {
    //Delta = 1.;

    //if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

    //AR-12Aug2016  if( proj_momentum >= 10.)
    {
        B0 = 7.5;
        A0 = 100. - B0*G4Log(3.0e7);

        xsection = A0 + B0*G4Log(proj_energy) - 11
                  + 103*G4Pow::GetInstance()->powA(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
    }
  }
  else // low energy pp = nn != np
  {
      if(pParticle == tParticle) // pp or nn      // nn to be pp
      {
        if( proj_momentum < 0.73 )
        {
          hnXsc = 23 + 50*( G4Pow::GetInstance()->powA( G4Log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hnXsc = 23 + 40*(G4Log(proj_momentum/0.73))*
                         (G4Log(proj_momentum/0.73));
        }
        else  // if( proj_momentum < 10.  )
        {
          hnXsc = 39.0 + 
              75*(proj_momentum - 1.2)/(G4Pow::GetInstance()->powA(proj_momentum,3.0) + 0.15);
        }
        xsection = hnXsc;
      }
      else  // pn to be np
      {
        if( proj_momentum < 0.8 )
        {
          hpXsc = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hpXsc = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/0.95),2.0);
        }
        else    // if( proj_momentum < 10.  )
        {
          hpXsc = 33.3+
              20.8*(G4Pow::GetInstance()->powA(proj_momentum,2.0)-1.35)/
                 (G4Pow::GetInstance()->powA(proj_momentum,2.50)+0.95);
        }
        xsection = hpXsc;
      }
  }
  xsection *= CLHEP::millibarn; // parametrised in mb
  G4cout<<"xsection = "<<xsection/CLHEP::millibarn<<" mb"<<G4endl;
  return xsection;
}

/////////////////////////////////////////////////////////////////
//
// The ratio el/ruth for Fresnel smooth nucleus profile

G4double G4NuclNuclDiffuseElastic::GetRatioGen(G4double theta)
{
  G4double sinThetaR  = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double dTheta     = 0.5*(theta - fRutherfordTheta);
  G4double sindTheta  = std::sin(dTheta);

  G4double prof       = Profile(theta);
  G4double prof2      = prof*prof;
  // G4double profmod    = std::abs(prof);
  G4double order      = std::sqrt(fProfileLambda/sinThetaR/CLHEP::pi)*2.*sindTheta;

  order = std::abs(order); // since sin changes sign!
  // G4cout<<"order = "<<order<<G4endl; 
 
  G4double cosFresnel = GetCint(order);  
  G4double sinFresnel = GetSint(order);  

  G4double out;

  if ( theta <= fRutherfordTheta )
  {
    out  = 1. + 0.5*( (0.5-cosFresnel)*(0.5-cosFresnel)+(0.5-sinFresnel)*(0.5-sinFresnel) )*prof2; 
    out += ( cosFresnel + sinFresnel - 1. )*prof;
  }
  else
  {
    out = 0.5*( (0.5-cosFresnel)*(0.5-cosFresnel)+(0.5-sinFresnel)*(0.5-sinFresnel) )*prof2;
  }

  return out;
}

///////////////////////////////////////////////////////////////////
//
// For the calculation of arg Gamma(z) one needs complex extension 
// of ln(Gamma(z))

G4complex G4NuclNuclDiffuseElastic::GammaLogarithm(G4complex zz)
{
  const G4double cof[6] = { 76.18009172947146,     -86.50532032941677,
                             24.01409824083091,      -1.231739572450155,
                              0.1208650973866179e-2, -0.5395239384953e-5  } ;
  G4int j;
  G4complex z = zz - 1.0;
  G4complex tmp = z + 5.5;
  tmp -= (z + 0.5) * std::log(tmp);
  G4complex ser = G4complex(1.000000000190015,0.);

  for ( j = 0; j <= 5; j++ )
  {
    z += 1.0;
    ser += cof[j]/z;
  }
  return -tmp + std::log(2.5066282746310005*ser);
}

/////////////////////////////////////////////////////////////
//
// Bessel J0 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

G4double G4NuclNuclDiffuseElastic::BesselJzero(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = std::fabs(value);

  if ( value < 8.0 && value > -8.0 )
  {
    value2 = value*value;

    fact1  = 57568490574.0 + value2*(-13362590354.0 
                           + value2*( 651619640.7 
                           + value2*(-11214424.18 
                           + value2*( 77392.33017 
                           + value2*(-184.9052456   ) ) ) ) );

    fact2  = 57568490411.0 + value2*( 1029532985.0 
                           + value2*( 9494680.718
                           + value2*(59272.64853
                           + value2*(267.8532712 
                           + value2*1.0               ) ) ) );

    bessel = fact1/fact2;
  } 
  else 
  {
    arg    = 8.0/modvalue;

    value2 = arg*arg;

    shift  = modvalue-0.785398164;

    fact1  = 1.0 + value2*(-0.1098628627e-2 
                 + value2*(0.2734510407e-4
                 + value2*(-0.2073370639e-5 
                 + value2*0.2093887211e-6    ) ) );

    fact2  = -0.1562499995e-1 + value2*(0.1430488765e-3
                              + value2*(-0.6911147651e-5 
                              + value2*(0.7621095161e-6
                              - value2*0.934945152e-7    ) ) );

    bessel = std::sqrt(0.636619772/modvalue)*(std::cos(shift)*fact1 - arg*std::sin(shift)*fact2 );
  }
  return bessel;
}

/////////////////////////////////////////////////////////////
//
// Bessel J1 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

G4double G4NuclNuclDiffuseElastic::BesselJone(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = std::fabs(value);

  if ( modvalue < 8.0 ) 
  {
    value2 = value*value;

    fact1  = value*(72362614232.0 + value2*(-7895059235.0 
                                  + value2*( 242396853.1
                                  + value2*(-2972611.439 
                                  + value2*( 15704.48260 
                                  + value2*(-30.16036606  ) ) ) ) ) );

    fact2  = 144725228442.0 + value2*(2300535178.0 
                            + value2*(18583304.74
                            + value2*(99447.43394 
                            + value2*(376.9991397 
                            + value2*1.0             ) ) ) );
    bessel = fact1/fact2;
  } 
  else 
  {
    arg    = 8.0/modvalue;

    value2 = arg*arg;

    shift  = modvalue - 2.356194491;

    fact1  = 1.0 + value2*( 0.183105e-2 
                 + value2*(-0.3516396496e-4
                 + value2*(0.2457520174e-5 
                 + value2*(-0.240337019e-6          ) ) ) );

    fact2 = 0.04687499995 + value2*(-0.2002690873e-3
                          + value2*( 0.8449199096e-5
                          + value2*(-0.88228987e-6
                          + value2*0.105787412e-6       ) ) );

    bessel = std::sqrt( 0.636619772/modvalue)*(std::cos(shift)*fact1 - arg*std::sin(shift)*fact2);

    if (value < 0.0) bessel = -bessel;
  }
  return bessel;
}

//
//
/////////////////////////////////////////////////////////////////////////////////
