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
//
// Physics model class G4DiffuseElasticV2 
//
//
// G4 Model: optical diffuse elastic scattering with 4-momentum balance
//                         
// 24-May-07 V. Grichine
//
// 21.10.15 V. Grichine 
//             Bug fixed in BuildAngleTable, improving accuracy for 
//             angle bins at high energies > 50 GeV for pions.
//
// 24.11.17 W. Pokorski, code cleanup and performance improvements
//

#include "G4DiffuseElasticV2.hh"
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

#include "G4Exp.hh"

#include "G4HadronicParameters.hh"

/////////////////////////////////////////////////////////////////////////
//


G4DiffuseElasticV2::G4DiffuseElasticV2() 
  : G4HadronElastic("DiffuseElasticV2"), fParticle(0)
{
  SetMinEnergy( 0.01*MeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  verboseLevel         = 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ      = 0.0*GeV;  
  lowEnergyLimitHE     = 0.0*GeV;  
  lowestEnergyLimit    = 0.0*keV;  
  plabLowLimit         = 20.0*MeV;

  theProton    = G4Proton::Proton();
  theNeutron   = G4Neutron::Neutron();

  fEnergyBin = 300;  // Increased from the original 200 to have no wider log-energy-bins up to 10 PeV  
  fAngleBin  = 200;

  fEnergyVector =  new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );

  fEnergyAngleVector = 0;
  fEnergySumVector = 0;
  
  fParticle      = 0;
  fWaveVector    = 0.;
  fAtomicWeight  = 0.;
  fAtomicNumber  = 0.;
  fNuclearRadius = 0.;
  fBeta          = 0.;
  fZommerfeld    = 0.;
  fAm = 0.;
  fAddCoulomb = false;
}

//////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4DiffuseElasticV2::~G4DiffuseElasticV2()
{
  if ( fEnergyVector ) 
  {
    delete fEnergyVector;
    fEnergyVector = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// Initialisation for given particle using element table of application

void G4DiffuseElasticV2::Initialise() 
{

  const G4ElementTable* theElementTable = G4Element::GetElementTable();

  std::size_t jEl, numOfEl = G4Element::GetNumberOfElements();

  for( jEl = 0; jEl < numOfEl; ++jEl) // application element loop
  {
    fAtomicNumber = (*theElementTable)[jEl]->GetZ();     // atomic number
    fAtomicWeight = G4NistManager::Instance()->GetAtomicMassAmu( static_cast< G4int >( fAtomicNumber ) );
    fNuclearRadius = CalculateNuclearRad(fAtomicWeight);

    if( verboseLevel > 0 ) 
    {   
      G4cout<<"G4DiffuseElasticV2::Initialise() the element: "
	    <<(*theElementTable)[jEl]->GetName()<<G4endl;
    }
    fElementNumberVector.push_back(fAtomicNumber);
    fElementNameVector.push_back((*theElementTable)[jEl]->GetName());

    BuildAngleTable();

    fEnergyAngleVectorBank.push_back(fEnergyAngleVector);
    fEnergySumVectorBank.push_back(fEnergySumVector);

  }  
  return;
}

////////////////////////////////////////////////////////////////////////////
//
// return differential elastic probability d(probability)/d(t) with 
// Coulomb correction. It is called from BuildAngleTable()

G4double 
G4DiffuseElasticV2::GetDiffElasticSumProbA( G4double theta )
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

  if ( fParticle == theProton )
  {
    diffuse = 0.63*fermi;
    gamma   = 0.3*fermi;
    delta   = 0.1*fermi*fermi;
    e1      = 0.3*fermi;
    e2      = 0.35*fermi;
  }
  else if ( fParticle == theNeutron )
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
  
  G4double lambda = 15; // 15 ok
  // G4double kgamma    = fWaveVector*gamma;   // wavek*delta;
  G4double kgamma    = lambda*(1.-G4Exp(-fWaveVector*gamma/lambda));   // wavek*delta;

  if( fAddCoulomb )  // add Coulomb correction
  {
    G4double sinHalfTheta  = std::sin(0.5*theta); 
    G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta;

    kgamma += 0.5*fZommerfeld/kr/(sinHalfTheta2+fAm); // correction at J0()
  }
  
  G4double kgamma2   = kgamma*kgamma;

  // G4double dk2t  = delta*fWaveVector*fWaveVector*theta; // delta*wavek*wavek*theta;
  // G4double dk2t2 = dk2t*dk2t;

  // G4double pikdt = pi*fWaveVector*diffuse*theta;// pi*wavek*diffuse*theta;
  G4double pikdt    = lambda*(1. - G4Exp( -pi*fWaveVector*diffuse*theta/lambda ) );   // wavek*delta;

  damp           = DampFactor( pikdt );
  damp2          = damp*damp;

  G4double mode2k2 = ( e1*e1 + e2*e2 )*fWaveVector*fWaveVector;  
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
G4DiffuseElasticV2::GetIntegrandFunction( G4double alpha )
{
  G4double result;

  result  = GetDiffElasticSumProbA(alpha) * 2 * CLHEP::pi * std::sin(alpha);
  
  return result;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////  Table preparation and reading ////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4DiffuseElasticV2::SampleInvariantT( const G4ParticleDefinition* aParticle, G4double p, 
                                               G4int Z, G4int A)
{
  fParticle = aParticle;
  G4double m1 = fParticle->GetPDGMass(), t;
  G4double totElab = std::sqrt(m1*m1+p*p);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv1(p,0.0,0.0,totElab);
  G4LorentzVector  lv(0.0,0.0,0.0,mass2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double momentumCMS = p1.mag();
  
  if( aParticle == theNeutron)
  {
    G4double Tmax = NeutronTuniform( Z );
    G4double pCMS2 = momentumCMS*momentumCMS;
    G4double Tkin = std::sqrt(pCMS2+m1*m1)-m1;

    if( Tkin <= Tmax )
    {
      t = 4.*pCMS2*G4UniformRand();
      return t;
    }
  }
  
  t = SampleTableT( aParticle,  momentumCMS, G4double(Z), G4double(A) ); // sample theta in cms

  return t;
}

///////////////////////////////////////////////////////

G4double G4DiffuseElasticV2::NeutronTuniform(G4int Z)
{
  G4double elZ  = G4double(Z);
  elZ -= 1.;
  G4double Tkin = 12.*G4Exp(-elZ/10.) + 1.;

  return Tkin;
}


////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4DiffuseElasticV2::SampleTableT( const G4ParticleDefinition* aParticle, G4double p, 
                                               G4double Z, G4double A)
{
  G4double alpha = SampleTableThetaCMS( aParticle,  p, Z, A); // sample theta in cms
  G4double t     = 2*p*p*( 1 - std::cos(alpha) );             // -t !!!
  
  return t;
}

////////////////////////////////////////////////////////////////////////////
//
// Return scattering angle2 sampled in cms according to precalculated table.


G4double 
G4DiffuseElasticV2::SampleTableThetaCMS(const G4ParticleDefinition* particle, 
                                       G4double momentum, G4double Z, G4double A)
{
  std::size_t iElement;
  G4int iMomentum;
  unsigned long iAngle = 0;  
  G4double randAngle, position, theta1, theta2, E1, E2, W1, W2, W;  
  G4double m1 = particle->GetPDGMass();

  for(iElement = 0; iElement < fElementNumberVector.size(); iElement++)
  {
    if( std::fabs(Z - fElementNumberVector[iElement]) < 0.5) break;
  }

  if ( iElement == fElementNumberVector.size() ) 
  {
    InitialiseOnFly(Z,A); // table preparation, if needed
  }
  
  fEnergyAngleVector = fEnergyAngleVectorBank[iElement];
  fEnergySumVector = fEnergySumVectorBank[iElement];

  
  G4double kinE = std::sqrt(momentum*momentum + m1*m1) - m1;
    
  iMomentum = G4int(fEnergyVector->FindBin(kinE,1000) + 1);
  
  position = (*(*fEnergySumVector)[iMomentum])[0]*G4UniformRand();

  for(iAngle = 0; iAngle < fAngleBin; ++iAngle)
    {
      if (position > (*(*fEnergySumVector)[iMomentum])[iAngle]) break;
    }

  
  if (iMomentum == fEnergyBin -1 || iMomentum == 0 )   // the table edges
    {     
      randAngle = GetScatteringAngle(iMomentum, iAngle, position);
    }
  else  // kinE inside between energy table edges
    {                  
      theta2  = GetScatteringAngle(iMomentum, iAngle, position);
      
      E2 = fEnergyVector->Energy(iMomentum);
      
      iMomentum--;
      
      theta1  = GetScatteringAngle(iMomentum, iAngle, position);
    
      E1 = fEnergyVector->Energy(iMomentum);
    
      W  = 1.0/(E2 - E1);
      W1 = (E2 - kinE)*W;
      W2 = (kinE - E1)*W;

      randAngle = W1*theta1 + W2*theta2;
    }



  if(randAngle < 0.) randAngle = 0.;

  return randAngle;
}

//////////////////////////////////////////////////////////////////////////////
//
// Initialisation for given particle on fly using new element number

void G4DiffuseElasticV2::InitialiseOnFly(G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // atomic number
  fAtomicWeight  = G4NistManager::Instance()->GetAtomicMassAmu( static_cast< G4int >( Z ) );

  fNuclearRadius = CalculateNuclearRad(fAtomicWeight);

  if( verboseLevel > 0 )    
  {
    G4cout<<"G4DiffuseElasticV2::InitialiseOnFly() the element with Z = "
	  <<Z<<"; and A = "<<A<<G4endl;
  }
  fElementNumberVector.push_back(fAtomicNumber);

  BuildAngleTable();

  fEnergyAngleVectorBank.push_back(fEnergyAngleVector);
  fEnergySumVectorBank.push_back(fEnergySumVector);
  
  return;
}

///////////////////////////////////////////////////////////////////////////////
//
// Build for given particle and element table of momentum, angle probability.
// For the moment in lab system. 

void G4DiffuseElasticV2::BuildAngleTable() 
{
  G4double partMom, kinE, a = 0., z = fParticle->GetPDGCharge(), m1 = fParticle->GetPDGMass();
  G4double alpha1, alpha2, alphaMax, alphaCoulomb, delta = 0., sum = 0.;

  G4Integrator<G4DiffuseElasticV2,G4double(G4DiffuseElasticV2::*)(G4double)> integral;
  
  fEnergyAngleVector = new std::vector<std::vector<G4double>*>;
  fEnergySumVector = new std::vector<std::vector<G4double>*>;
  
  for( G4int i = 0; i < fEnergyBin; ++i)
  {
    kinE        = fEnergyVector->Energy(i);
    partMom     = std::sqrt( kinE*(kinE + 2*m1) );

    fWaveVector = partMom/hbarc;

    G4double kR     = fWaveVector*fNuclearRadius;
    G4double kRmax  = 18.6; // 10.6; 10.6, 18, 10.174; ~ 3 maxima of J1 or 15., 25.
    G4double kRcoul = 1.9; // 1.2; 1.4, 2.5; // on the first slope of J1

    alphaMax = kRmax/kR;

    if ( alphaMax >= CLHEP::pi ) alphaMax = CLHEP::pi;   // vmg21.10.15 

    alphaCoulomb = kRcoul/kR;

    if( z )
    {
      a           = partMom/m1;         // beta*gamma for m1
      fBeta       = a/std::sqrt(1+a*a);
      fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
      fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
      fAddCoulomb = true;
    }

    std::vector<G4double>* angleVector = new std::vector<G4double>(fAngleBin);
    std::vector<G4double>* sumVector = new std::vector<G4double>(fAngleBin);

    
    G4double delth = alphaMax/fAngleBin;
        
    sum = 0.;
 
    for(G4int j = (G4int)fAngleBin-1; j >= 0; --j)
    {
      alpha1 = delth*j;
      alpha2 = alpha1 + delth;
      
      if( fAddCoulomb && ( alpha2 < alphaCoulomb)) fAddCoulomb = false;

      delta = integral.Legendre10(this, &G4DiffuseElasticV2::GetIntegrandFunction, alpha1, alpha2);

      sum += delta;
      
      (*angleVector)[j] = alpha1;
      (*sumVector)[j] = sum; 
      
    }
    fEnergyAngleVector->push_back(angleVector);
    fEnergySumVector->push_back(sumVector);
      
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4DiffuseElasticV2::GetScatteringAngle( G4int iMomentum, unsigned long iAngle, G4double position )
{
  G4double x1, x2, y1, y2, randAngle = 0;
  
  if( iAngle == 0 )
  {
    randAngle = (*(*fEnergyAngleVector)[iMomentum])[iAngle];
  }
  else
  {
    if ( iAngle >= (*fEnergyAngleVector)[iMomentum]->size() )
    {
      iAngle = (*fEnergyAngleVector)[iMomentum]->size() - 1;
    }
    
    y1 = (*(*fEnergySumVector)[iMomentum])[iAngle-1];
    y2 = (*(*fEnergySumVector)[iMomentum])[iAngle];

    x1 = (*(*fEnergyAngleVector)[iMomentum])[iAngle-1];
    x2 = (*(*fEnergyAngleVector)[iMomentum])[iAngle];

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
// Return scattering angle in lab system (target at rest) knowing theta in CMS



G4double 
G4DiffuseElasticV2::ThetaCMStoThetaLab( const G4DynamicParticle* aParticle, 
                                        G4double tmass, G4double thetaCMS)
{
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
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
G4DiffuseElasticV2::ThetaLabToThetaCMS( const G4DynamicParticle* aParticle, 
                                        G4double tmass, G4double thetaLab)
{
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = aParticle->GetTotalMomentum();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,tmass);   

  lv += lv1;

  G4ThreeVector bst = lv.boostVector();

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

