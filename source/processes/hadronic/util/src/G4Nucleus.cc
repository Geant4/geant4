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
// original by H.P. Wellisch
// modified by J.L. Chuma, TRIUMF, 19-Nov-1996
// last modified: 27-Mar-1997
// J.P.Wellisch: 23-Apr-97: minor simplifications
// modified by J.L.Chuma 24-Jul-97  to set the total momentum in Cinema and
//                                  EvaporationEffects
// modified by J.L.Chuma 21-Oct-97  put std::abs() around the totalE^2-mass^2
//                                  in calculation of total momentum in
//                                  Cinema and EvaporationEffects
// Chr. Volcker, 10-Nov-1997: new methods and class variables.
// HPW added utilities for low energy neutron transport. (12.04.1998)
// M.G. Pia, 2 Oct 1998: modified GetFermiMomentum to avoid memory leaks
// G.Folger, spring 2010:  add integer A/Z interface
// A. Ribon, summer 2015:  migrated to G4Exp and G4Log
// A. Ribon, autumn 2021:  extended to hypernuclei
 
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4HyperNucleiProperties.hh"
#include "G4HadronicParameters.hh"


G4Nucleus::G4Nucleus()
  : theA(0), theZ(0), theL(0), aEff(0.0), zEff(0)
{
  pnBlackTrackEnergy = 0.0;
  dtaBlackTrackEnergy = 0.0;
  pnBlackTrackEnergyfromAnnihilation = 0.0;
  dtaBlackTrackEnergyfromAnnihilation = 0.0;
  excitationEnergy = 0.0;
  momentum = G4ThreeVector(0.,0.,0.);
  fermiMomentum = 1.52*hbarc/fermi;
  theTemp = 293.16*kelvin;
  fIsotope = 0;
}

G4Nucleus::G4Nucleus( const G4double A, const G4double Z, const G4int numberOfLambdas )
{
  SetParameters( A, Z, std::max(numberOfLambdas, 0) );
  pnBlackTrackEnergy = 0.0;
  dtaBlackTrackEnergy = 0.0;
  pnBlackTrackEnergyfromAnnihilation = 0.0;
  dtaBlackTrackEnergyfromAnnihilation = 0.0;
  excitationEnergy = 0.0;
  momentum = G4ThreeVector(0.,0.,0.);
  fermiMomentum = 1.52*hbarc/fermi;
  theTemp = 293.16*kelvin;
  fIsotope = 0;
}

G4Nucleus::G4Nucleus( const G4int A, const G4int Z, const G4int numberOfLambdas )
{
  SetParameters( A, Z, std::max(numberOfLambdas, 0) );
  pnBlackTrackEnergy = 0.0;
  dtaBlackTrackEnergy = 0.0;
  pnBlackTrackEnergyfromAnnihilation = 0.0;
  dtaBlackTrackEnergyfromAnnihilation = 0.0;
  excitationEnergy = 0.0;
  momentum = G4ThreeVector(0.,0.,0.);
  fermiMomentum = 1.52*hbarc/fermi;
  theTemp = 293.16*kelvin;
  fIsotope = 0;
}

G4Nucleus::G4Nucleus( const G4Material *aMaterial )
{
  ChooseParameters( aMaterial );
  pnBlackTrackEnergy = 0.0;
  dtaBlackTrackEnergy = 0.0;
  pnBlackTrackEnergyfromAnnihilation = 0.0;
  dtaBlackTrackEnergyfromAnnihilation = 0.0;
  excitationEnergy = 0.0;
  momentum = G4ThreeVector(0.,0.,0.);
  fermiMomentum = 1.52*hbarc/fermi;
  theTemp = aMaterial->GetTemperature();
  fIsotope = 0;
}

G4Nucleus::~G4Nucleus() {}


//-------------------------------------------------------------------------------------------------
// SVT (Sampling of the Velocity of the Target nucleus) method, L. Thulliez (CEA-Saclay) 2021/05/04
//-------------------------------------------------------------------------------------------------
G4ReactionProduct 
G4Nucleus::GetBiasedThermalNucleus(G4double aMass, G4ThreeVector aVelocity, G4double temp) const
{
  // If E_neutron <= E_threshold, Then apply the Sampling ot the Velocity of the Target (SVT) method;
  // Else consider the target nucleus being without motion.
  G4double E_threshold = G4HadronicParameters::Instance()->GetNeutronKineticEnergyThresholdForSVT();
  if ( E_threshold == -1. ) {
    E_threshold = 400.0*8.617333262E-11*temp;
  }
  G4double E_neutron = 0.5*aVelocity.mag2()*G4Neutron::Neutron()->GetPDGMass(); // E=0.5*m*v2

  G4ReactionProduct result;
  result.SetMass(aMass*G4Neutron::Neutron()->GetPDGMass());

  if ( E_neutron <= E_threshold ) { 
    
    // Beta = sqrt(m/2kT)
    G4double beta = std::sqrt(result.GetMass()/(2.*8.617333262E-11*temp)); // kT E-5[eV] mass E-11[MeV] => beta in [m/s]-1
	
    // Neutron speed vn
    G4double vN_norm = aVelocity.mag();
    G4double vN_norm2 = vN_norm*vN_norm;
    G4double y = beta*vN_norm;

    // Normalize neutron velocity
    aVelocity = (1./vN_norm)*aVelocity;

    // Sample target speed
    G4double x2;
    G4double randThreshold;
    G4double vT_norm, vT_norm2, mu; //theta, val1, val2,
    G4double acceptThreshold;
    G4double vRelativeSpeed;
    G4double cdf0 = 2./(2.+std::sqrt(CLHEP::pi)*y);

    do {
      // Sample the target velocity vT in the laboratory frame
      if ( G4UniformRand() < cdf0 ) {
        // Sample in C45 from https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf
        x2 = -std::log(G4UniformRand()*G4UniformRand());
      } else {
        // Sample in C61 from https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf
        G4double ampl = std::cos(CLHEP::pi/2.0 * G4UniformRand());
        x2 = -std::log(G4UniformRand()) - std::log(G4UniformRand())*ampl*ampl;
      }

      vT_norm = std::sqrt(x2)/beta;
      vT_norm2 = vT_norm*vT_norm;
		
      // Sample cosine between the incident neutron and the target in the laboratory frame
      mu = 2*G4UniformRand() - 1;

      // Define acceptance threshold
      vRelativeSpeed = std::sqrt(vN_norm2 + vT_norm2 - 2*vN_norm*vT_norm*mu);
      acceptThreshold = vRelativeSpeed/(vN_norm + vT_norm);
      randThreshold = G4UniformRand();
    } while ( randThreshold >= acceptThreshold );

    DoKinematicsOfThermalNucleus(mu, vT_norm, aVelocity, result);
    
  } else { // target nucleus considered as being without motion
    
    result.SetMomentum(0., 0., 0.);
    result.SetKineticEnergy(0.);

  }

  return result;
}


void
G4Nucleus::DoKinematicsOfThermalNucleus(const G4double mu, const G4double vT_norm, const G4ThreeVector& aVelocity,
                                        G4ReactionProduct& result) const {

  // Get target nucleus direction from the neutron direction and the relative angle between target nucleus and neutron (mu)
  G4double cosTh = mu;
  G4ThreeVector uNorm = aVelocity;
  	
  G4double sinTh = std::sqrt(1. - cosTh*cosTh);
  	
  // Sample randomly the phi angle between the neutron veloicty and the target velocity
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double sinPhi = std::sin(phi);
  G4double cosPhi = std::cos(phi);

  // Find orthogonal vector to aVelocity - solve equation xx' + yy' + zz' = 0
  G4ThreeVector ortho(1., 1., 1.);
  if      ( uNorm[0] )  ortho[0] = -(uNorm[1]+uNorm[2])/uNorm[0];
  else if ( uNorm[1] )  ortho[1] = -(uNorm[0]+uNorm[2])/uNorm[1];
  else if ( uNorm[2] )  ortho[2] = -(uNorm[0]+uNorm[1])/uNorm[2];

  // Normalize the vector
  ortho = (1/ortho.mag())*ortho;
  	
  // Find vector to draw a plan perpendicular to uNorm (i.e neutron velocity) with vectors ortho & orthoComp
  G4ThreeVector orthoComp( uNorm[1]*ortho[2] - ortho[1]*uNorm[2],
                           uNorm[2]*ortho[0] - ortho[2]*uNorm[0],
                           uNorm[0]*ortho[1] - ortho[0]*uNorm[1] );

  // Find the direction of the target velocity in the laboratory frame
  G4ThreeVector directionTarget( cosTh*uNorm[0] + sinTh*(cosPhi*orthoComp[0] + sinPhi*ortho[0]),
                                 cosTh*uNorm[1] + sinTh*(cosPhi*orthoComp[1] + sinPhi*ortho[1]),
                                 cosTh*uNorm[2] + sinTh*(cosPhi*orthoComp[2] + sinPhi*ortho[2]) );
  	
  // Normalize directionTarget
  directionTarget = ( 1./directionTarget.mag() )*directionTarget;

  // Set momentum
  G4double px = result.GetMass()*vT_norm*directionTarget[0];
  G4double py = result.GetMass()*vT_norm*directionTarget[1];
  G4double pz = result.GetMass()*vT_norm*directionTarget[2];
  result.SetMomentum(px, py, pz);

  G4double tMom = std::sqrt(px*px+py*py+pz*pz);
  G4double tEtot = std::sqrt( (tMom+result.GetMass())*(tMom+result.GetMass())
  		              - 2.*tMom*result.GetMass() );

  if ( tEtot/result.GetMass() - 1. > 0.001 ) {
    // use relativistic energy for higher energies
    result.SetTotalEnergy(tEtot);
  } else {
    // use p**2/2M for lower energies (to preserve precision?)
    result.SetKineticEnergy(tMom*tMom/(2.*result.GetMass()));
  }

}


G4ReactionProduct
G4Nucleus::GetThermalNucleus(G4double targetMass, G4double temp) const
{
  G4double currentTemp = temp;
  if (currentTemp < 0) currentTemp = theTemp;
  G4ReactionProduct theTarget;    
  theTarget.SetMass(targetMass*G4Neutron::Neutron()->GetPDGMass());
  G4double px, py, pz;
  px = GetThermalPz(theTarget.GetMass(), currentTemp);
  py = GetThermalPz(theTarget.GetMass(), currentTemp);
  pz = GetThermalPz(theTarget.GetMass(), currentTemp);
  theTarget.SetMomentum(px, py, pz);
  G4double tMom = std::sqrt(px*px+py*py+pz*pz);
  G4double tEtot = std::sqrt((tMom+theTarget.GetMass())*
                             (tMom+theTarget.GetMass())-
                              2.*tMom*theTarget.GetMass());
  //  if(1-tEtot/theTarget.GetMass()>0.001)  this line incorrect (Bug report 1911) 
  if (tEtot/theTarget.GetMass() - 1. > 0.001) {
    // use relativistic energy for higher energies
    theTarget.SetTotalEnergy(tEtot);

  } else {
    // use p**2/2M for lower energies (to preserve precision?) 
    theTarget.SetKineticEnergy(tMom*tMom/(2.*theTarget.GetMass()));
  }    
  return theTarget;
}

 
void
G4Nucleus::ChooseParameters(const G4Material* aMaterial)
{
  G4double random = G4UniformRand();
  G4double sum = aMaterial->GetTotNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  G4double running(0);
  //  G4Element* element(0);
  const G4Element* element = (*theElementVector)[aMaterial->GetNumberOfElements()-1];

  for (unsigned int i = 0; i < aMaterial->GetNumberOfElements(); ++i) {
    running += aMaterial->GetVecNbOfAtomsPerVolume()[i];
    if (running > random*sum) {
      element = (*theElementVector)[i];
      break;
    }
  }

  if (element->GetNumberOfIsotopes() > 0) {
    G4double randomAbundance = G4UniformRand();
    G4double sumAbundance = element->GetRelativeAbundanceVector()[0];
    unsigned int iso=0;
    while (iso < element->GetNumberOfIsotopes() &&  /* Loop checking, 02.11.2015, A.Ribon */
           sumAbundance < randomAbundance) {
      ++iso;
      sumAbundance += element->GetRelativeAbundanceVector()[iso];
    }
    theA=element->GetIsotope(iso)->GetN();
    theZ=element->GetIsotope(iso)->GetZ();
    theL=0;
    aEff=theA;
    zEff=theZ;
  } else {   
    aEff = element->GetN();
    zEff = element->GetZ();
    theZ = G4int(zEff + 0.5);
    theA = G4int(aEff + 0.5);
    theL=0;
  }
}


void
G4Nucleus::SetParameters( const G4double A, const G4double Z, const G4int numberOfLambdas )
{
  theZ = G4lrint(Z);
  theA = G4lrint(A);
  theL = std::max(numberOfLambdas, 0);
  if (theA<1 || theZ<0 || theZ>theA) {
    throw G4HadronicException(__FILE__, __LINE__,
            "G4Nucleus::SetParameters called with non-physical parameters");
  }
  aEff = A;  // atomic weight
  zEff = Z;  // atomic number
  fIsotope = 0;
}


void
G4Nucleus::SetParameters( const G4int A, const G4int Z, const G4int numberOfLambdas )
{
  theZ = Z;
  theA = A;
  theL = std::max(numberOfLambdas, 0);
  if( theA<1 || theZ<0 || theZ>theA )
    {
      throw G4HadronicException(__FILE__, __LINE__,
				"G4Nucleus::SetParameters called with non-physical parameters");
    }
  aEff = A;  // atomic weight
  zEff = Z;  // atomic number
  fIsotope = 0;
}


G4DynamicParticle *
G4Nucleus::ReturnTargetParticle() const
{
  // choose a proton or a neutron (or a lamba if a hypernucleus) as the target particle
  G4DynamicParticle *targetParticle = new G4DynamicParticle;
  const G4double rnd = G4UniformRand();
  if ( rnd < zEff/aEff ) {
    targetParticle->SetDefinition( G4Proton::Proton() );
  } else if ( rnd < (zEff + theL*1.0)/aEff ) {
    targetParticle->SetDefinition( G4Lambda::Lambda() );
  } else {
    targetParticle->SetDefinition( G4Neutron::Neutron() );
  }
  return targetParticle;
}

 
G4double
G4Nucleus::AtomicMass( const G4double A, const G4double Z, const G4int numberOfLambdas ) const
{
  // Now returns (atomic mass - electron masses)
  if ( numberOfLambdas > 0 ) {
    return G4HyperNucleiProperties::GetNuclearMass(G4int(A), G4int(Z), numberOfLambdas);
  } else {
    return G4NucleiProperties::GetNuclearMass(A, Z);
  }
}

 
G4double
G4Nucleus::AtomicMass( const G4int A, const G4int Z, const G4int numberOfLambdas ) const
{
  // Now returns (atomic mass - electron masses)
  if ( numberOfLambdas > 0 ) {
    return G4HyperNucleiProperties::GetNuclearMass(A, Z, numberOfLambdas);
  } else {
    return G4NucleiProperties::GetNuclearMass(A, Z);
  }
}
 
 
G4double
G4Nucleus::GetThermalPz( const G4double mass, const G4double temp ) const
{
  G4double result = G4RandGauss::shoot();
  result *= std::sqrt(k_Boltzmann*temp*mass); // Das ist impuls (Pz),
                                         // nichtrelativistische rechnung
                                         // Maxwell verteilung angenommen
  return result;
}
 

G4double 
G4Nucleus::EvaporationEffects( G4double kineticEnergy )
{
  // derived from original FORTRAN code EXNU by H. Fesefeldt (10-Dec-1986)
  //
  // Nuclear evaporation as function of atomic number
  // and kinetic energy (MeV) of primary particle
  //
  // returns kinetic energy (MeV)
  //
  if( aEff < 1.5 )
  {
    pnBlackTrackEnergy = dtaBlackTrackEnergy = 0.0;
    return 0.0;
  }
  G4double ek = kineticEnergy/GeV;
  G4float ekin = std::min( 4.0, std::max( 0.1, ek ) );
  const G4float atno = std::min( 120., aEff ); 
  const G4float gfa = 2.0*((aEff-1.0)/70.)*G4Exp(-(aEff-1.0)/70.);
  //
  // 0.35 value at 1 GeV
  // 0.05 value at 0.1 GeV
  //
  G4float cfa = std::max( 0.15, 0.35 + ((0.35-0.05)/2.3)*G4Log(ekin) );
  G4float exnu = 7.716 * cfa * G4Exp(-cfa)
    * ((atno-1.0)/120.)*G4Exp(-(atno-1.0)/120.);
  G4float fpdiv = std::max( 0.5, 1.0-0.25*ekin*ekin );
  //
  // pnBlackTrackEnergy  is the kinetic energy (in GeV) available for
  //                     proton/neutron black track particles
  // dtaBlackTrackEnergy is the kinetic energy (in GeV) available for
  //                     deuteron/triton/alpha black track particles
  //
  pnBlackTrackEnergy = exnu*fpdiv;
  dtaBlackTrackEnergy = exnu*(1.0-fpdiv);
  
  if( G4int(zEff+0.1) != 82 )
  { 
    G4double ran1 = -6.0;
    G4double ran2 = -6.0;
    for( G4int i=0; i<12; ++i )
    {
      ran1 += G4UniformRand();
      ran2 += G4UniformRand();
    }
    pnBlackTrackEnergy *= 1.0 + ran1*gfa;
    dtaBlackTrackEnergy *= 1.0 + ran2*gfa;
  }
  pnBlackTrackEnergy = std::max( 0.0, pnBlackTrackEnergy );
  dtaBlackTrackEnergy = std::max( 0.0, dtaBlackTrackEnergy );
  while( pnBlackTrackEnergy+dtaBlackTrackEnergy >= ek )  /* Loop checking, 02.11.2015, A.Ribon */
  {
    pnBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
    dtaBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
  }
  //G4cout << "EvaporationEffects "<<kineticEnergy<<" "
  //       <<pnBlackTrackEnergy+dtaBlackTrackEnergy<< G4endl;
  return (pnBlackTrackEnergy+dtaBlackTrackEnergy)*GeV;
}

 
G4double
G4Nucleus::AnnihilationEvaporationEffects(G4double kineticEnergy, G4double ekOrg)
{
  // Nuclear evaporation as a function of atomic number and kinetic 
  // energy (MeV) of primary particle.  Modified for annihilation effects. 
  //
  if( aEff < 1.5 || ekOrg < 0.)
  {
    pnBlackTrackEnergyfromAnnihilation = 0.0;
    dtaBlackTrackEnergyfromAnnihilation = 0.0;
    return 0.0;
  }
  G4double ek = kineticEnergy/GeV;
  G4float ekin = std::min( 4.0, std::max( 0.1, ek ) );
  const G4float atno = std::min( 120., aEff ); 
  const G4float gfa = 2.0*((aEff-1.0)/70.)*G4Exp(-(aEff-1.0)/70.);

  G4float cfa = std::max( 0.15, 0.35 + ((0.35-0.05)/2.3)*G4Log(ekin) );
  G4float exnu = 7.716 * cfa * G4Exp(-cfa)
    * ((atno-1.0)/120.)*G4Exp(-(atno-1.0)/120.);
  G4float fpdiv = std::max( 0.5, 1.0-0.25*ekin*ekin );

  pnBlackTrackEnergyfromAnnihilation = exnu*fpdiv;
  dtaBlackTrackEnergyfromAnnihilation = exnu*(1.0-fpdiv);
  
  G4double ran1 = -6.0;
  G4double ran2 = -6.0;
  for( G4int i=0; i<12; ++i ) {
    ran1 += G4UniformRand();
    ran2 += G4UniformRand();
  }
  pnBlackTrackEnergyfromAnnihilation *= 1.0 + ran1*gfa;
  dtaBlackTrackEnergyfromAnnihilation *= 1.0 + ran2*gfa;

  pnBlackTrackEnergyfromAnnihilation = std::max( 0.0, pnBlackTrackEnergyfromAnnihilation);
  dtaBlackTrackEnergyfromAnnihilation = std::max( 0.0, dtaBlackTrackEnergyfromAnnihilation);
  G4double blackSum = pnBlackTrackEnergyfromAnnihilation+dtaBlackTrackEnergyfromAnnihilation;
  if (blackSum >= ekOrg/GeV) {
    pnBlackTrackEnergyfromAnnihilation *= ekOrg/GeV/blackSum;
    dtaBlackTrackEnergyfromAnnihilation *= ekOrg/GeV/blackSum;
  }

  return (pnBlackTrackEnergyfromAnnihilation+dtaBlackTrackEnergyfromAnnihilation)*GeV;
}

 
G4double 
G4Nucleus::Cinema( G4double kineticEnergy )
{
  // derived from original FORTRAN code CINEMA by H. Fesefeldt (14-Oct-1987)
  //
  // input: kineticEnergy (MeV)
  // returns modified kinetic energy (MeV)
  //
  static const G4double expxu =  82.;           // upper bound for arg. of exp
  static const G4double expxl = -expxu;         // lower bound for arg. of exp
  
  G4double ek = kineticEnergy/GeV;
  G4double ekLog = G4Log( ek );
  G4double aLog = G4Log( aEff );
  G4double em = std::min( 1.0, 0.2390 + 0.0408*aLog*aLog );
  G4double temp1 = -ek * std::min( 0.15, 0.0019*aLog*aLog*aLog );
  G4double temp2 = G4Exp( std::max( expxl, std::min( expxu, -(ekLog-em)*(ekLog-em)*2.0 ) ) );
  G4double result = 0.0;
  if( std::abs( temp1 ) < 1.0 )
  {
    if( temp2 > 1.0e-10 )result = temp1*temp2;
  }
  else result = temp1*temp2;
  if( result < -ek )result = -ek;
  return result*GeV;
}


G4ThreeVector G4Nucleus::GetFermiMomentum()
{
  // chv: .. we assume zero temperature!
  
  // momentum is equally distributed in each phasespace volume dpx, dpy, dpz.
  G4double ranflat1=
    G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);   
  G4double ranflat2=
    G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);   
  G4double ranflat3=
    G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);   
  G4double ranmax = (ranflat1>ranflat2? ranflat1: ranflat2);
  ranmax = (ranmax>ranflat3? ranmax : ranflat3);
  
  // Isotropic momentum distribution
  G4double costheta = 2.*G4UniformRand() - 1.0;
  G4double sintheta = std::sqrt(1.0 - costheta*costheta);
  G4double phi = 2.0*pi*G4UniformRand();
  
  G4double pz=costheta*ranmax;
  G4double px=sintheta*std::cos(phi)*ranmax;
  G4double py=sintheta*std::sin(phi)*ranmax;
  G4ThreeVector p(px,py,pz);
  return p;
}
 

G4ReactionProductVector* G4Nucleus::Fragmentate()
{
  // needs implementation!
  return nullptr;
}
 

void G4Nucleus::AddMomentum(const G4ThreeVector aMomentum)
{
  momentum+=(aMomentum);
}

 
void G4Nucleus::AddExcitationEnergy( G4double anEnergy )
{
  excitationEnergy+=anEnergy;
}

 /* end of file */

