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
 // A. Ribon, 6 August 2015: migrated to G4Exp and G4Log.
 
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"

#include "G4Exp.hh"
#include "G4Log.hh"

 
G4Nucleus::G4Nucleus()
  : theA(0), theZ(0), aEff(0.0), zEff(0)
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

G4Nucleus::G4Nucleus( const G4double A, const G4double Z )
{
  SetParameters( A, Z );
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

G4Nucleus::G4Nucleus( const G4int A, const G4int Z )
{
  SetParameters( A, Z );
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

G4ReactionProduct G4Nucleus::
GetBiasedThermalNucleus(G4double aMass, G4ThreeVector aVelocity, G4double temp) const
{
  G4double velMag = aVelocity.mag();
  G4ReactionProduct result;
  G4double value = 0;
  G4double random = 1;
  G4double norm = 3.*std::sqrt(k_Boltzmann*temp*aMass*G4Neutron::Neutron()->GetPDGMass());
  norm /= G4Neutron::Neutron()->GetPDGMass();
  norm *= 5.;
  norm += velMag;
  norm /= velMag;
  const G4int maxNumberOfLoops = 1000000;
  G4int loopCounter = -1;
  while ( (value/norm<random) && ++loopCounter < maxNumberOfLoops )  /* Loop checking, 02.11.2015, A.Ribon */
  {
     result = GetThermalNucleus(aMass, temp);
     G4ThreeVector targetVelocity = 1./result.GetMass()*result.GetMomentum();
     value = (targetVelocity+aVelocity).mag()/velMag;
     random = G4UniformRand();
  }
  if ( loopCounter >= maxNumberOfLoops ) {
    G4ExceptionDescription ed;
    ed << " Failed sampling after maxNumberOfLoops attempts : forced exit! " << G4endl;
    G4Exception( " G4Nucleus::GetBiasedThermalNucleus ", "HAD_NUCLEUS_001", JustWarning, ed );
    result = GetThermalNucleus(aMass, temp);    
  }
  return result;
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
  G4Element* element = (*theElementVector)[aMaterial->GetNumberOfElements()-1];

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
    aEff=theA;
    zEff=theZ;
  } else {   
    aEff = element->GetN();
    zEff = element->GetZ();
    theZ = G4int(zEff + 0.5);
    theA = G4int(aEff + 0.5);   
  }
}


void
G4Nucleus::SetParameters(G4double A, G4double Z)
{
  theZ = G4lrint(Z);
  theA = G4lrint(A);   
  if (theA<1 || theZ<0 || theZ>theA) {
    throw G4HadronicException(__FILE__, __LINE__,
            "G4Nucleus::SetParameters called with non-physical parameters");
  }
  aEff = A;  // atomic weight
  zEff = Z;  // atomic number
  fIsotope = 0;
}

void
G4Nucleus::SetParameters(G4int A, const G4int Z )
{
  theZ = Z;
  theA = A;   
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
    // choose a proton or a neutron as the target particle
    
    G4DynamicParticle *targetParticle = new G4DynamicParticle;
    if( G4UniformRand() < zEff/aEff )
      targetParticle->SetDefinition( G4Proton::Proton() );
    else
      targetParticle->SetDefinition( G4Neutron::Neutron() );
    return targetParticle;
  }
 
 G4double
  G4Nucleus::AtomicMass( const G4double A, const G4double Z ) const
  {
    // Now returns (atomic mass - electron masses) 
    return G4NucleiProperties::GetNuclearMass(A, Z);
  }
 
 G4double
  G4Nucleus::AtomicMass( const G4int A, const G4int Z ) const
  {
    // Now returns (atomic mass - electron masses) 
    return G4NucleiProperties::GetNuclearMass(A, Z);
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
//    G4cout << "EvaporationEffects "<<kineticEnergy<<" "
//           <<pnBlackTrackEnergy+dtaBlackTrackEnergy<<endl;
    return (pnBlackTrackEnergy+dtaBlackTrackEnergy)*GeV;
  }
 
 G4double G4Nucleus::AnnihilationEvaporationEffects(G4double kineticEnergy, G4double ekOrg)
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

 //
 // methods for class G4Nucleus  ... by Christian Volcker
 //

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
    return NULL;
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

