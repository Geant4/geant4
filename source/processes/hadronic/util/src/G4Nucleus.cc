// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Nucleus.cc,v 1.1 1999-01-07 16:13:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 19-Nov-1996
 // last modified: 27-Mar-1997
 // J.P.Wellisch: 23-Apr-97: minor simplifications
 // modified by J.L.Chuma 24-Jul-97  to set the total momentum in Cinema and
 //                                  EvaporationEffects
 // modified by J.L.Chuma 21-Oct-97  put abs() around the totalE^2-mass^2
 //                                  in calculation of total momentum in
 //                                  Cinema and EvaporationEffects
 // Chr. Volcker, 10-Nov-1997: new methods and class variables.
 // HPW added utilities for low energy neutron transport. (12.04.1998)
 // M.G. Pia, 2 Oct 1998: modified GetFermiMomentum to avoid memory leaks
 
#include "G4Nucleus.hh"
#include "Randomize.hh"
 
 G4ReactionProduct G4Nucleus::GetThermalNucleus(G4double targetMass) const
  {
    G4ReactionProduct theTarget;    
    theTarget.SetMass(targetMass*G4Neutron::Neutron()->GetPDGMass());
    G4double px, py, pz;
    px = GetThermalPz(theTarget.GetMass(), theTemp);
    py = GetThermalPz(theTarget.GetMass(), theTemp);
    pz = GetThermalPz(theTarget.GetMass(), theTemp);
    theTarget.SetMomentum(px, py, pz);
    G4double tMom = sqrt(px*px+py*py+pz*pz);
    G4double tEtot = sqrt((tMom+theTarget.GetMass())*
                          (tMom+theTarget.GetMass())-
                          2.*tMom*theTarget.GetMass());
    if(1-tEtot/theTarget.GetMass()>0.001)
    {
      theTarget.SetTotalEnergy(tEtot);
    }
    else
    {
      theTarget.SetKineticEnergy(tMom*tMom/(2.*theTarget.GetMass()));
    }    
    return theTarget;
  }
 
 void
  G4Nucleus::ChooseParameters( const G4Material *aMaterial )
  {
    G4double random = G4UniformRand();
    G4double sum = 0;
    const G4ElementVector *theElementVector = aMaterial->GetElementVector();
    for( G4int i=0; i<aMaterial->GetNumberOfElements(); ++i )
    {
      sum += aMaterial->GetAtomicNumDensityVector()[i];
      if( sum > random ) {
        aEff = (*theElementVector)(i)->GetA()*mole/g;
        zEff = (*theElementVector)(i)->GetZ();
        break;
      }
    }
  }
 
 void
  G4Nucleus::SetParameters( const G4double A, const G4double Z )
  {
    G4int myZ = G4int(Z + 0.5);
    G4int myA = G4int(A + 0.5);   
    if( myA<1 || myZ<0 || myZ>myA )
      G4Exception("G4Nucleus::SetParameters called with non-physical parameters");
    aEff = A;  // atomic weight
    zEff = Z;  // atomic number
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
    // derived from original FORTRAN code ATOMAS by H. Fesefeldt (2-Dec-1986)
    //
    // Computes atomic mass in MeV
    // units for A example:  A = material->GetA()/(g/mole);
    //
    // Note:  can't just use aEff and zEff since the Nuclear Reaction
    //        function needs to calculate atomic mass for various values of A and Z

    const G4double electron_mass = G4Electron::Electron()->GetPDGMass()/MeV;
    const G4double proton_mass = G4Proton::Proton()->GetPDGMass()/MeV;
    const G4double neutron_mass = G4Neutron::Neutron()->GetPDGMass()/MeV;
    const G4double deuteron_mass = G4Deuteron::Deuteron()->GetPDGMass()/MeV;
    const G4double alpha_mass = G4Alpha::Alpha()->GetPDGMass()/MeV;
    
    G4int myZ = G4int(Z + 0.5);
    G4int myA = G4int(A + 0.5);
    if( myZ < 0 )return 0.0;
    if( myZ > myA )return 0.0;
    if( myA == 1 )
    {
      if( myZ == 0 )return neutron_mass*MeV;
      if( myZ == 1 )return proton_mass*MeV + electron_mass*MeV;   // hydrogen
    }
    else if( myA == 2 && myZ == 1 )
    {
      return deuteron_mass*MeV;
    }
    else if( myA == 4 && myZ == 2 )
    {
      return alpha_mass*MeV;
    }
    //
    // Weitzsaecker's Mass formula
    //
    G4double mass =
      (A-Z)*neutron_mass + Z*proton_mass + Z*electron_mass
      - 15.67*A                                          // nuclear volume
      + 17.23*pow(A,2./3.)                               // surface energy
      + 93.15*pow(A/2.-Z,2.)/A                           // asymmetry
      + 0.6984523*pow(Z,2.)*pow(A,-1./3.);               // coulomb
    G4int ipp = (myA - myZ)%2;            // pairing
    G4int izz = myZ%2;
    if( ipp == izz )mass += (ipp+izz-1) * 12.0 * pow(A,-0.5);
    return mass*MeV;
  }
 
 G4double
  G4Nucleus::GetThermalPz( const G4double mass, const G4double temp ) const
  {
    G4double result = 0.0;
    for( int i=0; i<12 ; ++i )
      result += G4UniformRand() - 0.5;
    result *= sqrt(k_Boltzmann*temp*mass); // Das ist impuls (Pz),
                                           // nichtrelativistische rechnung
                                           // Maxwell verteilung angenommen
    if ( G4UniformRand()<0.5 ) result =-result;
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
    G4float ekin = min( 4.0, max( 0.1, ek ) );
    const G4float atno = min( 120., aEff ); 
    const G4float gfa = 2.0*((aEff-1.0)/70.)*exp(-(aEff-1.0)/70.);
    //
    // 0.35 value at 1 GeV
    // 0.05 value at 0.1 GeV
    //
    G4float cfa = max( 0.15, 0.35 + ((0.35-0.05)/2.3)*log(ekin) );
    G4float exnu = 7.716 * cfa * exp(-cfa)
      * ((atno-1.0)/120.)*exp(-(atno-1.0)/120.);
    G4float fpdiv = max( 0.5, 1.0-0.25*ekin*ekin );
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
      //G4double ran1 = RandGauss::shoot();
      //G4double ran2 = RandGauss::shoot();
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
    pnBlackTrackEnergy = max( 0.0, pnBlackTrackEnergy );
    dtaBlackTrackEnergy = max( 0.0, dtaBlackTrackEnergy );
    while( pnBlackTrackEnergy+dtaBlackTrackEnergy >= ek )
    {
      pnBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
      dtaBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
    }
    return (pnBlackTrackEnergy+dtaBlackTrackEnergy)*GeV;
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
    G4double ekLog = log( ek );
    G4double aLog = log( aEff );
    G4double em = min( 1.0, 0.2390 + 0.0408*aLog*aLog );
    G4double temp1 = -ek * min( 0.15, 0.0019*aLog*aLog*aLog );
    G4double temp2 = exp( max( expxl, min( expxu, -(ekLog-em)*(ekLog-em)*2.0 ) ) );
    G4double result = 0.0;
    if( abs( temp1 ) < 1.0 )
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
    G4double ranflat1=RandFlat::shoot((HepDouble)0.,(HepDouble)fermiMomentum);   
    G4double ranflat2=RandFlat::shoot((HepDouble)0.,(HepDouble)fermiMomentum);   
    G4double ranflat3=RandFlat::shoot((HepDouble)0.,(HepDouble)fermiMomentum);   
    G4double ranmax = (ranflat1>ranflat2? ranflat1: ranflat2);
    ranmax = (ranmax>ranflat3? ranmax : ranflat3);
    
    // - random decay angle
    G4double theta=RandFlat::shoot((HepDouble)0.,(HepDouble)pi);  // isotropic decay angle theta
    G4double phi  =RandFlat::shoot((HepDouble)0.,(HepDouble)2*pi);  // isotropic decay angle phi
    
    // - setup ThreeVector
    G4double pz=cos(theta)*ranmax;
    G4double px=sin(theta)*cos(phi)*ranmax;
    G4double py=sin(theta)*sin(phi)*ranmax;
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

