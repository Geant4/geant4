//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RToEConvForElectron.cc,v 1.2 2002-12-16 11:15:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4RToEConvForElectron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "g4std/iomanip"
#include "g4std/strstream"

G4RToEConvForElectron::G4RToEConvForElectron() : G4VRangeToEnergyConverter()
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("e-");
  if (theParticle ==0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4RToEConvForElectron::G4RToEConvForElectron() ";
      G4cout << " Electron is not defined !!" << G4endl;
    }
#endif
  } 
}

G4RToEConvForElectron::~G4RToEConvForElectron()
{ 
}


// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************
G4double G4RToEConvForElectron::ComputeLoss(G4double AtomicNumber,
					    G4double KineticEnergy) const
{
  static G4double Z;  
  static G4double taul, ionpot, ionpotlog;
  const  G4double cbr1=0.02, cbr2=-5.7e-5, cbr3=1., cbr4=0.072;
  const  G4double Tlow=10.*keV, Thigh=1.*GeV;
  static G4double bremfactor= 0.1 ;
 
  G4double Mass = theParticle->GetPDGMass();       
  //  calculate dE/dx for electrons
  if( abs(AtomicNumber-Z)>0.1 ) {
    Z = AtomicNumber;
    taul = Tlow/Mass;
    ionpot = 1.6e-5*MeV*exp(0.9*log(Z))/Mass;
    ionpotlog = log(ionpot);
  } 


  G4double tau = KineticEnergy/Mass;
  G4double dEdx;

  if(tau<taul) {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double f = 1.-beta2+log(tsq/2.)
                  +(0.5+0.25*tsq+(1.+2.*taul)*log(0.5))/(t1*t1);
    dEdx = (log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*sqrt(taul);
    dEdx = clow/sqrt(KineticEnergy/Mass);

  } else {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 1.-beta2+log(tsq/2.)
                   +(0.5+0.25*tsq+(1.+2.*tau)*log(0.5))/(t1*t1);
    dEdx = (log(2.*tau+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;

    // loss from bremsstrahlung follows
    G4double cbrem = (cbr1+cbr2*Z)
                       *(cbr3+cbr4*log(KineticEnergy/Thigh));
    cbrem = Z*(Z+1.)*cbrem*tau/beta2;

    cbrem *= bremfactor ;

    dEdx += twopi_mc2_rcl2*cbrem;
  }

  return dEdx;
}


void G4RToEConvForElectron::BuildRangeVector(const G4Material* aMaterial,
					     G4double       maxEnergy,
					     G4double       aMass,
					     G4PhysicsLogVector* rangeVector)
{
  //  create range vector for a material
  const G4double tlim = 10.*keV;
  const G4int maxnbint = 100;

  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  G4int NumEl = aMaterial->GetNumberOfElements();

  // calculate parameters of the low energy part first
  size_t i;
  G4double loss=0.;
  for (i=0; i<size_t(NumEl); i++) {
    G4bool isOut;
    G4int IndEl = (*elementVector)[i]->GetIndex();
    loss += atomicNumDensityVector[i]*
           (*theLossTable)[IndEl]->GetValue(tlim,isOut);
  }
  G4double taulim = tlim/aMass;
  G4double clim = sqrt(taulim)*loss;
  G4double taumax = maxEnergy/aMass;

  // now the range vector can be filled
  for ( i=0; i<size_t(TotBin); i++) {
    G4double LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    G4double tau = LowEdgeEnergy/aMass;

    if ( tau <= taulim ) {
      G4double Value = 2.*aMass*tau*sqrt(tau)/(3.*clim);
      rangeVector->PutValue(i,Value);
    } else {
      G4double rangelim = 2.*aMass*taulim*sqrt(taulim)/(3.*clim);
      G4double ltaulow = log(taulim);
      G4double ltauhigh = log(tau);
      G4double ltaumax = log(taumax);
      G4int    nbin = G4int(maxnbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow));
      if( nbin < 1 ) nbin = 1;
      G4double Value = RangeLogSimpson(elementVector, atomicNumDensityVector,
                                       aMass,
                                       ltaulow,       ltauhigh) 
	             + rangelim;
      rangeVector->PutValue(i,Value);
    }
  }
} 
