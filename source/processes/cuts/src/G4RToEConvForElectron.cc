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
// $Id: G4RToEConvForElectron.cc,v 1.5 2006/06/29 19:30:22 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
  if( std::abs(AtomicNumber-Z)>0.1 ) {
    Z = AtomicNumber;
    taul = Tlow/Mass;
    ionpot = 1.6e-5*MeV*std::exp(0.9*std::log(Z))/Mass;
    ionpotlog = std::log(ionpot);
  } 


  G4double tau = KineticEnergy/Mass;
  G4double dEdx;

  if(tau<taul) {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double f = 1.-beta2+std::log(tsq/2.)
                  +(0.5+0.25*tsq+(1.+2.*taul)*std::log(0.5))/(t1*t1);
    dEdx = (std::log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*std::sqrt(taul);
    dEdx = clow/std::sqrt(KineticEnergy/Mass);

  } else {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 1.-beta2+std::log(tsq/2.)
                   +(0.5+0.25*tsq+(1.+2.*tau)*std::log(0.5))/(t1*t1);
    dEdx = (std::log(2.*tau+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;

    // loss from bremsstrahlung follows
    G4double cbrem = (cbr1+cbr2*Z)
                       *(cbr3+cbr4*std::log(KineticEnergy/Thigh));
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
  G4double clim = std::sqrt(taulim)*loss;
  G4double taumax = maxEnergy/aMass;

  // now the range vector can be filled
  for ( i=0; i<size_t(TotBin); i++) {
    G4double LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    G4double tau = LowEdgeEnergy/aMass;

    if ( tau <= taulim ) {
      G4double Value = 2.*aMass*tau*std::sqrt(tau)/(3.*clim);
      rangeVector->PutValue(i,Value);
    } else {
      G4double rangelim = 2.*aMass*taulim*std::sqrt(taulim)/(3.*clim);
      G4double ltaulow = std::log(taulim);
      G4double ltauhigh = std::log(tau);
      G4double ltaumax = std::log(taumax);
      G4int    nbin = G4int(maxnbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow));
      if( nbin < 1 ) nbin = 1;
      G4double Value = RangeLogSimpson( NumEl, elementVector,
					atomicNumDensityVector, aMass,
					ltaulow, ltauhigh, nbin) 
	             + rangelim;
      rangeVector->PutValue(i,Value);
    }
  }
} 
