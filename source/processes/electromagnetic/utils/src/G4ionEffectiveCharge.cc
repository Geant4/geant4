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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionEffectiveCharge
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
// 12.09.2004 Set low energy limit to 1 keV (V.Ivanchenko) 
// 25.01.2005 Add protection - min Charge 0.1 eplus (V.Ivanchenko) 
// 28.04.2006 Set upper energy limit to 50 MeV (V.Ivanchenko) 
// 23.05.2006 Set upper energy limit to Z*10 MeV (V.Ivanchenko) 
// 15.08.2006 Add protection for not defined material (V.Ivanchenko) 
// 27-09-2007 Use Fermi energy from material, optimazed formulas (V.Ivanchenko)
//

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionEffectiveCharge.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionEffectiveCharge::G4ionEffectiveCharge()
{
  chargeCorrection = 1.0;
  energyHighLimit  = 20.0*MeV;
  energyLowLimit   = 1.0*keV;
  energyBohr       = 25.*keV;
  massFactor       = amu_c2/(proton_mass_c2*keV);
  minCharge        = 1.0;
  lastPart         = 0;
  lastMat          = 0;
  lastKinEnergy    = 0.0;
  effCharge        = eplus;
  nist = G4NistManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionEffectiveCharge::~G4ionEffectiveCharge()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionEffectiveCharge::EffectiveCharge(const G4ParticleDefinition* p,
                                               const G4Material* material,
					       G4double kineticEnergy)
{
  if(p == lastPart && material == lastMat && kineticEnergy == lastKinEnergy)
    return effCharge;

  lastPart      = p;
  lastMat       = material;
  lastKinEnergy = kineticEnergy;

  G4double mass   = p->GetPDGMass();
  G4double charge = p->GetPDGCharge();
  G4double Zi     = charge/eplus;

  chargeCorrection = 1.0;
  effCharge = charge;

  // The aproximation of ion effective charge from:
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Fast ions or hadrons
  G4double reducedEnergy = kineticEnergy * proton_mass_c2/mass ;

  //G4cout << "e= " << reducedEnergy << " Zi= " << Zi << "  " << material->GetName() << G4endl;

  if( reducedEnergy > Zi*energyHighLimit || Zi < 1.5 || !material) return charge;

  G4double z    = material->GetIonisation()->GetZeffective();
  reducedEnergy = std::max(reducedEnergy,energyLowLimit);

  // Helium ion case
  if( Zi < 2.5 ) {

    static G4double c[6] = {0.2865,  0.1266, -0.001429,
			    0.02402,-0.01135, 0.001475} ;

    G4double Q = std::max(0.0,std::log(reducedEnergy*massFactor));
    G4double x = c[0];
    G4double y = 1.0;
    for (G4int i=1; i<6; i++) {
      y *= Q;
      x += y * c[i] ;
    }
    G4double ex;
    if(x < 0.2) ex = x * (1 - 0.5*x);
    else        ex = 1. - std::exp(-x);

    G4double tq = 7.6 - Q;
    G4double tq2= tq*tq;
    G4double tt = ( 0.007 + 0.00005 * z );
    if(tq2 < 0.2) tt *= (1.0 - tq2 + 0.5*tq2*tq2);
    else          tt *= std::exp(-tq2);

    effCharge = charge*(1.0 + tt) * std::sqrt(ex);

    // Heavy ion case
  } else {
    
    G4double y;
    //    = nist->GetZ13(z);
    //G4double z23  = y*y;
    G4double zi13 = nist->GetZ13(Zi);
    G4double zi23 = zi13*zi13;
    //    G4double e = std::max(reducedEnergy,energyBohr/z23);
    //G4double e = reducedEnergy;

    // v1 is ion velocity in vF unit
    G4double eF   = material->GetIonisation()->GetFermiEnergy();
    G4double v1sq = reducedEnergy/eF;
    G4double vFsq = eF/energyBohr;
    G4double vF   = std::sqrt(eF/energyBohr);

    // Faster than Fermi velocity
    if ( v1sq > 1.0 ) {
      y = vF * std::sqrt(v1sq) * ( 1.0 + 0.2/v1sq ) / zi23 ;

      // Slower than Fermi velocity
    } else {
      y = 0.692308 * vF * (1.0 + 0.666666*v1sq + v1sq*v1sq/15.0) / zi23 ;
    }

    G4double q;
    G4double y3 = std::pow(y, 0.3) ;
    // G4cout<<"y= "<<y<<" y3= "<<y3<<" v1= "<<v1<<" vF= "<<vF<<G4endl; 
    q = 1.0 - std::exp( 0.803*y3 - 1.3167*y3*y3 - 0.38157*y - 0.008983*y*y ) ;
    
    //y *= 0.77;
    //y *= (0.75 + 0.52/Zi);

    //if( y < 0.2 ) q = y*(1.0 - 0.5*y);
    //else          q = 1.0 - std::exp(-y);

    G4double qmin = minCharge/Zi;
    if(q < qmin) q = qmin;
  
    effCharge = q*charge;

    /*
    G4double x1 = 1.0*effCharge*(1.0 - 0.132*std::log(y))/(y*std::sqrt(z));
    G4double x2 = 0.1*effCharge*effCharge*energyBohr/reducedEnergy;

    chargeCorrection = 1.0 + x1 - x2;

    G4cout << "x1= "<<x1<<" x2= "<< x2<<" corr= "<<chargeCorrection<<G4endl;
    */
    
    G4double tq = 7.6 - std::log(reducedEnergy/keV);
    G4double tq2= tq*tq;
    G4double sq = ( 0.18 + 0.0015 * z ) / (Zi*Zi);
    if(tq2 < 0.2) sq *= (1.0 - tq2 + 0.5*tq2*tq2);
    else          sq *= std::exp(-tq2);
    sq += 1.0;
    //    G4cout << "sq= " << sq << G4endl;

    // Screen length according to
    // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
    // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

    G4double lambda = 10.0 * vF / (zi13 * (6.0 + q));
    if(q < 0.2) lambda *= (1.0 - 0.66666667*q - q*q/9.0);
    else        lambda *= std::pow(1.0-q, 0.666666);

    G4double lambda2 = lambda*lambda;

    G4double xx = (0.5/q - 0.5)/vFsq;
    if(lambda2 < 0.2) xx *= lambda2*(1.0 - 0.5*lambda2);
    else              xx *= std::log(1.0 + lambda2); 

    chargeCorrection = sq * (1.0 + xx);
    
  }
  //  G4cout << "G4ionEffectiveCharge: charge= " << charge << " q= " << q 
  //         << " chargeCor= " << chargeCorrection 
  //	   << " e(MeV)= " << kineticEnergy/MeV << G4endl;
  return effCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


