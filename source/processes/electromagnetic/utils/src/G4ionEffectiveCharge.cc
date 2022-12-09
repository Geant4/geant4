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
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionEffectiveCharge::G4ionEffectiveCharge()
{
  chargeCorrection = 1.0;
  energyHighLimit  = 20.0*CLHEP::MeV;
  energyLowLimit   = 1.0*CLHEP::keV;
  energyBohr       = 25.*CLHEP::keV;
  massFactor       = CLHEP::amu_c2/(CLHEP::proton_mass_c2*CLHEP::keV);
  minCharge        = 1.0;
  lastKinEnergy    = 0.0;
  effCharge        = CLHEP::eplus;
  inveplus         = 1.0/CLHEP::eplus;
  g4calc = G4Pow::GetInstance();
}

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

  G4double mass = p->GetPDGMass();
  effCharge = p->GetPDGCharge();
  G4int Zi = G4lrint(effCharge*inveplus);
  chargeCorrection = 1.0;
  if(Zi <= 1) { return effCharge; }

  // The aproximation of ion effective charge from:
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Fast ions or hadrons
  G4double reducedEnergy = kineticEnergy * CLHEP::proton_mass_c2/mass;

  //G4cout << "e= " << reducedEnergy << " Zi= " << Zi << "  " 
  //<< material->GetName() << G4endl;

  if(reducedEnergy > effCharge*energyHighLimit ) {
    return effCharge;
  }
  G4double z = material->GetIonisation()->GetZeffective();
  reducedEnergy = std::max(reducedEnergy,energyLowLimit);

  // Helium ion case
  if( Zi <= 2 ) {

    static const G4double c[6] = 
      {0.2865,0.1266,-0.001429,0.02402,-0.01135,0.001475};

    G4double Q = std::max(0.0,G4Log(reducedEnergy*massFactor));
    G4double x = c[0];
    G4double y = 1.0;
    for (G4int i=1; i<6; ++i) {
      y *= Q;
      x += y * c[i] ;
    }
    G4double ex = (x < 0.2) ? x * (1 - 0.5*x) : 1. - G4Exp(-x);

    G4double tq = 7.6 - Q;
    G4double tq2= tq*tq;
    G4double tt = ( 0.007 + 0.00005 * z );
    if(tq2 < 0.2) { tt *= (1.0 - tq2 + 0.5*tq2*tq2); }
    else          { tt *= G4Exp(-tq2); }

    effCharge *= (1.0 + tt) * std::sqrt(ex);

    // Heavy ion case
  } else {
    
    G4double zi13 = g4calc->Z13(Zi);
    G4double zi23 = zi13*zi13;

    // v1 is ion velocity in vF unit
    G4double eF   = material->GetIonisation()->GetFermiEnergy();
    G4double v1sq = reducedEnergy/eF;
    G4double vFsq = eF/energyBohr;
    G4double vF   = std::sqrt(eF/energyBohr);

    G4double y = ( v1sq > 1.0 ) 
      // Faster than Fermi velocity
      ? vF * std::sqrt(v1sq) * ( 1.0 + 0.2/v1sq ) / zi23
      // Slower than Fermi velocity
      : 0.692308 * vF * (1.0 + 0.666666*v1sq + v1sq*v1sq/15.0) / zi23;

    G4double y3 = G4Exp(0.3*G4Log(y));
    // G4cout<<"y= "<<y<<" y3= "<<y3<<" v1= "<<v1<<" vF= "<<vF<<G4endl; 
    G4double q = std::max(1.0 - G4Exp( 0.803*y3 - 1.3167*y3*y3 - 0.38157*y
                                     - 0.008983*y*y), minCharge/effCharge);
    
    // compute charge correction
    G4double tq = 7.6 - G4Log(reducedEnergy/CLHEP::keV);
    G4double tq2= tq*tq;
    G4double sq = 1.0 + ( 0.18 + 0.0015 * z )*G4Exp(-tq2)/ (Zi*Zi); 
    //    G4cout << "sq= " << sq << G4endl;

    // Screen length according to
    // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
    // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

    G4double lambda = 10.0 * vF *g4calc->A23(1.0 - q)/ (zi13 * (6.0 + q));
    G4double lambda2 = lambda*lambda;
    G4double xx = (0.5/q - 0.5)*G4Log(1.0 + lambda2)/vFsq;

    effCharge *= q;
    chargeCorrection = sq * (1.0 + xx);
  }
  //  G4cout << "G4ionEffectiveCharge: charge= " << charge << " q= " << q 
  //         << " chargeCor= " << chargeCorrection 
  //           << " e(MeV)= " << kineticEnergy/MeV << G4endl;
  return effCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
