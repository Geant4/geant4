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
// Implementation of formulas in analogy to NASA technical paper 3621 by 
// Tripathi, et al.
// 
// 26-Dec-2006 Isotope dependence added by D. Wright
// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element
//

#include "G4TripathiCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

G4TripathiCrossSection::G4TripathiCrossSection()
 : G4VCrossSectionDataSet("Tripathi")
{}

G4TripathiCrossSection::~G4TripathiCrossSection() 
{}

G4bool 
G4TripathiCrossSection::IsElementApplicable(const G4DynamicParticle* aPart, 
					    G4int, const G4Material*)
{
  G4bool result = false;
  if ( (aPart->GetDefinition()->GetBaryonNumber()>2.5) &&
       ( aPart->GetKineticEnergy()/aPart->GetDefinition()->GetBaryonNumber()<1*GeV) ) {
    result = true;
  }
  return result;
}

G4double G4TripathiCrossSection::
GetElementCrossSection(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result = 0.;
  G4double targetAtomicNumber = G4NistManager::Instance()->GetAtomicMassAmu(ZZ);
  G4double nTargetProtons = ZZ;
  
  G4double kineticEnergy = aPart->GetKineticEnergy()/MeV;
  G4double nProjProtons = aPart->GetDefinition()->GetPDGCharge();
  G4double projectileAtomicNumber = 
    aPart->GetDefinition()->GetBaryonNumber();

  static const G4double nuleonRadius=1.1E-15;
  static const G4double myNuleonRadius=1.36E-15;
  
  // needs target mass
  G4double targetMass = 
     G4ParticleTable::GetParticleTable()->GetIonTable()
	 ->GetIonMass(G4lrint(nTargetProtons), G4lrint(targetAtomicNumber));
  G4LorentzVector pTarget(0,0,0,targetMass); 
  G4LorentzVector pProjectile(aPart->Get4Momentum());
  pTarget = pTarget+pProjectile;
  G4double E_cm = (pTarget.mag()-targetMass-pProjectile.m())/MeV;
  if(E_cm <= DBL_MIN) { return result; }
  // done
  G4double r_rms_p = 0.6 * myNuleonRadius * 
                                   G4Pow::GetInstance()->powA(projectileAtomicNumber, 1./3.);
  G4double r_rms_t = 0.6 * myNuleonRadius * 
                                   G4Pow::GetInstance()->powA(targetAtomicNumber, 1./3.);
  
  // done
  G4double r_p = 1.29*r_rms_p/nuleonRadius ;
  G4double r_t = 1.29*r_rms_t/nuleonRadius;
  
  // done
  G4double Radius = r_p + r_t + 
           1.2*(G4Pow::GetInstance()->powA(targetAtomicNumber, 1./3.) + 
            G4Pow::GetInstance()->powA(projectileAtomicNumber, 1./3.))/G4Pow::GetInstance()->powA(E_cm, 1./3.);

  //done
  G4double B = 1.44*nProjProtons*nTargetProtons/Radius;
  if(E_cm <= B) return result; 
  // done
  G4double Energy = kineticEnergy/projectileAtomicNumber;

  // done
  //
  // Note that this correction to G4TripathiCrossSection is just to accurately
  // reflect Tripathi's algorithm.  However, if you're using alpha 
  // particles/protons consider using the more accurate 
  // G4TripathiLightCrossSection, which Tripathi developed specifically for 
  // light systems.
  //

  G4double D;
  if (nProjProtons==1 && projectileAtomicNumber==1)
  {
    D = 2.05;
  }
  else if (nProjProtons==2 && projectileAtomicNumber==4)
  {
    D = 2.77-(8.0E-3*targetAtomicNumber)+
          (1.8E-5*targetAtomicNumber*targetAtomicNumber)
                   - 0.8/(1+G4Exp((250.-Energy)/75.));
  }
  else
  {
  //
  // This is the original value used in the G4TripathiCrossSection 
  // implementation, and was used for all projectile/target conditions.  
  // I'm not touching this, although judging from Tripathi's paper, this is 
  // valid for cases where the nucleon density changes little with A.
  // 
    D = 1.75;
  }
  // done
  G4double C_E = D * (1-G4Exp(-Energy/40.)) - 
       0.292*G4Exp(-Energy/792.)*std::cos(0.229*G4Pow::GetInstance()->powA(Energy, 0.453));
  
  // done
  G4double S = G4Pow::GetInstance()->powA(projectileAtomicNumber, 1./3.)*
               G4Pow::GetInstance()->powA(targetAtomicNumber, 1./3.)/
               (G4Pow::GetInstance()->powA(projectileAtomicNumber, 1./3.) + 
               G4Pow::GetInstance()->powA(targetAtomicNumber, 1./3.)); 
  
  // done
  G4double deltaE = 1.85*S + 0.16*S/G4Pow::GetInstance()->powA(E_cm,1./3.) - C_E +
                    0.91*(targetAtomicNumber-2.*nTargetProtons)*nProjProtons/
		    (targetAtomicNumber*projectileAtomicNumber);
  
  // done 
  result = pi * nuleonRadius*nuleonRadius * 
           G4Pow::GetInstance()->powA(( G4Pow::GetInstance()->powA(targetAtomicNumber, 1./3.) + 
	         G4Pow::GetInstance()->powA(projectileAtomicNumber, 1./3.) + deltaE),2.) * 
                 (1-B/E_cm);
  
  if(result < 0.) { result = 0.; }
  return result*m2;

}

