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
//

#include "G4TripathiCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"

G4TripathiCrossSection::G4TripathiCrossSection()
 : G4VCrossSectionDataSet("G4TripathiCrossSection")
{
  //  G4cout <<"New G4TripathiCrossSection " << this << G4endl;
}

G4TripathiCrossSection::~G4TripathiCrossSection() 
{}

G4double G4TripathiCrossSection::
GetZandACrossSection(const G4DynamicParticle* aPart, G4int ZZ, G4int AA, 
                     G4double /*temperature*/) 
{
  G4double result = 0.;
  const G4double targetAtomicNumber = AA;
  const G4double nTargetProtons = ZZ;
  
  const G4double kineticEnergy = aPart->GetKineticEnergy()/MeV;
  const G4double nProjProtons = aPart->GetDefinition()->GetPDGCharge();
  const G4double projectileAtomicNumber = 
                             aPart->GetDefinition()->GetBaryonNumber();

  const G4double nuleonRadius=1.1E-15;
  const G4double myNuleonRadius=1.36E-15;
  
  // needs target mass
  G4double targetMass = 
     G4ParticleTable::GetParticleTable()->GetIonTable()
	 ->GetIonMass(G4lrint(nTargetProtons), G4lrint(targetAtomicNumber));
  G4LorentzVector pTarget(0,0,0,targetMass); 
  G4LorentzVector pProjectile(aPart->Get4Momentum());
  pTarget = pTarget+pProjectile;
  G4double E_cm = (pTarget.mag()-targetMass-pProjectile.m())/MeV;
  if(E_cm <= DBL_MIN) return result;  
  // done
  G4double r_rms_p = 0.6 * myNuleonRadius * 
                                   std::pow(projectileAtomicNumber, 1./3.);
  G4double r_rms_t = 0.6 * myNuleonRadius * 
                                   std::pow(targetAtomicNumber, 1./3.);
  
  // done
  G4double r_p = 1.29*r_rms_p/nuleonRadius ;
  G4double r_t = 1.29*r_rms_t/nuleonRadius;
  
  // done
  G4double Radius = r_p + r_t + 
           1.2*(std::pow(targetAtomicNumber, 1./3.) + 
            std::pow(projectileAtomicNumber, 1./3.))/std::pow(E_cm, 1./3.);

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
                   - 0.8/(1+std::exp((250.-Energy)/75.));
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
  G4double C_E = D * (1-std::exp(-Energy/40.)) - 
       0.292*std::exp(-Energy/792.)*std::cos(0.229*std::pow(Energy, 0.453));
  
  // done
  G4double S = std::pow(projectileAtomicNumber, 1./3.)*
               std::pow(targetAtomicNumber, 1./3.)/
               (std::pow(projectileAtomicNumber, 1./3.) + 
               std::pow(targetAtomicNumber, 1./3.)); 
  
  // done
  G4double deltaE = 1.85*S + 0.16*S/std::pow(E_cm,1./3.) - C_E +
                    0.91*(targetAtomicNumber-2.*nTargetProtons)*nProjProtons/
		    (targetAtomicNumber*projectileAtomicNumber);
  
  // done 
  result = pi * nuleonRadius*nuleonRadius * 
           std::pow(( std::pow(targetAtomicNumber, 1./3.) + 
	         std::pow(projectileAtomicNumber, 1./3.) + deltaE),2.) * 
                 (1-B/E_cm);
  
  if(result < 0.) result = 0.;
  return result*m2;

}


G4double G4TripathiCrossSection::
GetCrossSection(const G4DynamicParticle* aPart, const G4Element* anEle, 
    G4double temperature)
{
  G4int nIso = anEle->GetNumberOfIsotopes();
  G4double xsection = 0;
     
  if (nIso) {
    G4double sig;
    G4IsotopeVector* isoVector = anEle->GetIsotopeVector();
    G4double* abundVector = anEle->GetRelativeAbundanceVector();
    G4int ZZ;
    G4int AA;
     
    for (G4int i = 0; i < nIso; i++) {
      ZZ = (*isoVector)[i]->GetZ();
      AA = (*isoVector)[i]->GetN();
      sig = GetZandACrossSection(aPart, ZZ, AA, temperature);
      xsection += sig*abundVector[i];
    }
   
  } else {
    G4int ZZ = G4lrint(anEle->GetZ());
    G4int AA = G4lrint(anEle->GetN());
    xsection = GetZandACrossSection(aPart, ZZ, AA, temperature);
  }

  return xsection;
}
