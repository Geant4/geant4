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
// 18-Sep-2003 First version is written by T. Koi
// 10-Nov-2003 Bug fix at Cal. ke_per_n and D T. Koi
// 12-Nov-2003 Add energy check at lower side T. Koi
// 26-Dec-2006 Add isotope dependence D. Wright
// 14-Mar-2011 Moved constructor, destructor and virtual methods to source by V.Ivanchenko
// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element

#include "G4IonsKoxCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"

G4IonsKoxCrossSection::G4IonsKoxCrossSection()
  : G4VCrossSectionDataSet("IonsKox"),
// lowerLimit ( 10*MeV ), 
    r0 ( 1.1*fermi ), rc ( 1.3*fermi )
{}

G4IonsKoxCrossSection::~G4IonsKoxCrossSection()
{}

void
G4IonsKoxCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4IonsKoxCrossSection calculates the total reaction cross\n"
          << "section for nucleus-nucleus scattering using the Kox\n"
          << "parameterization.  It is valid for projectiles and targets\n"
          << "of all Z, at projectile energies up to 10 GeV/n.  If the\n"
          << "projectile energy is less than 10 MeV/n, a zero cross section\n"
          << "is returned.\n";
}

G4bool G4IonsKoxCrossSection::IsElementApplicable(const G4DynamicParticle* aDP, 
						  G4int, const G4Material*)
{
  return (1 <= aDP->GetDefinition()->GetBaryonNumber());
}

G4double 
G4IonsKoxCrossSection::GetElementCrossSection(
  const G4DynamicParticle* aParticle, G4int ZZ, const G4Material*)
{
   G4double xsection = 0.0;

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
   G4int Zp = G4int(aParticle->GetDefinition()->GetPDGCharge() / eplus + 0.5); 
   G4double ke_per_N = aParticle->GetKineticEnergy() / Ap;

   // Apply energy check, if less than lower limit then 0 value is returned
   // if (  ke_per_N < lowerLimit ) return xsection;

   G4int At = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(ZZ));
   G4int Zt = ZZ;  

   G4double one_third = 1.0 / 3.0;

   G4double cubicrAt = G4Pow::GetInstance()->powA ( G4double(At) , G4double(one_third) );  
   G4double cubicrAp = G4Pow::GetInstance()->powA ( G4double(Ap) , G4double(one_third) );  

   // rc divide fermi
   G4double Bc = Zt * Zp / ( (rc/fermi) * (cubicrAp+cubicrAt) );

   G4double targ_mass = G4NucleiProperties::GetNuclearMass(At, Zt);
   G4double proj_mass = aParticle->GetMass();
   G4double proj_momentum = aParticle->GetMomentum().mag(); 

   G4double Ecm = calEcm ( proj_mass , targ_mass , proj_momentum ); 
   if( Ecm <= Bc) return xsection;

   G4double Rvol = r0 * (  cubicrAp + cubicrAt );

//   G4double ke_per_N = aParticle->GetKineticEnergy() / Ap;
   G4double c = calCeValue ( ke_per_N / MeV  );  

   G4double a = 1.85;
   G4double Rsurf = r0 * (a*cubicrAp * cubicrAt/(cubicrAp + cubicrAt) - c); 
   G4double D = 5.0 * ( At - 2 * Zt ) * Zp / ( Ap * At );
   Rsurf = Rsurf + D * fermi;  // multiply D by fermi 

   G4double Rint = Rvol + Rsurf;
   xsection = pi * Rint * Rint * ( 1 - Bc / ( Ecm / MeV ) );
  
   return xsection; 
}

G4double
G4IonsKoxCrossSection::calEcm(G4double mp, G4double mt, G4double Plab)
{
   G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
   G4double Ecm = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
   G4double Pcm = Plab * mt / Ecm;
   G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;
   return KEcm;
}


G4double G4IonsKoxCrossSection::calCeValue(const G4double ke)
{
   // Calculate c value 
   // This value is indepenent from projectile and target particle 
   // ke is projectile kinetic energy per nucleon in the Lab system with MeV unit 
   // fitting function is made by T. Koi 
   // There are no data below 30 MeV/n in Kox et al., 

   G4double Ce; 
   G4double log10_ke = std::log10 ( ke );   
   if (log10_ke > 1.5) 
   {
      Ce = - 10.0 / G4Pow::GetInstance()->powA ( G4double(log10_ke) , G4double(5) ) + 2.0;
   }
   else
   {
      Ce = (-10.0/G4Pow::GetInstance()->powA(G4double(1.5), G4double(5) ) + 2.0) /
           G4Pow::GetInstance()->powA(G4double(1.5), G4double(3)) * G4Pow::GetInstance()->powA(G4double(log10_ke), G4double(3) );

   }
   return Ce;
}

