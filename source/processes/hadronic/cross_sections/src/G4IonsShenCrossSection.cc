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
// 12-Nov-2003 Add energy check at lower side T. Koi
// 15-Nov-2006 Above 10GeV/n Cross Section become constant T. Koi (SLAC/SCCS)
// 23-Dec-2006 Isotope dependence adde by D. Wright
// 14-Mar-2011 Moved constructor, destructor and virtual methods to source by V.Ivanchenko
// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element
//

#include "G4IonsShenCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"

G4IonsShenCrossSection::G4IonsShenCrossSection()
  : G4VCrossSectionDataSet("IonsShen"),
    upperLimit( 10*GeV ),
//  lowerLimit( 10*MeV ),
    r0 ( 1.1 )
{}

G4IonsShenCrossSection::~G4IonsShenCrossSection()
{}

void
G4IonsShenCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4IonsShenCrossSection calculates the total reaction cross\n"
          << "section for nucleus-nucleus scattering using the Shen\n"
          << "parameterization.  It is valid for projectiles and targets of\n"
          << "all Z, and projectile energies up to 1 TeV/n.  Above 10 GeV/n"
          << "the cross section is constant.  Below 10 MeV/n zero cross\n"
          << "is returned.\n";
}
   
G4bool G4IonsShenCrossSection::IsElementApplicable(const G4DynamicParticle* aDP, 
						   G4int, const G4Material*)
{
  return (1 <= aDP->GetDefinition()->GetBaryonNumber());
}

G4double 
G4IonsShenCrossSection::GetElementCrossSection(const G4DynamicParticle* aParticle, 
					       G4int Z, 
					       const G4Material*)
{
  G4int A = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));
  return GetIsoCrossSection(aParticle, Z, A);
}

G4double G4IonsShenCrossSection::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
						    G4int Zt, G4int At,  
						    const G4Isotope*,
						    const G4Element*,
						    const G4Material*)

{
   G4double xsection = 0.0;

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
   G4int Zp = G4lrint(aParticle->GetDefinition()->GetPDGCharge()/eplus); 
   G4double ke_per_N = aParticle->GetKineticEnergy() / Ap; 
   if ( ke_per_N > upperLimit ) { ke_per_N = upperLimit; }

   // Apply energy check, if less than lower limit then 0 value is returned
   //if (  ke_per_N < lowerLimit ) { return xsection; }

   G4Pow* g4pow = G4Pow::GetInstance();
  
   G4double cubicrAt = g4pow->Z13(At);
   G4double cubicrAp = g4pow->Z13(Ap);
 
   G4double Rt = 1.12 * cubicrAt - 0.94 * ( 1.0 / cubicrAt );
   G4double Rp = 1.12 * cubicrAp - 0.94 * ( 1.0 / cubicrAp );

   G4double r = Rt + Rp + 3.2;   // in fm
   G4double b = 1.0;   // in MeV/fm
   G4double targ_mass = G4NucleiProperties::GetNuclearMass(At, Zt);

   G4double proj_mass = aParticle->GetMass();
   G4double proj_momentum = aParticle->GetMomentum().mag();

   G4double Ecm = calEcmValue (proj_mass, targ_mass, proj_momentum); 

   G4double B = 1.44 * Zt * Zp / r - b * Rt * Rp / ( Rt + Rp ); 
   if(Ecm <= B) { return xsection; }

   G4double c = calCeValue ( ke_per_N / MeV  );  

   G4double R1 = r0 * (cubicrAt + cubicrAp + 1.85*cubicrAt*cubicrAp/(cubicrAt + cubicrAp) - c); 

   G4double R2 = 1.0 * ( At - 2 * Zt ) * Zp / ( Ap * At );


   G4double R3 = (0.176 / g4pow->A13(Ecm)) * cubicrAt * cubicrAp /(cubicrAt + cubicrAp);

   G4double R = R1 + R2 + R3;

   xsection = 10 * pi * R * R * ( 1 - B / Ecm );   
   xsection = xsection * millibarn;   // mulitply xsection by millibarn

   return xsection; 
}

G4double
G4IonsShenCrossSection::calEcmValue(const G4double mp, const G4double mt,
                                    const G4double Plab)
{
   G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
   G4double Ecm = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
   G4double Pcm = Plab * mt / Ecm;
   G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;
   return KEcm;
}


G4double G4IonsShenCrossSection::calCeValue(const G4double ke)
{
  // Calculate c value 
  // This value is indepenent from projectile and target particle 
  // ke is projectile kinetic energy per nucleon in the Lab system 
  // with MeV unit 
  // fitting function is made by T. Koi 
  // There are no data below 30 MeV/n in Kox et al., 

   G4double Ce; 
   G4double log10_ke = std::log10 ( ke );   
   if (log10_ke > 1.5) 
   {
     Ce = -10.0/std::pow(G4double(log10_ke), G4double(5)) + 2.0;
   }
   else
   {
     Ce = (-10.0/std::pow(G4double(1.5), G4double(5) ) + 2.0) /
         std::pow(G4double(1.5) , G4double(3)) * std::pow(G4double(log10_ke), G4double(3));
   }
   return Ce;
}

