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
// 23-Dec-2006 Isotope dependence added by D. Wright
// 14-Mar-2011 Moved constructor, destructor and virtual methods to source by V.Ivanchenko
// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element
//

#include "G4IonsSihverCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"

G4IonsSihverCrossSection::G4IonsSihverCrossSection()
  : G4VCrossSectionDataSet("IonsSihver"), square_r0 ( (1.36*fermi) * (1.36*fermi) )
{}

G4IonsSihverCrossSection::~G4IonsSihverCrossSection()
{}

void
G4IonsSihverCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4IonsSihverCrossSection calculates the total reaction cross\n"
          << "section for nucleus-nucleus scattering using the Sihver\n"
          << "parameterization.  It is valid for projectiles and targets of\n"
          << "all Z, and for all projectile energies above 100 MeV/n.\n"; 
}
   
G4bool 
G4IonsSihverCrossSection::IsElementApplicable(const G4DynamicParticle* aDP, 
					      G4int, const G4Material*)
{
  G4int BaryonNumber = aDP->GetDefinition()->GetBaryonNumber();
  G4double KineticEnergy = aDP->GetKineticEnergy(); 
  if ( KineticEnergy / BaryonNumber >= 100*MeV && BaryonNumber > 1 ) { return true; }
  return false;
}

G4double 
G4IonsSihverCrossSection::GetElementCrossSection(
      const G4DynamicParticle* aParticle, G4int Z, const G4Material*)
{
  G4double xsection = 0.0;
  G4int At = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));

  G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
 
  G4Pow* g4pow = G4Pow::GetInstance();
  
  G4double cubicrAt = g4pow->Z13(At);
  G4double cubicrAp = g4pow->Z13(Ap);
 
  G4double b0 = 1.581 - 0.876 * (1.0/cubicrAp + 1.0/cubicrAt);

  G4double xr = cubicrAp + cubicrAt - b0 * (1.0/cubicrAp + 1.0/cubicrAt);
  xsection = pi * square_r0 * xr * xr;
  
  return xsection; 
}

