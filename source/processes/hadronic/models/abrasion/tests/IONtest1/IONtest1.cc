// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TripathiLightCrossSection.cc
//
// Version:		
// Date:		
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4WilsonRadius.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "globals.hh"

#include <iomanip>
///////////////////////////////////////////////////////////////////////////////
//
void CalculateCrossSections
  (G4ParticleDefinition *theProjectile, G4Element *theTarget)
{
  G4TripathiLightCrossSection theTripathiLightCrossSection;
  G4double A = theProjectile->GetBaryonNumber();
  G4DynamicParticle dynamicParticle =
    G4DynamicParticle (theProjectile,G4ThreeVector(0.0,0.0,1.0),A*MeV);
  G4double mult = pow(10.0,1.0/10.0);
  
  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"CROSS-SECTIONS FOR " <<theProjectile->GetParticleName()
         <<" INCIDENT UPON " <<theTarget->GetName()
         <<" (AT = " <<theTarget->GetN()
         <<", ZT = " <<theTarget->GetZ()
         <<")"
         <<G4endl;
  G4bool isApplicable = theTripathiLightCrossSection.IsApplicable
           (&dynamicParticle, theTarget);
  G4cout <<"IsApplicable RETURNS : " <<isApplicable <<G4endl;
  if (isApplicable)
  {
    G4cout <<"     Energy/nuc  Cross-section  Cross-section" <<G4endl;
    G4cout <<"      [MeV/nuc]          [mm2]        [mbarn]" <<G4endl;
    G4cout <<"---------------------------------------------" <<G4endl;
    for (G4double E=1.0; E<1.0E+06; E*=mult)
    {
      dynamicParticle.SetKineticEnergy(A*E);
      G4double xs = theTripathiLightCrossSection.GetCrossSection
               (&dynamicParticle, theTarget, 0.0);
      G4cout <<std::setw(15) <<E
             <<std::setw(15) <<xs
             <<std::setw(15) <<xs/millibarn
             <<G4endl;
    }
    G4cout <<"---------------------------------------------" <<G4endl;
    G4cout <<G4endl;
  }

  return;
}
///////////////////////////////////////////////////////////////////////////////
//
void CheckCrossSectionSource
  (G4ParticleDefinition *theProjectile, G4Element *theTarget)
{
  G4GeneralSpaceNNCrossSection generalCrossSection;
  generalCrossSection.SetVerboseLevel(2);
  G4double crossSection = 0.0;
  G4double E            = 100.0 * MeV * theProjectile->GetBaryonNumber();
  G4DynamicParticle dynamicParticle =
    G4DynamicParticle(theProjectile,G4ThreeVector(0.0,0.0,1.0),E);
  crossSection =
    generalCrossSection.GetCrossSection(&dynamicParticle, theTarget, 0.0);
  E               = 1100.0 * MeV * theProjectile->GetBaryonNumber();
  dynamicParticle =
    G4DynamicParticle(theProjectile,G4ThreeVector(0.0,0.0,1.0),E);
  crossSection =
    generalCrossSection.GetCrossSection(&dynamicParticle, theTarget, 0.0);

  return;
}
///////////////////////////////////////////////////////////////////////////////
//
int main()
{
//
//
// Provide simple banner at start of output.
//
  G4cout <<"COMMENCING IONtest1 ...." <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
//
//
// Check for correct generation of nuclear radii.
//
  G4WilsonRadius theWilsonRadius;
  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"CHECKING G4WilsonRadius (radii are in millimetres)" <<G4endl;
  G4cout <<" Nucleon number     RMS radius NUCFRG2 radius" <<G4endl;
  G4cout <<"---------------------------------------------" <<G4endl;
  for (G4double A=1.0; A<241.0; A+=1.0)
  {
    G4cout <<std::setw(15) <<A 
           <<std::setw(15) <<theWilsonRadius.GetWilsonRMSRadius(A)
           <<std::setw(15) <<theWilsonRadius.GetWilsonRadius(A)
           <<G4endl;
  }
  G4cout <<"---------------------------------------------" <<G4endl;
  G4cout <<G4endl;
//
//
// Define projetiles and targets for which cross-sections are required.
//
  G4ParticleDefinition *neutron  = G4Neutron::NeutronDefinition();
  G4ParticleDefinition *proton   = G4Proton::ProtonDefinition();
  G4ParticleDefinition *deuteron = G4Deuteron::DeuteronDefinition();
  G4ParticleDefinition *he3      = G4He3::He3Definition();
  G4ParticleDefinition *alpha    = G4Alpha::AlphaDefinition();
  G4ParticleDefinition *B        =
    G4ParticleTable::GetParticleTable()->GetIon(5, 10, 0.0);
  G4ParticleDefinition *C        =
    G4ParticleTable::GetParticleTable()->GetIon(6, 12, 0.0);
  G4ParticleDefinition *Fe        =
    G4ParticleTable::GetParticleTable()->GetIon(26, 56, 0.0);
  
  G4Isotope *hydrogenI  = new G4Isotope("Hydrogen", 1, 1, 1.0);
  G4Element *hydrogen   = new G4Element("Hydrogen", "H", 1);
  hydrogen->AddIsotope(hydrogenI, 1.0);
  G4Isotope *deuteriumI = new G4Isotope("Deuterium", 1, 2, 2.0);
  G4Element *deuterium  = new G4Element("Deuterium", "H", 1);
  deuterium->AddIsotope(deuteriumI, 1.0);
  G4Isotope *helium4I   = new G4Isotope("Helium4", 2, 4, 4.0);
  G4Element *helium4    = new G4Element("Helium4", "He", 1);
  helium4->AddIsotope(helium4I, 1.0);
  G4Isotope *helium3I   = new G4Isotope("Helium3", 2, 3, 3.0);
  G4Element *helium3    = new G4Element("Helium3", "He", 1);
  helium3->AddIsotope(helium3I, 1.0);
  G4Isotope *lithium6I  = new G4Isotope("Lithium6", 3, 6, 6.0);
  G4Element *lithium6   = new G4Element("Lithium6", "Li", 1);
  lithium6->AddIsotope (lithium6I, 1.0);
  G4Isotope *lithium7I  = new G4Isotope("Lithium7", 3, 7, 7.0);
  G4Element *lithium7   = new G4Element("Lithium7", "Li", 1);
  lithium7->AddIsotope (lithium7I, 1.0);
  G4Element *carbon     = new G4Element("Carbon", "C", 6.0, 12.011*g/mole);
  G4Element *aluminium  = new G4Element("Aluminium", "Al", 13.0, 26.98154*g/mole);
  G4Element *beryllium  = new G4Element("Beryllium", "Be", 4.0, 9.01218*g/mole);
  G4Element *nitrogen   = new G4Element("Nitrogen", "N", 7.0, 14.0067*g/mole);
  G4Element *iron       = new G4Element("Iron", "Fe", 26.0, 55.847*g/mole);
  G4Element *tantalum   = new G4Element("Tantalum", "Ta", 73.0, 180.9479*g/mole);
  G4Element *gold       = new G4Element("Gold", "Au", 79.0, 196.9665*g/mole);
//
//
// Now calculate and displace cross-sections for different combinations of
// projectile and target over the energy range 1.0 MeV/nuc to 1.0E+06 MeV/nuc.
//  
  CalculateCrossSections (neutron, deuterium);
  CalculateCrossSections (neutron, helium4);
  CalculateCrossSections (proton, deuterium);
  CalculateCrossSections (proton, helium3);
  CalculateCrossSections (proton, helium4);
  CalculateCrossSections (proton, lithium6);
  CalculateCrossSections (proton, lithium7);
  CalculateCrossSections (deuteron, deuterium);
  CalculateCrossSections (deuteron, helium4);
  CalculateCrossSections (deuteron, carbon);
  CalculateCrossSections (he3, beryllium);
  CalculateCrossSections (he3, carbon);
  CalculateCrossSections (alpha, helium4);
  CalculateCrossSections (alpha, beryllium);
  CalculateCrossSections (alpha, nitrogen);
  CalculateCrossSections (alpha, aluminium);
  CalculateCrossSections (alpha, iron);
  CalculateCrossSections (alpha, tantalum);
  CalculateCrossSections (alpha, gold);
  CalculateCrossSections (C, helium3);
  CalculateCrossSections (Fe, helium4);
  CalculateCrossSections (C, carbon);
  CalculateCrossSections (C, nitrogen);
  CalculateCrossSections (C, aluminium);
  CalculateCrossSections (C, iron);
  CalculateCrossSections (C, tantalum);
  CalculateCrossSections (C, gold);
  CalculateCrossSections (Fe, carbon);
  CalculateCrossSections (Fe, nitrogen);
  CalculateCrossSections (Fe, aluminium);
  CalculateCrossSections (Fe, iron);
  CalculateCrossSections (Fe, tantalum);
  CalculateCrossSections (Fe, gold);
//
//
// Declare the G4GeneralSpaceNNCrossSection object and check proper selection of
// model.
//
  G4cout <<"CHECKING FOR PERFORMANCE OF G4GeneralSpaceNNCrossSection" <<G4endl;
  G4cout <<"------------------------------------------------------" <<G4endl;
  CheckCrossSectionSource(proton, hydrogen);
  CheckCrossSectionSource(proton, deuterium);
  CheckCrossSectionSource(proton, helium4);
  CheckCrossSectionSource(proton, lithium7);
  CheckCrossSectionSource(proton, beryllium);
  CheckCrossSectionSource(proton, carbon);
  CheckCrossSectionSource(proton, nitrogen);
  CheckCrossSectionSource(proton, iron);

  CheckCrossSectionSource(he3, hydrogen);
  CheckCrossSectionSource(alpha, hydrogen);
  CheckCrossSectionSource(B, hydrogen);
  CheckCrossSectionSource(C, hydrogen);
  CheckCrossSectionSource(Fe, hydrogen);

  CheckCrossSectionSource (proton, deuterium);
  CheckCrossSectionSource (proton, helium3);
  CheckCrossSectionSource (proton, helium4);
  CheckCrossSectionSource (proton, lithium6);
  CheckCrossSectionSource (proton, lithium7);
  CheckCrossSectionSource (deuteron, deuterium);
  CheckCrossSectionSource (deuteron, helium4);
  CheckCrossSectionSource (deuteron, carbon);
  CheckCrossSectionSource (he3, beryllium);
  CheckCrossSectionSource (he3, carbon);
  CheckCrossSectionSource (alpha, helium4);
  CheckCrossSectionSource (alpha, beryllium);
  CheckCrossSectionSource (alpha, nitrogen);
  CheckCrossSectionSource (alpha, aluminium);
  CheckCrossSectionSource (alpha, iron);
  CheckCrossSectionSource (alpha, tantalum);
  CheckCrossSectionSource (alpha, gold);
  CheckCrossSectionSource (C, helium3);
  CheckCrossSectionSource (Fe, helium4);
  CheckCrossSectionSource (C, carbon);
  CheckCrossSectionSource (C, nitrogen);
  CheckCrossSectionSource (C, aluminium);
  CheckCrossSectionSource (C, iron);
  CheckCrossSectionSource (C, tantalum);
  CheckCrossSectionSource (C, gold);
  CheckCrossSectionSource (Fe, carbon);
  CheckCrossSectionSource (Fe, nitrogen);
  CheckCrossSectionSource (Fe, aluminium);
  CheckCrossSectionSource (Fe, iron);
  CheckCrossSectionSource (Fe, tantalum);
  CheckCrossSectionSource (Fe, gold);
  G4cout <<"------------------------------------------------------" <<G4endl;
  
  delete hydrogen;
  delete deuterium;
  delete helium4;
  delete helium3;
  delete lithium6;
  delete lithium7;
  delete carbon;
  delete aluminium;
  delete beryllium;
  delete nitrogen;
  delete iron;
  delete tantalum;
  delete gold;
  
  G4cout <<"TEST IONtest1 COMPLETE" <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
  
  return 0;
}
