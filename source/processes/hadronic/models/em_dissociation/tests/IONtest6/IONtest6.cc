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
#include "G4EMDissociationCrossSection.hh"
#include "G4EMDissociationSpectrum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "globals.hh"

#include <iomanip>
///////////////////////////////////////////////////////////////////////////////
//
void CalculateVirtualPhotonSpectrum
  (G4ParticleDefinition *theProjectile, G4Element *theTarget)
{
//
//
// Get details of the projectile and target.
//
  G4double AP  = theProjectile->GetBaryonNumber();
  G4double ZP  = theProjectile->GetPDGCharge();
  G4double MP  = theProjectile->GetPDGMass();
  G4double AT  = theTarget->GetN();
  G4double ZT  = theTarget->GetZ();
//
//
// Instantiate the object used to determine the virtual photon spectrum.
//
  G4EMDissociationSpectrum thePhotonSpectrum;
  
  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"VIRTUAL PHOTON SPECTRUM FOR " <<theProjectile->GetParticleName()
         <<" INCIDENT UPON " <<theTarget->GetName()
         <<" (AT = " <<theTarget->GetN()
         <<", ZT = " <<theTarget->GetZ()
         <<")"
         <<G4endl;
  for (G4double Eg=10.0*MeV; Eg<=30.0*MeV; Eg+=5.0*MeV)
  {
    G4cout <<G4endl;
    G4cout <<"VIRTUAL GAMMAS OF ENERGY " <<Eg/MeV <<" MeV" <<G4endl;
    G4cout <<G4endl;
    
    G4cout <<G4endl;
    G4cout <<"     Energy/nuc           beta           bmin"
           <<"        flux E1        flux E2" <<G4endl;
    G4cout <<"      [MeV/nuc]             []        [fermi]"
           <<"     [/mm2-MeV]     [/mm2-MeV]" <<G4endl;
    G4cout <<"---------------------------------------------"
           <<"------------------------------" <<G4endl;
    G4double mult = pow(10.0,1.0/5.0);
    for (G4double E=1.0*MeV; E<1.0E+06*mult*MeV; E*=mult)
    {
      G4double TotE = MP+E*AP;
      G4double b    = sqrt(1.0 - MP*MP/TotE/TotE);
      G4double bmin = thePhotonSpectrum.GetClosestApproach (AP, ZP, AT, ZT, b);
      G4double sE1  = thePhotonSpectrum.GetGeneralE1Spectrum(Eg, b, bmin);
      G4double sE2  = thePhotonSpectrum.GetGeneralE2Spectrum(Eg, b, bmin);
      G4cout.setf(std::ios::scientific);
      G4cout.precision(5);
      G4cout <<std::setw(15) <<E/MeV;
      G4cout.setf(std::ios::fixed);
      G4cout.precision(4);
      G4cout.setf(std::ios::showpoint);
      G4cout <<std::setw(15) <<b
             <<std::setw(15) <<bmin/fermi;
      G4cout.setf(std::ios::scientific);
      G4cout.precision(5);
      G4cout <<std::setw(15) <<sE1*ZT*ZT
             <<std::setw(15) <<sE2*ZT*ZT
             <<G4endl;
    }
    G4cout <<"---------------------------------------------"
           <<"------------------------------" <<G4endl;
  }
  G4cout <<G4endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////
//
void CalculateCrossSections
  (G4ParticleDefinition *theProjectile, G4Element *theTarget)
{
  G4EMDissociationCrossSection theEMDissociationCrossSection;
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
  G4bool isApplicable = theEMDissociationCrossSection.IsApplicable
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
      G4double xs = theEMDissociationCrossSection.GetCrossSection
               (&dynamicParticle, theTarget, 0.0);
      G4cout.setf(std::ios::scientific);
      G4cout.precision(5);
      G4cout <<std::setw(15) <<E
             <<std::setw(15) <<xs;
      G4cout.setf(std::ios::fixed);
      G4cout.precision(5);
      G4cout.setf(std::ios::showpoint);
      G4cout <<std::setw(15) <<xs/millibarn
             <<G4endl;
    }
    G4cout <<"---------------------------------------------" <<G4endl;
    G4cout <<G4endl;
  }

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
  G4cout <<"COMMENCING IONtest6 ...." <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
//
//
// Define projetiles and targets for which cross-sections are required.
//
  G4ParticleDefinition *B        =
    G4ParticleTable::GetParticleTable()->GetIon(5, 10, 0.0);
  G4ParticleDefinition *C        =
    G4ParticleTable::GetParticleTable()->GetIon(6, 12, 0.0);
  G4ParticleDefinition *Fe        =
    G4ParticleTable::GetParticleTable()->GetIon(26, 56, 0.0);
  
  G4Element *carbon     = new G4Element("Carbon", "C", 6.0, 12.011*g/mole);
  G4Element *aluminium  = new G4Element("Aluminium", "Al", 13.0, 26.98154*g/mole);
  G4Element *beryllium  = new G4Element("Beryllium", "Be", 4.0, 9.01218*g/mole);
  G4Element *nitrogen   = new G4Element("Nitrogen", "N", 7.0, 14.0067*g/mole);
  G4Element *iron       = new G4Element("Iron", "Fe", 26.0, 55.847*g/mole);
  G4Element *tantalum   = new G4Element("Tantalum", "Ta", 73.0, 180.9479*g/mole);
  G4Element *gold       = new G4Element("Gold", "Au", 79.0, 196.9665*g/mole);
//
//
// Now calculate and virtual photon spectrum for different combinations of
// projectile and target over the energy range 1.0 MeV/nuc to 1.0E+06 MeV/nuc.
//  
  CalculateVirtualPhotonSpectrum (B, carbon);
  CalculateVirtualPhotonSpectrum (B, nitrogen);
  CalculateVirtualPhotonSpectrum (B, aluminium);
  CalculateVirtualPhotonSpectrum (B, iron);
  CalculateVirtualPhotonSpectrum (B, tantalum);
  CalculateVirtualPhotonSpectrum (B, gold);
  CalculateVirtualPhotonSpectrum (C, carbon);
  CalculateVirtualPhotonSpectrum (C, nitrogen);
  CalculateVirtualPhotonSpectrum (C, aluminium);
  CalculateVirtualPhotonSpectrum (C, iron);
  CalculateVirtualPhotonSpectrum (C, tantalum);
  CalculateVirtualPhotonSpectrum (C, gold);
  CalculateVirtualPhotonSpectrum (Fe, carbon);
  CalculateVirtualPhotonSpectrum (Fe, nitrogen);
  CalculateVirtualPhotonSpectrum (Fe, aluminium);
  CalculateVirtualPhotonSpectrum (Fe, iron);
  CalculateVirtualPhotonSpectrum (Fe, tantalum);
  CalculateVirtualPhotonSpectrum (Fe, gold);
  
  CalculateCrossSections (B, carbon);
  CalculateCrossSections (B, nitrogen);
  CalculateCrossSections (B, aluminium);
  CalculateCrossSections (B, iron);
  CalculateCrossSections (B, tantalum);
  CalculateCrossSections (B, gold);
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
  
  delete carbon;
  delete aluminium;
  delete beryllium;
  delete nitrogen;
  delete iron;
  delete tantalum;
  delete gold;
  
  G4cout <<"TEST IONtest6 COMPLETE" <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
  
  return 0;
}
