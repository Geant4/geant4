// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              IONtest2.cc
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
#include "G4NuclearAbrasionGeometry.hh"
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
void EvaluateNuclearAbrasionGeometry
  (G4ParticleDefinition *theProjectile, G4Element *theTarget)
{
  G4double AP  = theProjectile->GetBaryonNumber();
  G4double AT  = theTarget->GetN();
  G4WilsonRadius wilsonRadius;
  G4double rP  = wilsonRadius.GetWilsonRadius(AP);
  G4double rT  = wilsonRadius.GetWilsonRadius(AT);
  G4double rPT = rP + rT;

  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"EVALUATING FOR " <<theProjectile->GetParticleName()
         <<" INCIDENT UPON " <<theTarget->GetName()
         <<" (AT = " <<AT
         <<", ZT = " <<theTarget->GetZ()
         <<")"
         <<G4endl;
  G4cout <<"Projectile radius = " <<rP/fermi <<" fm" <<G4endl;
  G4cout <<"Target radius     = " <<rT/fermi <<" fm" <<G4endl;
  G4cout <<G4endl;

  G4cout <<"     radial thd   Impact param              F              P"
         <<"            ExP            ExT"
         <<G4endl;
  G4cout <<"             []           [fm]             []             []"
         <<"          [MeV]          [MeV]"
         <<G4endl;
  G4cout <<"------------------------------------------------------------"
         <<"------------------------------"
         <<G4endl;
  for (G4double thd=0.2; thd<0.5; thd+=0.05)
  { 
    for (G4double r=0.0; r<rPT; r+=rPT/10.0)
    {
      G4NuclearAbrasionGeometry abrasionGeometry(AP, AT, r);
      abrasionGeometry.SetPeripheralThreshold(thd);
      G4double F   = abrasionGeometry.F();
      G4double P   = abrasionGeometry.P();
      G4double ExP = abrasionGeometry.GetExcitationEnergyOfProjectile();
      G4double ExT = abrasionGeometry.GetExcitationEnergyOfTarget();
      G4cout <<std::setw(15) <<abrasionGeometry.GetPeripheralThreshold()
             <<std::setw(15) <<r/fermi
             <<std::setw(15) <<F
             <<std::setw(15) <<P
             <<std::setw(15) <<ExP
             <<std::setw(15) <<ExT
             <<G4endl;
    }
  }
  G4cout <<"------------------------------------------------------------"
         <<"------------------------------"
         <<G4endl;

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
  G4cout <<"COMMENCING IONtest2 ...." <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
//
//
// Define projetiles and targets for which cross-sections are required.
//
  G4ParticleDefinition *alpha    = G4Alpha::AlphaDefinition();
  G4ParticleDefinition *B        =
    G4ParticleTable::GetParticleTable()->GetIon(5, 10, 0.0);
  G4ParticleDefinition *C        =
    G4ParticleTable::GetParticleTable()->GetIon(6, 12, 0.0);
  G4ParticleDefinition *Fe        =
    G4ParticleTable::GetParticleTable()->GetIon(26, 56, 0.0);
  
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
  EvaluateNuclearAbrasionGeometry (alpha, helium4);
  EvaluateNuclearAbrasionGeometry (alpha, lithium7);
  EvaluateNuclearAbrasionGeometry (alpha, carbon);
  EvaluateNuclearAbrasionGeometry (alpha, aluminium);
  EvaluateNuclearAbrasionGeometry (alpha, beryllium);
  EvaluateNuclearAbrasionGeometry (alpha, nitrogen);
  EvaluateNuclearAbrasionGeometry (alpha, iron);
  EvaluateNuclearAbrasionGeometry (alpha, tantalum);
  EvaluateNuclearAbrasionGeometry (alpha, gold);

  EvaluateNuclearAbrasionGeometry (B, helium4);
  EvaluateNuclearAbrasionGeometry (B, lithium7);
  EvaluateNuclearAbrasionGeometry (B, carbon);
  EvaluateNuclearAbrasionGeometry (B, aluminium);
  EvaluateNuclearAbrasionGeometry (B, beryllium);
  EvaluateNuclearAbrasionGeometry (B, nitrogen);
  EvaluateNuclearAbrasionGeometry (B, iron);
  EvaluateNuclearAbrasionGeometry (B, tantalum);
  EvaluateNuclearAbrasionGeometry (B, gold);
  
  EvaluateNuclearAbrasionGeometry (C, helium4);
  EvaluateNuclearAbrasionGeometry (C, lithium7);
  EvaluateNuclearAbrasionGeometry (C, carbon);
  EvaluateNuclearAbrasionGeometry (C, aluminium);
  EvaluateNuclearAbrasionGeometry (C, beryllium);
  EvaluateNuclearAbrasionGeometry (C, nitrogen);
  EvaluateNuclearAbrasionGeometry (C, iron);
  EvaluateNuclearAbrasionGeometry (C, tantalum);
  EvaluateNuclearAbrasionGeometry (C, gold);
  
  EvaluateNuclearAbrasionGeometry (Fe, helium4);
  EvaluateNuclearAbrasionGeometry (Fe, lithium7);
  EvaluateNuclearAbrasionGeometry (Fe, carbon);
  EvaluateNuclearAbrasionGeometry (Fe, aluminium);
  EvaluateNuclearAbrasionGeometry (Fe, beryllium);
  EvaluateNuclearAbrasionGeometry (Fe, nitrogen);
  EvaluateNuclearAbrasionGeometry (Fe, iron);
  EvaluateNuclearAbrasionGeometry (Fe, tantalum);
  EvaluateNuclearAbrasionGeometry (Fe, gold);
    
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
  
  G4cout <<"TEST IONtest2 COMPLETE" <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
  
  return 0;
}
