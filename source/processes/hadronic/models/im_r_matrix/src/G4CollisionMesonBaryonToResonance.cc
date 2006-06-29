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
// $Id: G4CollisionMesonBaryonToResonance.cc,v 1.3 2006-06-29 20:37:28 gunter Exp $ //

#include "globals.hh"
#include "G4CollisionMesonBaryonToResonance.hh"
#include "G4ConcreteMesonBaryonToResonance.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4XAqmElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4XResonance.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4KaonPlus.hh"


G4CollisionMesonBaryonToResonance::G4CollisionMesonBaryonToResonance()
{ 
  // the particle types used are indicative for the particle class; 
  // i.e. pip stands for pions, proton for nucleon, etc..  
  
  G4ParticleDefinition * aProton = G4Proton::ProtonDefinition();
  G4ParticleDefinition * aPionp = G4PionPlus::PionPlusDefinition();
  
  G4ParticleDefinition * aDeltap = G4ParticleTable::GetParticleTable()->FindParticle(2214); // D+
  G4ParticleDefinition * aD1600 = G4ParticleTable::GetParticleTable()->FindParticle(32214); // D+
  G4ParticleDefinition * aD1620 = G4ParticleTable::GetParticleTable()->FindParticle(2122); // D+
  G4ParticleDefinition * aD1700 = G4ParticleTable::GetParticleTable()->FindParticle(12214); // D+
  G4ParticleDefinition * aD1900 = G4ParticleTable::GetParticleTable()->FindParticle(12122); // D+
  G4ParticleDefinition * aD1905 = G4ParticleTable::GetParticleTable()->FindParticle(2126); // D+
  G4ParticleDefinition * aD1910 = G4ParticleTable::GetParticleTable()->FindParticle(22122); // D+
  G4ParticleDefinition * aD1920 = G4ParticleTable::GetParticleTable()->FindParticle(22214); // D+
  G4ParticleDefinition * aD1930 = G4ParticleTable::GetParticleTable()->FindParticle(12126); // D+
  G4ParticleDefinition * aD1950 = G4ParticleTable::GetParticleTable()->FindParticle(2218); // D+

  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aDeltap, "D1232_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1600, "D1600_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1620, "D1620_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1700, "D1700_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1900, "D1900_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1905, "D1905_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1910, "D1910_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1920, "D1920_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1930, "D1930_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aD1950, "D1950_Npi"));  


  G4ParticleDefinition * aN1440 = G4ParticleTable::GetParticleTable()->FindParticle(12112); // N+
  G4ParticleDefinition * aN1520 = G4ParticleTable::GetParticleTable()->FindParticle(2124); // N+
  G4ParticleDefinition * aN1535 = G4ParticleTable::GetParticleTable()->FindParticle(22212); // N+
  G4ParticleDefinition * aN1650 = G4ParticleTable::GetParticleTable()->FindParticle(32212); // N+
  G4ParticleDefinition * aN1675 = G4ParticleTable::GetParticleTable()->FindParticle(2216); // N+
  G4ParticleDefinition * aN1680 = G4ParticleTable::GetParticleTable()->FindParticle(12216); // N+
  G4ParticleDefinition * aN1700 = G4ParticleTable::GetParticleTable()->FindParticle(22124); // N+
  G4ParticleDefinition * aN1710 = G4ParticleTable::GetParticleTable()->FindParticle(42212); // N+
  G4ParticleDefinition * aN1720 = G4ParticleTable::GetParticleTable()->FindParticle(32124); // N+
  G4ParticleDefinition * aN1900 = G4ParticleTable::GetParticleTable()->FindParticle(42124); // N+
  G4ParticleDefinition * aN1990 = G4ParticleTable::GetParticleTable()->FindParticle(12218); // N+
  G4ParticleDefinition * aN2090 = G4ParticleTable::GetParticleTable()->FindParticle(52214); // N+
  G4ParticleDefinition * aN2190 = G4ParticleTable::GetParticleTable()->FindParticle(2128); // N+
  G4ParticleDefinition * aN2220 = G4ParticleTable::GetParticleTable()->FindParticle(100002210); // N+
  G4ParticleDefinition * aN2250 = G4ParticleTable::GetParticleTable()->FindParticle(100012210); // N+

  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1440, "N1440_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1520, "N1520_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1535, "N1535_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1650, "N1650_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1675, "N1675_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1680, "N1680_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1700, "N1700_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1710, "N1710_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1720, "N1720_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1900, "N1900_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN1990, "N1990_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN2090, "N2090_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN2190, "N2190_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN2220, "N2220_Npi"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aPionp, aN2250, "N2250_Npi"));  

/*
  G4ParticleDefinition * aKaon = G4KaonPlus::KaonPlus();
  
// sig=0  G4ParticleDefinition * aL1405 = G4ParticleTable::GetParticleTable()->FindParticle(13122);
  G4ParticleDefinition * aL1520 = G4ParticleTable::GetParticleTable()->FindParticle(3124);
  G4ParticleDefinition * aL1600 = G4ParticleTable::GetParticleTable()->FindParticle(23122);
  G4ParticleDefinition * aL1670 = G4ParticleTable::GetParticleTable()->FindParticle(33122);
  G4ParticleDefinition * aL1690 = G4ParticleTable::GetParticleTable()->FindParticle(13124);
  G4ParticleDefinition * aL1800 = G4ParticleTable::GetParticleTable()->FindParticle(43122);
  G4ParticleDefinition * aL1810 = G4ParticleTable::GetParticleTable()->FindParticle(53122);
  G4ParticleDefinition * aL1820 = G4ParticleTable::GetParticleTable()->FindParticle(3126);
  G4ParticleDefinition * aL1830 = G4ParticleTable::GetParticleTable()->FindParticle(13126);
  G4ParticleDefinition * aL1890 = G4ParticleTable::GetParticleTable()->FindParticle(23124);
  G4ParticleDefinition * aL2100 = G4ParticleTable::GetParticleTable()->FindParticle(3128);
  G4ParticleDefinition * aL2110 = G4ParticleTable::GetParticleTable()->FindParticle(23126);
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1520, "L1520_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1600, "L1600_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1690, "L1690_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1670, "L1670_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1800, "L1800_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1810, "L1810_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1820, "L1820_NKbar"));
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1830, "L1830_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL1890, "L1890_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL2100, "L2100_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aL2110, "L2110_NKbar"));  

//@  G4ParticleDefinition * aSigma = G4ParticleTable::GetParticleTable()->FindParticle(3222);
//@  G4ParticleDefinition * aS1385 = G4ParticleTable::GetParticleTable()->FindParticle(3224);
  G4ParticleDefinition * aS1660 = G4ParticleTable::GetParticleTable()->FindParticle(13222);
  G4ParticleDefinition * aS1670 = G4ParticleTable::GetParticleTable()->FindParticle(13224);
  G4ParticleDefinition * aS1750 = G4ParticleTable::GetParticleTable()->FindParticle(23222);
  G4ParticleDefinition * aS1775 = G4ParticleTable::GetParticleTable()->FindParticle(3226);
  G4ParticleDefinition * aS1915 = G4ParticleTable::GetParticleTable()->FindParticle(13226);
  G4ParticleDefinition * aS1940 = G4ParticleTable::GetParticleTable()->FindParticle(23224);
  G4ParticleDefinition * aS2030 = G4ParticleTable::GetParticleTable()->FindParticle(3228);
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1660, "S1660_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1670, "S1670_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1750, "S1750_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1775, "S1775_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1915, "S1915_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS1940, "S1940_NKbar"));  
  G4CollisionComposite::AddComponent(new G4ConcreteMesonBaryonToResonance(aProton, aKaon, aS2030, "S2030_NKbar"));  
*/
}
