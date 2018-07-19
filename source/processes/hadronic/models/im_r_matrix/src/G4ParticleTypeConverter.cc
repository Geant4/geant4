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

#include "G4ParticleTypeConverter.hh"
#include "G4HadronicException.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4Lambda.hh"

G4ParticleTypeConverter::G4ParticleTypeConverter()
{
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4Proton::ProtonDefinition(), NUCLEON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4Neutron::NeutronDefinition(),  NUCLEON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4Gamma::GammaDefinition(),  GAMMA) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4PionPlus::PionPlusDefinition(),  PION) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4PionZero::PionZeroDefinition(),  PION) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4PionMinus::PionMinusDefinition(),  PION) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4KaonMinus::KaonMinusDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4KaonPlus::KaonPlusDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4KaonZero::KaonZeroDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4AntiKaonZero::AntiKaonZeroDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4KaonZeroShort::KaonZeroShortDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4KaonZeroLong::KaonZeroLongDefinition(),  KAON) );
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4Lambda::LambdaDefinition(),  Lambda) );

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(221),  ETA) ); // eta
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(113),  RHO) ); // rho0
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(213),  RHO) ); // rho+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(-213),  RHO) ); // rho- ?????
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(223),  omega) ); // omega

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2224),  D1232) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2214),  D1232) ); // D+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2114),  D1232) ); // D0
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1114),  D1232) ); // D-
  
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(31114),  D1600) ); // D-
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32114),  D1600) ); // D0
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32214),  D1600) ); // D+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32224),  D1600) ); // D++
  
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1112),  D1620) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1212),  D1620) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2122),  D1620) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2222),  D1620) ); // D++
  
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(11114),  D1700) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12114),  D1700) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12214),  D1700) ); // D++
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12224),  D1700) ); // D++
  
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(11112),  D1900) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(11212),  D1900) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12122),  D1900) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12222),  D1900) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1116),  D1905) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1216),  D1905) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2126),  D1905) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2226),  D1905) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(21112),  D1910) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(21212),  D1910) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22122),  D1910) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22222),  D1910) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(21114),  D1920) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22114),  D1920) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22214),  D1920) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22224),  D1920) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(11116),  D1930) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(11216),  D1930) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12126),  D1930) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12226),  D1930) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1118),  D1950) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2118),  D1950) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2218),  D1950) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2228),  D1950) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12112),  N1440) ); // N0
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12212),  N1440) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1214),  N1520) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2124),  N1520) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22112),  N1535) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22212),  N1535) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32112),  N1650) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32212),  N1650) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2116),  N1675) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2216),  N1675) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12116),  N1680) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12216),  N1680) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(22124),  N1700) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(21214),  N1700) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(42212),  N1710) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(42112),  N1710) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(32124),  N1720) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(31214),  N1720) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(42124),  N1900) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(41214),  N1900) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12218),  N1990) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(12118),  N1990) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(52214),  N2090) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(52114),  N2090) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(2128),  N2190) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(1218),  N2190) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(100002210),  N2220) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(100002110),  N2220) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(100012210),  N2250) ); // N+
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(100012110),  N2250) ); // N+

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13122),  L1405) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3124),  L1520) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23122),  L1600) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(33122),  L1670) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13124),  L1690) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(43122),  L1800) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(53122),  L1810) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3126),  L1820) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13126),  L1830) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23124),  L1890) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3128),  L2100) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23126),  L2110) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3222),  Sigma) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3212),  Sigma) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3112),  Sigma) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3224),  S1385) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3114),  S1385) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3214),  S1385) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13222),  S1660) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13112),  S1660) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13212),  S1660) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13224),  S1670) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13114),  S1670) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13214),  S1670) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23222),  S1750) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23112),  S1750) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23212),  S1750) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3226),  S1775) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3116),  S1775) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3216),  S1775) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13226),  S1915) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13116),  S1915) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13216),  S1915) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23224),  S1940) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23114),  S1940) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23214),  S1940) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3228),  S2030) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3118),  S2030) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3218),  S2030) ); 

  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3324),  X1530) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(3314),  X1530) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23324),  X1690) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(23314),  X1690) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13324),  X1820) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13314),  X1820) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(33324),  X1950) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(33314),  X1950) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13326),  X2030) ); 
  defMap.push_back(std::pair<const G4ParticleDefinition*, GenericType>(G4ParticleTable::GetParticleTable()->FindParticle(13316),  X2030) ); 


//   uMap[NUCLEON] =   1;
//   uMap[N1440]   =   2;
//   uMap[D1232]   =  12;
//   uMap[D1600]   =  13;
//   uMap[GAMMA]   = 100;
//   uMap[PION]    = 101;
}


G4ParticleTypeConverter::GenericType G4ParticleTypeConverter::GetGenericType(const G4ParticleDefinition* const aParticleDef) const
{
  for(size_t i=0;i<defMap.size(); i++)
  {
    if(defMap[i].first == aParticleDef) return defMap[i].second;
  }
  
//GF  G4cerr << "Unknown Particle : " << aParticleDef->GetParticleName() << G4endl;
  return UNKNOWN;
//  throw G4HadronicException(__FILE__, __LINE__, "G4ParticleTypeConverter: unknown particle type!");
}
 
G4ParticleTypeConverter::GenericType G4ParticleTypeConverter::GetGenericType(const G4KineticTrack& aTrack) const
{ 
  return GetGenericType(aTrack.GetDefinition()); 
}

G4int G4ParticleTypeConverter::GetUrqmdItyp(G4ParticleTypeConverter::GenericType ) const
{
  //if (uMap.find(gType)!=uMap.end())
  // hpw return uMap.operator[](gType); 
  //else
  throw G4HadronicException(__FILE__, __LINE__, "G4ParticleTypeConverter: unknown particle type!");
  return 0;
}

const G4ParticleDefinition* G4ParticleTypeConverter::FindIso3State(const G4ParticleTypeConverter::GenericType gType,
								   const G4int isospin3) const
{
  MapIterator iter;
  for (iter = defMap.begin(); iter!=defMap.end(); ++iter) {
    G4ParticleTypeConverter::GenericType foo = (*iter).second;
    if (gType==foo) {
//      G4cout << "convtype " << ((*iter).first)->GetParticleName() << G4endl;
//      G4cout << "conviso3 " << ((*iter).first)->GetPDGiIsospin3() << G4endl;
      if (((*iter).first)->GetPDGiIsospin3()==isospin3)
	return (*iter).first;
    }
  }
//  G4cout << "FindIso3State: can't find " << static_cast<G4int>(gType) << " with iso3 " << isospin3 << G4endl;
  return nullptr;
}
