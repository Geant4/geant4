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
//  Calculation of the total, elastic and inelastic cross-sections
//  of anti-nucleon and anti-nucleus interactions with nuclei
//  based on Glauber approach and V. Grishine formulaes for
//  interpolations (ref. V.M.Grichine, Eur.Phys.J., C62(2009) 399;
//  NIM, B267 (2009) 2460) and our parametrization of hadron-nucleon
//  cross-sections
// 
// 
//   Created by A.Galoyan and V. Uzhinsky, 18.11.2010


#include "G4ComponentAntiNuclNuclearXS.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicException.hh"


/////////////////////////////////////////////////////////////////////////////

G4ComponentAntiNuclNuclearXS::G4ComponentAntiNuclNuclearXS() 
: G4VComponentCrossSection("AntiAGlauber"),
  fRadiusEff(0.0),
  fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0),
  fAntiHadronNucleonTotXsc(0.0), fAntiHadronNucleonElXsc(0.0),
  Elab(0.0), S(0.0), SqrtS(0) 
{
  theAProton   = G4AntiProton::AntiProton();
  theANeutron  = G4AntiNeutron::AntiNeutron();
  theADeuteron = G4AntiDeuteron::AntiDeuteron();
  theATriton   = G4AntiTriton::AntiTriton();
  theAAlpha    = G4AntiAlpha::AntiAlpha();
  theAHe3      = G4AntiHe3::AntiHe3();
  Mn     = 0.93827231;           // GeV
  b0     = 11.92;                // GeV^(-2)
  b2     = 0.3036;               // GeV^(-2)
  SqrtS0 = 20.74;                // GeV
  S0     = 33.0625;              // GeV^2
  R0     = 1.0;                  // default value (V.Ivanchenko)
}


/////////////////////////////////////////////////////////////////////////////

G4ComponentAntiNuclNuclearXS::~G4ComponentAntiNuclNuclearXS()
{
}


/////////////////////////////////////////////////////////////////////////////
//
// Calculation of total CrossSection of Anti-Nucleus - Nucleus 

G4double G4ComponentAntiNuclNuclearXS::GetTotalElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double A)
{
  if ( aParticle == nullptr ) { 
    G4ExceptionDescription ed;
    ed << "anti-nucleus with nullptr particle definition: " << aParticle << G4endl; 
    G4Exception( "G4ComponentAntiNuclNuclearXS::GetTotalElementCrossSection", 
                 "antiNuclNuclearXS001", JustWarning, ed );
    return 0.0;
  }
  
  const G4ParticleDefinition* theParticle = aParticle;
  G4double sigmaTotal = GetAntiHadronNucleonTotCrSc(theParticle,kinEnergy);

  // calculation of squared radius of  NN-collision
  G4int i(-1), j(-1);
  if      ( theParticle == theAProton  ||
	    theParticle == theANeutron )  { i=0; } 
  else if ( theParticle == theADeuteron ) { i=1; }
  else if ( theParticle == theATriton   ) { i=2; }
  else if ( theParticle == theAHe3      ) { i=3; }
  else if ( theParticle == theAAlpha    ) { i=4; } 
  else {};

  if ( i < 0  && ( ! theParticle->IsAntiHypernucleus() ) ) { 
    G4ExceptionDescription ed;
    ed << "Unknown anti-nucleus : " << theParticle->GetParticleName() << G4endl
       << "Target (Z, A)=(" << Z << "," << A << ")" << G4endl;
    G4Exception( "G4ComponentAntiNuclNuclearXS::GetTotalElementCrossSection", 
                 "antiNuclNuclearXS002", JustWarning, ed );
  }

  G4int intA = static_cast<G4int>( A );

  if      ( Z == 1  &&  intA == 1 ) { j=0; }
  else if ( Z == 1  &&  intA == 2 ) { j=1; }
  else if ( Z == 1  &&  intA == 3 ) { j=2; }
  else if ( Z == 2  &&  intA == 3 ) { j=3; }
  else if ( Z == 2  &&  intA == 4 ) { j=4; }
  else {}

  if ( i <  0  &&  j >= 0 ) { fRadiusEff = ReffTot[4][j]; }  // Treat all anti-hypernuclei as anti-alpha
  if ( i == 0  &&  j == 0 ) return sigmaTotal * millibarn;   // Pbar/Nbar + P 
  if ( i >= 0  &&  j >= 0 ) { fRadiusEff = ReffTot[i][j]; }  // Light anti-nuclei + Light nuclei

  if ( j < 0 ) {
    if      ( i  == 0 ) { fRadiusEff = 1.34 * theG4Pow->powZ(intA, 0.23)  // Anti-proton/Anti-neutron + Nucleus
                                     + 1.35 / theG4Pow->Z13(intA); } 
    else if ( i  == 1 ) { fRadiusEff = 1.46 * theG4Pow->powZ(intA, 0.21)  // Anti-deuteron + Nucleus
                                     + 1.45 / theG4Pow->Z13(intA); }
    else if ( i  == 2 ) { fRadiusEff = 1.40 * theG4Pow->powZ(intA, 0.21)  // Anti-tritium + Nucleus
                                     + 1.63 / theG4Pow->Z13(intA); }
    else if ( i  == 3 ) { fRadiusEff = 1.40 * theG4Pow->powZ(intA, 0.21)  // Anti-He3 + Nucleus
                                     + 1.63 / theG4Pow->Z13(intA); }
    else if ( i  == 4 ) { fRadiusEff = 1.35 * theG4Pow->powZ(intA, 0.21)  // Anti-alpha + Nucleus
                                     + 1.10 / theG4Pow->Z13(intA); }
    else if ( i  <  0 ) { fRadiusEff = 1.35 * theG4Pow->powZ(intA, 0.21)  // Anti-hypernucleus + Nucleus
	                             + 1.10 / theG4Pow->Z13(intA); }      // is treated as Anti-alpha + Nucleus
    else {}
  }

  G4double R2   = fRadiusEff*fRadiusEff;
  G4double ApAt = std::abs(theParticle->GetBaryonNumber()) * A;

  G4double xsection = millibarn*2.*pi*R2*10.*G4Log(1.+(ApAt*sigmaTotal/(2.*pi*R2*10.)));  //mb
  fTotalXsc = xsection;

  return fTotalXsc; 
}


/////////////////////////////////////////////////////////////////////////////
// 
// Calculation of total CrossSection of Anti-Nucleus - Nucleus 

G4double G4ComponentAntiNuclNuclearXS::GetTotalIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A )
{ 
  return GetTotalElementCrossSection(aParticle, kinEnergy, Z, (G4double) A);
}


/////////////////////////////////////////////////////////////////////////////
// Calculation of inelastic CrossSection of Anti-Nucleus - Nucleus

G4double G4ComponentAntiNuclNuclearXS::GetInelasticElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double A)
{
  if ( aParticle == nullptr ) {
    G4ExceptionDescription ed;
    ed << "anti-nucleus with nullptr particle definition: " << aParticle << G4endl; 
    G4Exception( "G4ComponentAntiNuclNuclearXS::GetInelasticElementCrossSection", 
                 "antiNuclNuclearXS003", JustWarning, ed );
    return 0.0;
  }
  
  const G4ParticleDefinition* theParticle = aParticle;
  G4double sigmaTotal   = GetAntiHadronNucleonTotCrSc(theParticle,kinEnergy);
  G4double sigmaElastic = GetAntiHadronNucleonElCrSc(theParticle,kinEnergy);
  
  // calculation of sqr of radius NN-collision
  G4int i(-1), j(-1);
  if      ( theParticle == theAProton  ||
	    theParticle == theANeutron )  { i=0; } 
  else if ( theParticle == theADeuteron ) { i=1; }
  else if ( theParticle == theATriton   ) { i=2; }
  else if ( theParticle == theAHe3      ) { i=3; }
  else if ( theParticle == theAAlpha    ) { i=4; }
  else {};

  if ( i < 0  && ( ! theParticle->IsAntiHypernucleus() ) ) { 
    G4ExceptionDescription ed;
    ed << "Unknown anti-nucleus : " << theParticle->GetParticleName() << G4endl
       << "Target (Z, A)=(" << Z << "," << A << ")" << G4endl;
    G4Exception( "G4ComponentAntiNuclNuclearXS::GetInelasticElementCrossSection", 
                 "antiNuclNuclearXS004", JustWarning, ed );
  }

  G4int intA = static_cast<G4int>( A );

  if      ( Z == 1  &&  intA == 1 ) { j=0; }
  else if ( Z == 1  &&  intA == 2 ) { j=1; }
  else if ( Z == 1  &&  intA == 3 ) { j=2; }
  else if ( Z == 2  &&  intA == 3 ) { j=3; }
  else if ( Z == 2  &&  intA == 4 ) { j=4; }
  else {}

  if ( i <  0  &&  j >= 0 ) { fRadiusEff = ReffInel[4][j]; }                 // Treat all anti-hypernuclei as anti-alpha
  if ( i == 0  &&  j == 0 ) return (sigmaTotal - sigmaElastic) * millibarn;  // Pbar/Nbar + P 
  if ( i >= 0  &&  j >= 0 ) { fRadiusEff = ReffInel[i][j]; }                 // Light anti-nuclei + Light nuclei

  if ( j < 0) {
    if      ( i  == 0 ) { fRadiusEff = 1.31*theG4Pow->powZ(intA, 0.22)  // Anti-proton/Anti-neutron + Nucleus
                                     + 0.90/theG4Pow->Z13(intA); }
    else if ( i  == 1 ) { fRadiusEff = 1.38*theG4Pow->powZ(intA, 0.21)  // Anti-deuteron + Nucleus
                                     + 1.55/theG4Pow->Z13(intA); }
    else if ( i  == 2 ) { fRadiusEff = 1.34*theG4Pow->powZ(intA, 0.21)  // Anti-tritium + Nucleus
                                     + 1.51/theG4Pow->Z13(intA); }
    else if ( i  == 3 ) { fRadiusEff = 1.34*theG4Pow->powZ(intA, 0.21)  // Anti-He3 + Nucleus
                                     + 1.51/theG4Pow->Z13(intA); }
    else if ( i  == 4 ) { fRadiusEff = 1.30*theG4Pow->powZ(intA, 0.21)  // Anti-alpha + Nucleus
                                     + 1.05/theG4Pow->Z13(intA); }
    else if ( i  <  0 ) { fRadiusEff = 1.30*theG4Pow->powZ(intA,0.21)   // Anti-hypernucleus + Nucleus
                                     + 1.05/theG4Pow->Z13(intA); }      // is treated as Anti-alpha + Nucleus
    else {}
  }

  G4double R2   = fRadiusEff*fRadiusEff;
  G4double ApAt = std::abs(theParticle->GetBaryonNumber()) * A;

  G4double inelxsection = millibarn*pi*R2*10.*G4Log(1.+(ApAt*sigmaTotal/(pi*R2*10.)));  //mb
  fInelasticXsc = inelxsection; 

  return fInelasticXsc;
}


/////////////////////////////////////////////////////////////////////////////
//
// Calculates Inelastic Anti-nucleus-Nucleus cross-section   

G4double G4ComponentAntiNuclNuclearXS::GetInelasticIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A)
{
  return GetInelasticElementCrossSection(aParticle, kinEnergy, Z, (G4double) A);
}


/////////////////////////////////////////////////////////////////////////////
//
// Calculates elastic Anti-nucleus-Nucleus cross-section  as Total - Inelastic 

G4double G4ComponentAntiNuclNuclearXS::GetElasticElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double A)
{
  fElasticXsc = GetTotalElementCrossSection(aParticle, kinEnergy, Z, A)-
                GetInelasticElementCrossSection(aParticle, kinEnergy, Z, A);
  if (fElasticXsc < 0.) fElasticXsc = 0.;
  return fElasticXsc;
}

 
/////////////////////////////////////////////////////////////////////////////
//
// Calculates elastic Anti-nucleus-Nucleus cross-section   

G4double G4ComponentAntiNuclNuclearXS::GetElasticIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A)
{ 
  return GetElasticElementCrossSection(aParticle, kinEnergy, Z, (G4double) A);
}


/////////////////////////////////////////////////////////////////////////////
// Calculation of  Antihadron - hadron Total Cross-section  

G4double G4ComponentAntiNuclNuclearXS::GetAntiHadronNucleonTotCrSc
(const G4ParticleDefinition* aParticle, G4double kinEnergy)
{
  G4double xsection, Pmass, Energy, momentum;
  const G4ParticleDefinition* theParticle = aParticle;
  Pmass=theParticle->GetPDGMass();
  Energy=Pmass+kinEnergy;
  momentum=std::sqrt(Energy*Energy-Pmass*Pmass)/std::abs(theParticle->GetBaryonNumber());
  G4double Plab = momentum / GeV;

  G4double   B, SigAss;
  G4double   C, d1, d2, d3;
  Elab     = std::sqrt(Mn*Mn + Plab*Plab);   // GeV
  S        = 2.*Mn*Mn + 2. *Mn*Elab;         // GeV^2
  SqrtS    = std::sqrt(S);                   // GeV 
  B        = b0+b2*G4Log(SqrtS/SqrtS0)*G4Log(SqrtS/SqrtS0); //GeV^(-2)
  SigAss   = 36.04 +0.304*G4Log(S/S0)*G4Log(S/S0);          //mb 
  R0       = std::sqrt(0.40874044*SigAss - B);                   //GeV^(-2)
  C        = 13.55;
  d1       = -4.47;
  d2       = 12.38;
  d3       = -12.43;

  xsection = SigAss * ( 1 + 1./(std::sqrt(S-4.*Mn*Mn)) / (theG4Pow->powN(R0, 3))
                        * C * ( 1 + d1/SqrtS + d2/(theG4Pow->powN(SqrtS, 2))
                                + d3/(theG4Pow->powN(SqrtS, 3)) ) );

  //xsection *= millibarn;
  fAntiHadronNucleonTotXsc = xsection;

  return fAntiHadronNucleonTotXsc;
}


// //////////////////////////////////////////////////////////////////////////
// Calculation of  Antihadron - hadron Elastic Cross-section  

G4double G4ComponentAntiNuclNuclearXS :: 
GetAntiHadronNucleonElCrSc(const G4ParticleDefinition* aParticle, G4double kinEnergy)
{
  G4double xsection;
  G4double   SigAss;
  G4double   C, d1, d2, d3;
  GetAntiHadronNucleonTotCrSc(aParticle,kinEnergy);
  SigAss   = 4.5 + 0.101*G4Log(S/S0)*G4Log(S/S0);            //mb
  C        = 59.27;
  d1       = -6.95;
  d2       = 23.54;
  d3       = -25.34;

  xsection = SigAss * ( 1 + 1. / (std::sqrt(S-4.*Mn*Mn)) / (theG4Pow->powN(R0, 3))
                        * C * ( 1 + d1/SqrtS + d2/(theG4Pow->powN(SqrtS, 2))
                                + d3/(theG4Pow->powN(SqrtS, 3)) ) );  

  //xsection *= millibarn;
  fAntiHadronNucleonElXsc = xsection;

  return fAntiHadronNucleonElXsc;
}


/////////////////////////////////////////////////////////////////////////////

void G4ComponentAntiNuclNuclearXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "The G4ComponentAntiNuclNuclearXS calculates total,\n"
          << "inelastic, elastic cross sections  of anti-nucleons and light \n"
          << "anti-nucleus interactions with nuclei using Glauber's approach.\n" 
          << "It uses parametrizations of antiproton-proton total and elastic \n"
          << "cross sections and Wood-Saxon distribution of nuclear density.\n"
          << "See details in Phys.Lett. B705 (2011) 235. \n";
}

