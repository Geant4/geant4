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
// $Id: G4QProbability.cc,v 1.3 2009-09-04 16:13:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QProbability ----------------
//      by Mikhail Kossov Oct, 2006
//   class for Pomeron & Reggeon amplitudes used by CHIPS
//   For comparison mirror member functions are taken from G4 class:
//   G4PomeronCrossSection
// ------------------------------------------------------------------
// Short description: Pomeron is one of the possible vacuum pole (the
// second is Oderon, but they are identical in the present model), by
// which particle exchang in the ellastic scattering process. Others
// are Reggeons and, possibly Instantons (for spin-flip reactions).
// Strings are cuts of Pomerons and Reggeons (optic theorem connects
// the amplitude of scattering at zero angle with the total inelastic
// cross-section). They describe inelastic processes at high energies.
// ------------------------------------------------------------------

#include "G4QProbability.hh"

G4QProbability::G4QProbability(G4int PDG)
{
  S0             = 1.*GeV*GeV;     // Must be a constant (just GeV^2 unit !)
  pom_Alpha      = 1.0808;         // Must be the same for all hadrons
  pom_Alphaprime = 0.25/GeV/GeV;     
  qex_Gamma      = 9./GeV/GeV;
  qex_R2         = 27./GeV/GeV;
  qex_Alphaprime = 1.5/GeV/GeV;     

  G4int aP = std::abs(PDG);
  if     (PDG==2212 || PDG==2112)                     InitForNucleon();
  else if(PDG==111 || aP==211)                        InitForPion();
  else if(PDG==130 || PDG==310 || aP==311 || aP==321) InitForKaon();
  else if(PDG==22)                                    InitForGamma();
  else if(PDG > 3000)                                 InitForHyperon();
  else if(PDG <-2000)                                 InitForAntiBaryon();
  else
  {
    G4cout<<"-Warning-G4QProbability is initialized for PDGCode="<<PDG<<" as Pion"<<G4endl;
    InitForPion();
  }
  pom_sqC=std::sqrt(pom_C);
}

G4double G4QProbability::GetCutPomProbability(const G4double s, const G4double imp2,
                                              const G4int nPom)
{
  static const G4int nft=11;
  static const G4int nf1=nft-1;
  static const G4double ft[nft]={1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3628800.};
  if(nPom<0) return 0.;
  G4double f=ft[nf1];
  if(nPom<nft) f=ft[nPom];
  else for(G4int i=nft; i<= nPom; i++) f*=i;         // Calculate factorial for high nPom
  G4double e=PomEikonal(s,imp2); e+=e;               // Doubled Eikonal
  return std::exp(-e)*std::pow(e,nPom)/pom_C/f;
}

G4double G4QProbability::GetCutQexProbability(const G4double s, const G4double imp2,
                                              const G4int nQex)
{
  static const G4int nft=11;
  static const G4int nf1=nft-1;
  static const G4double ft[nft]={1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3628800.};
  if(nQex<0) return 0.;
  G4double f=ft[nf1];
  if(nQex<nft) f=ft[nQex];
  else for(G4int i=nft; i<= nQex; i++) f*=i;         // Calculate factorial for high nPom
  G4double e=QexEikonal(s,imp2); e+=e;               // Doubled Eikonal
  return std::exp(-e)*std::pow(e,nQex)/f;
}

void G4QProbability::InitForNucleon()
{
  pom_Gamma = 2.16/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.4;
  pom_R2    = 3.30/GeV/GeV;
}

void G4QProbability::InitForHyperon()
{
  pom_Gamma = 2.16/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.4;
  pom_R2    = 3.30/GeV/GeV;       // ? M.K.@@ Just a guess
}
void G4QProbability::InitForAntiBaryon()
{
  pom_Gamma = 2.16/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.4;
  pom_R2    = 3.30/GeV/GeV;
}

void G4QProbability::InitForPion()
{
  pom_Gamma = 2.16/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.6;                // Only: for mesons it is bigger than for baryons
  pom_R2    = 2.36/GeV/GeV;
}

void G4QProbability::InitForKaon()
{
  pom_Gamma = 1.92/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.8;                // 1.7 (?)
  pom_R2    = 1.96/GeV/GeV;       // ? M.K.@@ Just a guess
}

void G4QProbability::InitForGamma()
{
  pom_Gamma = 2.16/GeV/GeV;       // ? M.K.@@ Must be taken from total cross-sections
  pom_C     = 1.7;
  pom_R2    = 2.16/GeV/GeV;       // ? M.K.@@ Just a guess
}

G4double G4QProbability::Expand(G4double z)
{
  G4double sum=1.;
  G4double current=1.;
  for(G4int j=2; j<21; j++)
  {
    current *= -z*(j-1)/j/j;
    sum+=current;
  }
  return sum;
}

G4double G4QProbability::GetQexTotProbability(const G4double s, const G4double imp2)
{
  G4double ExpPom=std::exp(-PomEikonal(s,imp2));
  G4double ExpQex=std::exp(-QexEikonal(s,imp2));
  G4double Amp=(ExpQex*(1.-ExpPom) + sqr(pom_sqC-1.)*ExpPom*(1.-ExpQex))/pom_C;
  return Amp+Amp;
}

G4double G4QProbability::GetQexElProbability(const G4double s, const G4double imp2)
{
  G4double ExpPom=std::exp(-PomEikonal(s,imp2));
  G4double ExpQex=std::exp(-QexEikonal(s,imp2));
  G4double Amp=(ExpQex*(1.-ExpPom) + sqr(pom_sqC-1.)*ExpPom*(1.-ExpQex))/pom_C;
  return Amp*Amp;
}

G4double G4QProbability::GetQexDubDiffProbability(const G4double s, const G4double imp2)
{
  G4double ExpPom=std::exp(-PomEikonal(s,imp2));
  G4double ExpQex=std::exp(-QexEikonal(s,imp2));
  G4double Amp=sqr(pom_sqC-1.)*(ExpQex*(1.-ExpPom) + ExpPom*(1.-ExpQex))/pom_C;
  return Amp*Amp;
}

G4double G4QProbability::GetQexSinDiffProbability(const G4double s, const G4double imp2)
{
  G4double ExpPom=std::exp(-PomEikonal(s,imp2));
  G4double ExpQex=std::exp(-QexEikonal(s,imp2));
  G4double Amp=(pom_sqC-1.)*(ExpQex*(1.-ExpPom) - (pom_sqC-1.)*ExpPom*(1.-ExpQex))/pom_C;
  return Amp*Amp;
}

G4double G4QProbability::GetQexDiffProbability(const G4double s, const G4double imp2)
{
  return GetQexDubDiffProbability(s,imp2)+2*GetQexSinDiffProbability(s,imp2);
}

G4double G4QProbability::GetQexInelProbability(const G4double s, const G4double imp2)
{
  G4double ExpPom=std::exp(-PomEikonal(s,imp2));
  G4double ExpQex=std::exp(-QexEikonal(s,imp2));
  G4double Amp=sqr(pom_sqC-1.)*(ExpQex*(1.-ExpPom) + ExpPom*(1.-ExpQex))/pom_C;
  return Amp+Amp-Amp*Amp;
}
