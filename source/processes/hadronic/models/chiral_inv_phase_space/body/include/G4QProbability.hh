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
#ifndef G4QProbability_h
#define G4QProbability_h 1
//
// $Id: G4QProbability.hh,v 1.4 2009-09-04 16:13:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QProbability ----------------
//      by Mikhail Kossov Oct, 2006
//   class for Pomeron & Reggeon amplitudes used by CHIPS
//   For comparison similar member functions are in the G4 class:
//   G4PomeronCrossSection
// ------------------------------------------------------------
// Short description: Pomeron is one of the possible vacuum pole (the
// second is Oderon, but they are identical in the present model), by
// which particle exchang in the ellastic scattering process. Others
// are Reggeons and, possibly Instantons (for spin-flip reactions).
// Strings are cuts of Pomerons and Reggeons (optic theorem connects
// the amplitude of scattering at zero angle with the total inelastic
// cross-section). They describe inelastic processes at high energies.
// ------------------------------------------------------------------

#include "globals.hh"

class G4QProbability
{
 public:
  G4QProbability(G4int PDGCode = 2212);
  ~G4QProbability(){;}
  void SetS0(G4double aS0)                        {S0 = aS0;}
  void SetPom_Gamma(G4double aPom_Gamma)          {pom_Gamma = aPom_Gamma;}
  void SetGamma(const G4double aGam)              {pom_Gamma=aGam/GeV/GeV;}// @@ Temporary?
  void SetPom_C(G4double aPom_C)                  {pom_C = aPom_C;}
  void SetPom_R2(G4double aPom_R2)                {pom_R2 = aPom_R2;}
  void SetPom_Alpha(G4double aPom_Alpha)          {pom_Alpha = aPom_Alpha;}
  void SetPom_Alphaprime(G4double aPom_Alphaprime){pom_Alphaprime = aPom_Alphaprime;}
  // Genegal (with low energies)
  G4double GetQexTotProbability(const G4double s, const G4double imp2);
  G4double GetQexCohProbability(const G4double s, const G4double imp2);
  G4double GetQexDiffProbability(const G4double s, const G4double imp2);
  G4double GetQexDubDiffProbability(const G4double s, const G4double imp2);
  G4double GetQexSinDiffProbability(const G4double s, const G4double imp2);//For each T & B
  G4double GetQexAbsProbability(const G4double s, const G4double imp2);
  G4double GetQexElProbability(const G4double s, const G4double imp2);
  G4double GetQexInelProbability(const G4double s, const G4double imp2);
  // Only Pomeron (high energies)
  G4double GetPomTotProbability(const G4double s, const G4double imp2)
         {return 2*(1.-std::exp(-PomEikonal(s,imp2)))/pom_C;}
  G4double GetPomCohProbability(const G4double s, const G4double imp2)
         {return sqr(1.-std::exp(-PomEikonal(s,imp2)))/pom_C;}
  G4double GetPomDiffProbability(const G4double s, const G4double imp2)
         {return ((pom_C-1.)/pom_C)*GetPomCohProbability(s,imp2);}
  G4double GetPomDubDiffProbability(const G4double s, const G4double imp2)
         {return (sqr(pom_sqC-1.)/pom_C)*GetPomCohProbability(s,imp2);}
  G4double GetPomSinDiffProbability(const G4double s, const G4double imp2) //For each T & B
         {return ((pom_sqC-1.)/pom_C)*GetPomCohProbability(s,imp2);}
  G4double GetPomAbsProbability(const G4double s, const G4double imp2)
         {return (1.-std::exp(-2*PomEikonal(s,imp2)))/pom_C;}
  G4double GetPomElProbability(const G4double s, const G4double imp2)
         {return GetPomCohProbability(s,imp2)/pom_C;}
  G4double GetPomInelProbability(const G4double s, const G4double imp2)
         {return GetPomDiffProbability(s,imp2) + GetPomAbsProbability(s,imp2);}

  G4double GetCutPomProbability(const G4double s, const G4double ip2, const G4int nPom);
  G4double GetCutQexProbability(const G4double s, const G4double ip2, const G4int nQex);
 private:   
  void InitForNucleon();
  void InitForHyperon();
  void InitForAntiBaryon();
  void InitForPion();
  void InitForKaon();
  void InitForGamma();
 
  G4double Expand(G4double z);
  G4double PowerQex(const G4double s)  {return qex_Gamma/(s/S0);} // qex_Alpha=0 (anti-p?)
  G4double PowerPom(const G4double s)  {return pom_Gamma*std::pow(s/S0, pom_Alpha-1.);}
  G4double SigQex(const G4double s)    {return 8*pi*hbarc_squared*PowerQex(s);}
  G4double SigPom(const G4double s)    {return 8*pi*hbarc_squared*PowerPom(s);}
  G4double LambdaQex(const G4double s) {return qex_R2+qex_Alphaprime*std::log(s/S0);}
  G4double LambdaPom(const G4double s) {return pom_R2+pom_Alphaprime*std::log(s/S0);}
  G4double ZQex(const G4double s)      {return 2*PowerQex(s)/LambdaQex(s);} // qex_C=1.
  G4double ZPom(const G4double s)      {return 2*pom_C*PowerPom(s)/LambdaPom(s);}
  G4double QexEikonal(const G4double s, const G4double imp2)
                         {return ZQex(s)*std::exp(-imp2/LambdaQex(s)/hbarc_squared/4)/2;}
  G4double PomEikonal(G4double s, G4double imp2)
                         {return ZPom(s)*std::exp(-imp2/LambdaPom(s)/hbarc_squared/4)/2;}
  // Body
  G4double S0;
  G4double pom_Gamma;
  G4double pom_C;
  G4double pom_sqC;
  G4double pom_R2;
  G4double pom_Alpha;
  G4double pom_Alphaprime;
  G4double qex_Gamma;
  G4double qex_R2;
  G4double qex_Alphaprime;
};
#endif
