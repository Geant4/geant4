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
// $Id: G4QProbability.hh,v 1.1 2009-07-31 12:43:28 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QProbability ----------------
//      by Mikhail Kossov Oct, 2006
//      class for a Pomeron used by Parton String Models
//   For comparison mirror member functions are taken from G4 class:
//   G4PomeronCrossSection
// ------------------------------------------------------------
// Short description: Pomeron is one of the possible vacuum pole (the
// second is Oderon, but they are identical in the present model), by
// which particle exchang in the ellastic scattering process. Strings,
// which appear as a cut of the Pomeron (optic theorem connects the
// amplitude of scattering at zero angle with the total inelastic
// cross-section), describe most of the processes at high energies.
// ------------------------------------------------------------------

#include "globals.hh"

class G4QProbability
{
 public:
  G4QProbability(G4int PDGCode = 211);
  ~G4QProbability(){;}
  void Pomeron_S(G4double apomeron_S)              {pomeron_S = apomeron_S;}
  void Pomeron_Gamma(G4double apomeron_Gamma)      {pomeron_Gamma = apomeron_Gamma;}
  void Pomeron_C(G4double apomeron_C)              {pomeron_C = apomeron_C;}
  void Pomeron_Rsquare(G4double apomeron_Rsquare)  {pomeron_Rsquare = apomeron_Rsquare;}
  void Pomeron_Alpha(G4double apomeron_Alpha)      {pomeron_Alpha = apomeron_Alpha;}
  void Pomeron_Alphaprime(G4double apom_Alphaprime){pomeron_Alphaprime = apom_Alphaprime;}
  void Pomeron_Gamma_Hard(G4double apom_Gamma_Hard){pomeron_Gamma_Hard = apom_Gamma_Hard;}
  void Pomeron_Alpha_Hard(G4double apom_Alpha_Hard){pomeron_Alpha_Hard = apom_Alpha_Hard;}
  G4double GetTotalCrossSection(const G4double s) {return SigP(s) * Expand(Z(s)/2);}
  G4double GetDiffractiveCrossSection(const G4double s)
                                         {return (pomeron_C-1.)*GetElasticCrossSection(s);}
  G4double GetElasticCrossSection(const G4double s)
                                  {return SigP(s)/pomeron_C*(Expand(Z(s)/2)-Expand(Z(s)));}
  G4double GetInelasticCrossSection(const G4double s)
                                {return GetTotalCrossSection(s)-GetElasticCrossSection(s);}
  G4double GetTotalProbability(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::GetTotalP: Eikonal="<<Eikonal(s,imp2)<<G4endl;
    return 2*(1.-std::exp(-Eikonal(s,imp2)))/pomeron_C;
  }

  G4double GetDiffractiveProbability(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::GetDiffractiveP: pomeron_C="<<pomeron_C<<G4endl;
    return ((pomeron_C-1.)/pomeron_C)*(GetTotalProbability(s,imp2) -
                                                     GetNondiffractiveProbability(s,imp2));
  }

  G4double GetNondiffractiveProbability(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::GetNondifP: s="<<s<<",imp2="<<imp2<<G4endl;
    //G4cout<<"G4QProbability::GetNondifP: Eikonal="<<Eikonal(s,imp2)<<G4endl;
    return (1.-std::exp(-2*Eikonal(s,imp2)))/pomeron_C;
  }
  G4double GetElasticProbability(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::GetElasticP: s="<<s<<",imp2="<<imp2<<G4endl;
    return GetTotalProbability(s,imp2)-GetInelasticProbability(s,imp2);
  }

  G4double GetInelasticProbability(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::GetInelasticP: s="<<s<<",imp2="<<imp2<<G4endl;
    //G4double DP=GetDiffractiveProbability(s,imp2);
    //G4double NDP=GetNondiffractiveProbability(s,imp2);
    //G4cout<<"G4QProbability::GetInelasticP: DP="<<DP<<", NDP="<<NDP<<G4endl;
    //return NDP+DP;
    return GetDiffractiveProbability(s,imp2) + GetNondiffractiveProbability(s,imp2);
  }

  G4double GetCutPomeronProbability(const G4double s,const G4double ip2, const G4int nPom);
  void SetGamma(const G4double agam) {pomeron_Gamma=agam/GeV/GeV;} // @@ Temporary
  G4double SoftEikonal(G4double s, G4double imp2)
                         {return Zsoft(s)/2*std::exp(-imp2/LambdaSoft(s)/hbarc_squared/4);}
  G4double HardEikonal(G4double s, G4double imp2)
                         {return Zhard(s)/2*std::exp(-imp2/LambdaHard(s)/hbarc_squared/4);}
 private: 
  G4double PowerSoft(const G4double s)
                            {return pomeron_Gamma*std::pow(s/pomeron_S, pomeron_Alpha-1.);}
  G4double PowerHard(const G4double s)
                  {return pomeron_Gamma_Hard*std::pow(s/pomeron_S, pomeron_Alpha_Hard-1.);}
  G4double LambdaSoft(const G4double s)
                         {return pomeron_Rsquare+pomeron_Alphaprime*std::log(s/pomeron_S);}
  G4double LambdaHard(const G4double){return pomeron_Rsquare;} // @@ ? M.K.
  G4double Zsoft(const G4double s)   {return 2*pomeron_C*PowerSoft(s)/LambdaSoft(s);}
  G4double Zhard(const G4double s)   {return 2*pomeron_C*PowerHard(s)/LambdaHard(s);}
  
  void InitForNucleon();
  void InitForHyperon();
  void InitForAntiBaryon();
  void InitForPion();
  void InitForKaon();
  void InitForGamma();
 
  G4double Expand(G4double z);
  G4double Z(const G4double s)
  {
    //G4cout<<"G4QProbability::Z: s="<<s<<G4endl;
    //G4cout<<"G4QProbability::Z: Lambda(s)="<<Lambda(s)<<", Power(s)"<<Power(s)<<G4endl;
    return 2*pomeron_C*Power(s)/Lambda(s);
  }
  G4double SigP(const G4double s) {return 8*pi*hbarc_squared*Power(s);}
  G4double Power(const G4double s)
  {
    //G4cout<<"G4QProbability::Power: s="<<s<<", pom_S="<<pomeron_S<<G4endl;
    return pomeron_Gamma*std::pow(s/pomeron_S, pomeron_Alpha-1);
  }
  G4double Lambda(const G4double s)
  {
    //G4cout<<"G4QProbability::Z: s="<<s<<", pomeron_S="<<pomeron_S<<G4endl;
    return pomeron_Rsquare+pomeron_Alphaprime*std::log(s/pomeron_S);
  }
  G4double Eikonal(const G4double s, const G4double imp2)
  {
    //G4cout<<"G4QProbability::Eikonal: s="<<s<<", imp2="<<imp2<<G4endl;
    //G4cout<<"G4QProbability::Eikonal: Lambda(s)="<<Lambda(s)<<",Z(s)"<<Z(s)<<G4endl;
    return std::exp(-imp2/(4*Lambda(s)*hbarc_squared))*Z(s)/2;
  }
  // Body
  G4double pomeron_S;
  G4double pomeron_Gamma;
  G4double pomeron_C;
  G4double pomeron_Rsquare;
  G4double pomeron_Alpha;
  G4double pomeron_Alphaprime;
  G4double pomeron_Gamma_Hard; // @@ ??
  G4double pomeron_Alpha_Hard; // @@ ??
};
#endif
