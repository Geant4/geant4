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
//
// Angular distribution according to:
// G. Mao et al., Phys. Rev. C57 (1998) 1938
// G. Mao et al., Phys. Rev. C53 (1996) 2933 
//
// Id: G4AngularDistribution.hh,v 1.16 2000/05/11 19:07:29 pia Exp $ //
//
// -------------------------------------------------------------------

#ifndef G4ANGULARDISTRIBUTION_HH
#define G4ANGULARDISTRIBUTION_HH

#include "globals.hh"
#include "G4VAngularDistribution.hh"


class G4AngularDistribution : public G4VAngularDistribution
{

public:

  // Constructors
  G4AngularDistribution(G4bool symmetrize);

  virtual ~G4AngularDistribution();

  virtual G4double CosTheta(G4double s, G4double m1, G4double m2) const;

protected:
public:        // for testing only...

  G4double DifferentialCrossSection(G4double sIn, G4double m1, G4double m2, G4double cosTheta) const;

  G4double Cross(G4double tpPion, G4double tpSigma, G4double tpOmega,
		 G4double tmPion, G4double tmSigma, G4double tmOmega,
		 G4double bMix_o1, G4double bMix_s1, G4double bMix_Omega,
		 G4double bMix_sm, G4double bMix_oL, G4double bMix_sL,
		 G4double bOmega_0, G4double bOmega_1, G4double bOmega_2,
		 G4double bOmega_3, G4double bOmega_m, G4double bOmega_L) const;

private: 

  G4bool sym;

  // Model parameters

  G4double mPion;
  G4double mSigma;
  G4double mOmega;

  G4double cmPion;
  G4double cmSigma;
  G4double cmOmega;

  G4double gPion;
  G4double gSigma;
  G4double gOmega;

  G4double mNucleon;

  // Variables for pion-Term

  G4double m42;
  G4double mPion2;          
  G4double cmPion2;
  G4double dPion1;
  G4double dPion2;
  G4double cm6gp;
  
  G4double cPion_3;
  G4double cPion_2;
  G4double cPion_1;
  G4double cPion_m;
  G4double cPion_L;
  G4double cPion_0;

  // Variables for sigma-Term 

  G4double mSigma2;
  G4double cmSigma2;
  G4double cmSigma4;
  G4double cmSigma6;
  G4double dSigma1;
  G4double dSigma2;
  G4double dSigma3;
  G4double cm2gs;     
  
  G4double cSigma_3;
  G4double cSigma_2;
  G4double cSigma_1;
  G4double cSigma_m;
  G4double cSigma_L;
  G4double cSigma_0;

  // Variables for omega-Term

  G4double mOmega2;
  G4double cmOmega2;
  G4double cmOmega4;
  G4double cmOmega6;
  G4double dOmega1;
  G4double dOmega2;
  G4double dOmega3;
  G4double sOmega1;
  
  G4double cm2go;
  
  G4double cOmega_3;
  G4double cOmega_2;
  G4double cOmega_1;
  G4double cOmega_m;
  G4double cOmega_L;

  // Variables for mix-Term

  G4double fac1;  
  G4double dMix1;
  G4double dMix2;
  G4double dMix3;
  G4double cMix_o1;
  G4double cMix_s1;
  G4double cMix_Omega;
  G4double cMix_sm; 
  G4double fac2;
  G4double fac3; 
  
  G4double cMix_oLc;
  G4double cMix_oLs;
  G4double cMix_sLc;
  G4double cMix_sLs;

};
#endif














