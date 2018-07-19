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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      File name:     G4PolarizationTransition
//
//      Author:        Jason Detwiler (jasondet@gmail.com)
// 
//      Creation date: Aug 2012
//
//      Description:   
//      Stores and manipulates the statistical tensor describing the nuclear
//      polarization (see Alder and Winther, "Electromagnetic Excitation" (1975),
//      Appendix F). Functions are implemented for generating angular correlations
//      in gamma decays, following Alder and Winther, Appendix G.
//      This code assumes no polarization will be detected and uses eqs (17-20).
//      Adding polarization would require using instead (13) and the more generic
//      form of the statstical tensor after decay described by equation (6)
//      Could be expanded to also generate e.g. gamma-beta and other
//      correlations as well.
//
// -------------------------------------------------------------------

#ifndef G4POLARIZATIONTRANSITION_HH
#define G4POLARIZATIONTRANSITION_HH

#include "globals.hh"
#include "G4LegendrePolynomial.hh"
#include "G4PolynomialPDF.hh"
#include "G4Pow.hh"

class G4NuclearPolarization;

class G4PolarizationTransition
{
  typedef std::vector< std::vector<G4complex> > POLAR;

  public:
    explicit G4PolarizationTransition();
    ~G4PolarizationTransition();

    void SampleGammaTransition(G4NuclearPolarization* np, 
			       G4int twoJ1, G4int twoJ2, 
                               G4int L0, G4int Lp, G4double mpRatio, 
			       G4double& cosTheta, G4double& phi);

    // generic static functions
    G4double FCoefficient(G4int K, G4int L, G4int Lprime, 
			  G4int twoJ2, G4int twoJ1) const;
    G4double F3Coefficient(G4int K, G4int K2, G4int K1, G4int L, 
			   G4int Lprime, G4int twoJ2, G4int twoJ1) const;

    // transition-specific functions
    G4double GammaTransFCoefficient(G4int K) const;
    G4double GammaTransF3Coefficient(G4int K, G4int K2, G4int K1) const;

    void DumpTransitionData(const POLAR& pol) const;

    inline void SetVerbose(G4int val) { fVerbose = val; };

  private:

    G4PolarizationTransition(const G4PolarizationTransition &right) = delete;
    const G4PolarizationTransition& operator=(const G4PolarizationTransition &right) = delete;

    // Gamma angle generation and decay: call these functions in this order!
    // All angles are in the same coordinate system: user may choose any axis
    G4double GenerateGammaCosTheta(const POLAR&);
    G4double GenerateGammaPhi(G4double& cosTheta, const POLAR&);

    inline G4double LnFactorial(int k) const { return G4Pow::GetInstance()->logfactorial(k); }

    G4int fVerbose;
    G4int fTwoJ1, fTwoJ2;
    G4int fLbar, fL;
    G4double fDelta;
    G4double kEps;
    G4PolynomialPDF kPolyPDF;
    G4LegendrePolynomial fgLegendrePolys;
};


#endif
