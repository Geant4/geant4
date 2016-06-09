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
// $Id$
//
//---------------------------------------------------------------------
//
// Geant4 class G4E1Probability
//
// by V. Lara (May 2003)
//
// Modifications:
// 18.05.2010 V.Ivanchenko trying to speedup the most slow method
//            by usage of G4Pow, integer A and introduction of const members
// 17.11.2010 V.Ivanchenko perform general cleanup and simplification
//            of integration method; low-limit of integration is defined
//            by gamma energy or is zero (was always zero before)
//

#include "G4E1Probability.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"

// Constructors and operators
//

G4E1Probability::G4E1Probability():G4VEmissionProbability()
{
  G4double x = CLHEP::pi*CLHEP::hbarc;
  normC = 1.0 / (x*x);
  theLevelDensityParameter = 0.125/MeV;
  fG4pow = G4Pow::GetInstance(); 
}

G4E1Probability::~G4E1Probability()
{}

// Calculate the emission probability
//

G4double G4E1Probability::EmissionProbDensity(const G4Fragment& frag, 
					      G4double gammaE)
{

  // Calculate the probability density here

  // From nuclear fragment properties and the excitation energy, calculate
  // the probability density for photon evaporation from U to U - gammaE
  // (U = nucleus excitation energy, gammaE = total evaporated photon
  // energy). Fragment = nuclear fragment BEFORE de-excitation

  G4double theProb = 0.0;

  G4int Afrag = frag.GetA_asInt();
  G4double Uexcite = frag.GetExcitationEnergy();
  G4double U = std::max(0.0,Uexcite-gammaE);

  if(gammaE < 0.0) { return theProb; }

  // Need a level density parameter.
  // For now, just use the constant approximation (not reliable near magic
  // nuclei).

  G4double aLevelDensityParam = Afrag*theLevelDensityParameter;

  //  G4double levelDensBef = std::exp(2*std::sqrt(aLevelDensityParam*Uexcite));
  //  G4double levelDensAft = std::exp(2*std::sqrt(aLevelDensityParam*(Uexcite-gammaE)));
  // VI reduce number of calls to exp 
  G4double levelDens = 
    std::exp(2*(std::sqrt(aLevelDensityParam*U)-std::sqrt(aLevelDensityParam*Uexcite)));
  // Now form the probability density

  // Define constants for the photoabsorption cross-section (the reverse
  // process of our de-excitation)

  G4double sigma0 = 2.5 * Afrag * millibarn;  // millibarns

  G4double Egdp   = (40.3 / fG4pow->powZ(Afrag,0.2) )*MeV;
  G4double GammaR = 0.30 * Egdp;
 
  // CD
  //cout<<"  PROB TESTS "<<G4endl;
  //cout<<" hbarc = "<<hbarc<<G4endl;
  //cout<<" pi = "<<pi<<G4endl;
  //cout<<" Uexcite, gammaE = "<<Uexcite<<"  "<<gammaE<<G4endl;
  //cout<<" Uexcite, gammaE = "<<Uexcite*MeV<<"  "<<gammaE*MeV<<G4endl;
  //cout<<" lev density param = "<<aLevelDensityParam<<G4endl;
  //cout<<" level densities = "<<levelDensBef<<"  "<<levelDensAft<<G4endl;
  //cout<<" sigma0 = "<<sigma0<<G4endl;
  //cout<<" Egdp, GammaR = "<<Egdp<<"  "<<GammaR<<G4endl;
  //cout<<" normC = "<<normC<<G4endl;

  // VI implementation 18.05.2010
  G4double gammaE2 = gammaE*gammaE;
  G4double gammaR2 = gammaE2*GammaR*GammaR;
  G4double egdp2   = gammaE2 - Egdp*Egdp;
  G4double sigmaAbs = sigma0*gammaR2/(egdp2*egdp2 + gammaR2); 
  theProb = normC * sigmaAbs * gammaE2 * levelDens;

  // old implementation
  //  G4double numerator = sigma0 * gammaE*gammaE * GammaR*GammaR;
  // G4double denominator = (gammaE*gammaE - Egdp*Egdp)*
  //         (gammaE*gammaE - Egdp*Egdp) + GammaR*GammaR*gammaE*gammaE;

  //G4double sigmaAbs = numerator/denominator; 
  //theProb = normC * sigmaAbs * gammaE2 * levelDensAft/levelDensBef;

  // CD
  //cout<<" sigmaAbs = "<<sigmaAbs<<G4endl;
  //cout<<" Probability = "<<theProb<<G4endl;

  return theProb;

}

G4double G4E1Probability::EmissionProbability(const G4Fragment& frag, 
                                              G4double gammaE)
{
  // From nuclear fragment properties and the excitation energy, calculate
  // the probability for photon evaporation down to last ground level.
  // fragment = nuclear fragment BEFORE de-excitation

  G4double upperLim = gammaE;
  G4double lowerLim = 0.0; 

  //G4cout << "G4E1Probability::EmissionProbability:  Emin= " << lowerLim
  //	 << " Emax= " << upperLim << G4endl;
  if( upperLim - lowerLim <= CLHEP::keV ) { return 0.0; } 

  // Need to integrate EmissionProbDensity from lowerLim to upperLim 
  // and multiply by factor 3 (?!)

  G4double integ = 3.0 * EmissionIntegration(frag,lowerLim,upperLim);

  return integ;

}

G4double G4E1Probability::EmissionIntegration(const G4Fragment& frag, 
					      G4double lowLim, G4double upLim)

{
  // Simple integration
  // VI replace by direct integration over 100 point

  const G4int numIters = 100;
  G4double Step = (upLim-lowLim)/G4double(numIters);

  G4double res = 0.0;
  G4double x = lowLim - 0.5*Step;

  for(G4int i = 0; i < numIters; ++i) {
    x += Step;
    res += EmissionProbDensity(frag, x);
  }

  if(res > 0.0) { res /= G4double(numIters); }
  else { res = 0.0; }

  return res;

}


