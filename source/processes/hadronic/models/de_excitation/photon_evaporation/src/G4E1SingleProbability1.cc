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
//  Class G4E1SingleProbability1.cc
//

#include "G4E1SingleProbability1.hh"
#include "G4ConstantLevelDensityParameter.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Constructors and operators
//

G4E1SingleProbability1::G4E1SingleProbability1()
{}

G4E1SingleProbability1::~G4E1SingleProbability1()
{}

// Calculate the emission probability
//

G4double G4E1SingleProbability1::EmissionProbDensity(const G4Fragment& frag, 
						     G4double exciteE)
{

  // Calculate the probability density here

  // From nuclear fragment properties and the excitation energy, calculate
  // the probability density for photon evaporation from U to U - exciteE
  // (U = nucleus excitation energy, exciteE = total evaporated photon
  // energy).
  // fragment = nuclear fragment BEFORE de-excitation

  G4double theProb = 0.0;

  G4int Afrag = frag.GetA_asInt();
  G4int Zfrag = frag.GetZ_asInt();
  G4double Uexcite = frag.GetExcitationEnergy();

  if( (Uexcite-exciteE) < 0.0 || exciteE < 0 || Uexcite <= 0) return theProb;

  // Need a level density parameter.
  // For now, just use the constant approximation (not reliable near magic
  // nuclei).

  G4ConstantLevelDensityParameter a;
  G4double aLevelDensityParam = a.LevelDensityParameter(Afrag,Zfrag,Uexcite);

  G4double levelDensBef = std::exp(2.0*std::sqrt(aLevelDensityParam*Uexcite));
  G4double levelDensAft = std::exp(2.0*std::sqrt(aLevelDensityParam*(Uexcite-exciteE)));

  // Now form the probability density

  // Define constants for the photoabsorption cross-section (the reverse
  // process of our de-excitation)

  G4double sigma0 = 2.5 * Afrag * millibarn;  // millibarns

  G4double Egdp = (40.3 / G4Pow::GetInstance()->powZ(Afrag,0.2) )*MeV;
  G4double GammaR = 0.30 * Egdp;
 
  const G4double normC = 1.0 / ((pi * hbarc)*(pi * hbarc));

  // CD
  //cout<<"  PROB TESTS "<<G4endl;
  //cout<<" hbarc = "<<hbarc<<G4endl;
  //cout<<" pi = "<<pi<<G4endl;
  //cout<<" Uexcite, exciteE = "<<Uexcite<<"  "<<exciteE<<G4endl;
  //cout<<" Uexcite, exciteE = "<<Uexcite*MeV<<"  "<<exciteE*MeV<<G4endl;
  //cout<<" lev density param = "<<aLevelDensityParam<<G4endl;
  //cout<<" level densities = "<<levelDensBef<<"  "<<levelDensAft<<G4endl;
  //cout<<" sigma0 = "<<sigma0<<G4endl;
  //cout<<" Egdp, GammaR = "<<Egdp<<"  "<<GammaR<<G4endl;
  //cout<<" normC = "<<normC<<G4endl;

  G4double numerator = sigma0 * exciteE*exciteE * GammaR*GammaR;
  G4double denominator = (exciteE*exciteE - Egdp*Egdp)*
           (exciteE*exciteE - Egdp*Egdp) + GammaR*GammaR*exciteE*exciteE;

  G4double sigmaAbs = numerator/denominator; 

  theProb = normC * sigmaAbs * exciteE*exciteE *
            levelDensAft/levelDensBef;

  // CD
  //cout<<" sigmaAbs = "<<sigmaAbs<<G4endl;
  //cout<<" Probability = "<<theProb<<G4endl;

  return theProb;

}

G4double G4E1SingleProbability1::EmissionProbability(const G4Fragment& frag, 
						     G4double exciteE)
{

  // From nuclear fragment properties and the excitation energy, calculate
  // the probability for photon evaporation down to the level
  // Uexcite-exciteE.
  // fragment = nuclear fragment BEFORE de-excitation

  G4double theProb = 0.0;

  G4double ScaleFactor = 1.0;     // playing with scale factors

  const G4double Uexcite = frag.GetExcitationEnergy();
  G4double Uafter = Uexcite - exciteE;

  G4double normC = 3.0;

  const G4double upperLim = Uexcite;
  const G4double lowerLim = Uafter;
  const G4int numIters = 25;

  // Need to integrate EmissionProbDensity from lowerLim to upperLim 
  // and multiply by normC

  G4double integ = normC *
           EmissionIntegration(frag,exciteE,lowerLim,upperLim,numIters);

  if(integ > 0.0) theProb = integ;

  return theProb * ScaleFactor;

}

G4double G4E1SingleProbability1::EmissionIntegration(const G4Fragment& frag, 
						     G4double ,
						     G4double lowLim, G4double upLim,
						     G4int numIters)

{

  // Simple Gaussian quadrature integration

  G4double x;
  const G4double root3 = 1.0/std::sqrt(3.0);

  G4double Step = (upLim-lowLim)/(2.0*numIters);
  G4double Delta = Step*root3;

  G4double mean = 0.0;

  G4double theInt = 0.0;

  for(G4int i = 0; i < numIters; i++) {

    x = (2*i + 1)/Step;
    G4double E1ProbDensityA = EmissionProbDensity(frag,x+Delta);
    G4double E1ProbDensityB = EmissionProbDensity(frag,x-Delta);

    mean += E1ProbDensityA + E1ProbDensityB;

  }

  if(mean*Step > 0.0) theInt = mean*Step;

  return theInt;

}



