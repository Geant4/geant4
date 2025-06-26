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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4RiGeAngularGenerator
//
// Authors:       Gerardo Depaola & Ricardo Pacheco
// 
// Creation date: 27 October 2024
//
// -------------------------------------------------------------------
//

#include "G4RiGeAngularGenerator.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include <CLHEP/Units/PhysicalConstants.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RiGeAngularGenerator::G4RiGeAngularGenerator()
  : G4VEmAngularDistribution("RiGeAngularGen")
{}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector& 
G4RiGeAngularGenerator::SampleDirection(const G4DynamicParticle* dp,
                                        G4double gEnergy, G4int, 
                                        const G4Material*)
{
  // Sample gamma angle (Z - axis along the parent particle).
  G4double cost = SampleCosTheta(dp->GetKineticEnergy(), gEnergy, 
                                 dp->GetDefinition()->GetPDGMass());
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  G4double phi  = CLHEP::twopi*G4UniformRand(); 

  fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi), cost);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4RiGeAngularGenerator::SampleCosTheta(G4double primKinEnergy, 
						G4double gEnergy,
						G4double mass)
{
  G4double gam  = 1.0 + primKinEnergy/mass;
  G4double rmax = gam*CLHEP::halfpi*std::min(1.0, gam*mass/gEnergy - 1.0);
  G4double rmax2= rmax*rmax;
  G4double x = G4UniformRand()*rmax2/(1.0 + rmax2);

  return std::cos(std::sqrt(x/(1.0 - x))/gam);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4RiGeAngularGenerator::Sample5DPairDirections(const G4DynamicParticle* dp,
					       G4ThreeVector& dirElectron,
					       G4ThreeVector& dirPositron,
					       const G4double gEnergy, const G4double q2, 
					       const G4double gMomentum,
					       const G4double muFinalMomentum,
					       const G4double muFinalEnergy,
					       const G4double* randNumbs,
					       const G4double* W)
{
  G4double muKinEnergy = dp->GetKineticEnergy();
  G4ThreeVector muMomentumVector = dp->GetMomentum();
  G4double muMomentum = muMomentumVector.mag();
  
  // Electron mass
  G4double eMass = CLHEP::electron_mass_c2;
  G4double eMass2 = eMass*eMass;

  // Muon mass
  G4double muMass = dp->GetDefinition()->GetPDGMass();
  
  G4double muEnergy = muKinEnergy + muMass;
  G4LorentzVector muFinalFourMomentum(muMomentumVector, muEnergy);
  
  G4double mint3 = 0.;
  G4double maxt3 = CLHEP::pi;
  G4double Cmin = std::cos(maxt3);
  G4double Cmax = std::cos(mint3);

  if (randNumbs[7] < W[0]) {
    G4double A1 = -(q2 - 2.*muEnergy*gEnergy);
    G4double B1 = -(2.*gMomentum*muMomentum);
    G4double tginterval = G4Log((A1 + B1)/(A1 - B1))/B1;
    
    G4double costg = (-A1 + (A1 - B1)*G4Exp(B1*tginterval*randNumbs[1]))/B1;
    G4double sintg = std::sqrt((1.0 - costg)*(1.0 + costg));
    G4double phig  = CLHEP::twopi*randNumbs[2];
    G4double sinpg = std::sin(phig);
    G4double cospg = std::cos(phig);  

    G4ThreeVector dirGamma;
    dirGamma.set(sintg*cospg, sintg*sinpg, costg);
    G4LorentzVector gFourMomentum(dirGamma*gMomentum, gEnergy);

    G4double cost5 = -1. + 2.*randNumbs[6];
    G4double phi5 = CLHEP::twopi*randNumbs[8];

    G4LorentzVector eFourMomentumMQ = eDP2(q2, eMass2, eMass2, cost5, phi5);
    G4LorentzVector pFourMomentumMQ = pDP2(eMass2, eFourMomentumMQ);

    G4LorentzVector eFourMomentum = eFourMomentumMQ.boost(gFourMomentum.boostVector());
    G4LorentzVector pFourMomentum = pFourMomentumMQ.boost(gFourMomentum.boostVector());

    dirElectron = eFourMomentum.vect().unit();
    dirPositron = pFourMomentum.vect().unit();

    G4double phi = CLHEP::twopi*randNumbs[3];
    PhiRotation(dirElectron, phi);
    PhiRotation(dirPositron, phi);

  } else if (randNumbs[7] >= W[0] && randNumbs[7] < W[1]) {
    G4double A3 = q2 + 2.*gEnergy*muFinalEnergy;
    G4double B3 = -2.*gMomentum*muFinalMomentum;
    
    G4double tQ3interval = G4Log((A3 + B3)/(A3 - B3))/B3;
    G4double tQMG = (-A3 + (A3 - B3)*G4Exp(B3*tQ3interval*randNumbs[0]))/B3;
    G4double phiQP = CLHEP::twopi*randNumbs[2];
    
    G4double sintQ3 = std::sqrt(1. - tQMG*tQMG);
    G4double cospQP = std::cos(phiQP);
    G4double sinpQP = std::sin(phiQP);
    
    G4double Ap = muMomentum*muMomentum + muFinalMomentum*muFinalMomentum + gMomentum*gMomentum;
    G4double A = Ap + 2.*muFinalMomentum*gMomentum*tQMG;
    G4double B = -2.*muMomentum*gMomentum*sintQ3*cospQP;
    G4double C = -2.*muMomentum*gMomentum*tQMG - 2.*muMomentum*muFinalMomentum; 

    G4double absB = std::abs(B);
    G4double t3interval = (1./(A + C + absB*mint3) - 1./(A + C + absB*maxt3))/absB;
    G4double t3 = (-(A + C) + 1./(1./(A + C + absB*mint3) - absB*t3interval*randNumbs[0]))/absB;
    G4double sint3 = std::sin(t3);
    G4double cost3 = std::cos(t3);

    G4double cost = -sint3*sintQ3*cospQP + cost3*tQMG;
    G4double sint = std::sqrt((1. + cost)*(1. - cost));
    G4double cosp = (sintQ3*cospQP*cost3 + sint3*tQMG)/sint;
    G4double sinp = sintQ3*sinpQP/sint;
    
    G4ThreeVector dirGamma;
    dirGamma.set(sint*cosp, sint*sinp, cost);
    G4LorentzVector gFourMomentum(dirGamma*gMomentum, gEnergy);

    G4double cost5 = -1. + 2.*randNumbs[6];
    G4double phi5 = CLHEP::twopi*randNumbs[8];

    G4LorentzVector eFourMomentumMQ = eDP2(q2, eMass2, eMass2, cost5, phi5);
    G4LorentzVector pFourMomentumMQ = pDP2(eMass2, eFourMomentumMQ);

    G4LorentzVector eFourMomentum = eFourMomentumMQ.boost(gFourMomentum.boostVector());
    G4LorentzVector pFourMomentum = pFourMomentumMQ.boost(gFourMomentum.boostVector());

    dirElectron = eFourMomentum.vect().unit();
    dirPositron = pFourMomentum.vect().unit();

    G4double phi = CLHEP::twopi*randNumbs[3];
    PhiRotation(dirElectron, phi);
    PhiRotation(dirPositron, phi);

  } else if (randNumbs[7] >= W[1] && randNumbs[7] < W[2]) {
    G4double phi3 = CLHEP::twopi*randNumbs[0];
    G4double phi5 = CLHEP::twopi*randNumbs[1];
    G4double phi6 = CLHEP::twopi*randNumbs[2];
    G4double minmuFinalEnergy = muMass;
    G4double muEnergyInterval = muEnergy - 2.*eMass - minmuFinalEnergy;
    G4double muFEnergy = minmuFinalEnergy + muEnergyInterval*randNumbs[3];
    
    G4double mineEnergy = eMass;
    G4double maxeEnergy = muEnergy - muFEnergy - eMass;
    G4double eEnergyInterval = maxeEnergy - mineEnergy;
    G4double eEnergy = mineEnergy + eEnergyInterval*randNumbs[4];
    
    G4double cosp3 = 1.;
    G4double sinp3 = 0.;
    G4double cosp5 = std::cos(phi5);
    G4double sinp5 = std::sin(phi5);
    G4double cosp6 = std::cos(phi6);
    G4double sinp6 = std::sin(phi6);

    G4double muFMomentum = std::sqrt(muFEnergy*muFEnergy - muMass*muMass);
    G4double eMomentum = std::sqrt(eEnergy*eEnergy - eMass*eMass);
    G4double pEnergy = muEnergy - muFEnergy - eEnergy;
    G4double pMomentum = std::sqrt(pEnergy*pEnergy - eMass*eMass);
    
    G4double A3 = -2.*muMass*muMass + 2.*muEnergy*muFEnergy;
    G4double B3 = -2.*std::sqrt(muEnergy*muEnergy - muMass*muMass)*muFMomentum;
    G4double cost3interval = G4Log((A3 + B3*Cmax)/(A3 + B3*Cmin))/B3;

    G4double expanCost3r6 = G4Exp(B3*cost3interval*randNumbs[5]);
    G4double cost3 = A3*(expanCost3r6 - 1.)/B3 + Cmin*expanCost3r6;
    G4double sint3 = std::sqrt((1. - cost3)*(1. + cost3));
    
    G4ThreeVector muFinalMomentumVector(muFMomentum*sint3, 0., muFMomentum*cost3);

    G4LorentzVector muFourMomentum(muMomentumVector, muEnergy);
    muFinalFourMomentum.set(muFinalMomentumVector, muFEnergy);
    G4LorentzVector auxVec1 = muFourMomentum - muFinalFourMomentum;
    G4double A5 = auxVec1.mag2() - 2.*eEnergy*(muEnergy - muFEnergy) +
      2.*muMomentumVector[2]*eMomentum - 2.*muFMomentum*eMomentum*cost3;
    G4double B5 = -2.*muFMomentum*eMomentum*(sint3*cosp3*cosp5 + sint3*sinp3*sinp5);
    G4double absA5 = std::abs(A5);
    G4double absB5 = std::abs(B5);
    G4double mint5 = 0.;
    G4double maxt5 = CLHEP::pi;
    G4double t5interval = G4Log((absA5 + absB5*maxt5)/(absA5 + absB5*mint5))/absB5;
    G4double argexp = absB5*t5interval*randNumbs[6] + G4Log(absA5 + absB5*mint5);
    G4double t5 = -absA5/absB5 + G4Exp(argexp)/absB5;
    G4double sint5 = std::sin(t5);
    G4double cost5 = std::cos(t5);

    dirElectron.set(sint5*cosp5, sint5*sinp5, cost5);
    G4ThreeVector eMomentumVector = eMomentum*dirElectron;

    G4ThreeVector auxVec2 = muMomentumVector - muFinalMomentumVector - eMomentumVector;
    G4double p1mp3mp52 = auxVec2.dot(auxVec2);
    G4double Bp = muFinalMomentum*(sint3*cosp3*cosp6 + sint3*sinp3*sinp6) +
      eMomentum*(sint5*cosp5*cosp6 + sint5*sinp5*sinp6);
    G4double Cp = -muMomentum + muFMomentum*cost3 + eMomentum*cost5;
    G4double A6 = p1mp3mp52 + pMomentum*pMomentum;
    G4double B6 = 2.*pMomentum*Bp;
    G4double C6 = 2.*pMomentum*Cp;
    G4double mint6 = 0.;
    G4double maxt6 = CLHEP::pi;
    G4double absA6C6 = std::abs(A6 + C6);
    G4double absB6 = std::abs(B6);
    G4double t6interval = (1./(absA6C6 + absB6*mint6) - 1./(absA6C6 + absB6*maxt6))/absB6;
    G4double t6 = (-absA6C6 + 1./(1./(absA6C6 + absB6*mint6) - absB6*t6interval*randNumbs[8]))/absB6;
    G4double sint6 = std::sin(t6);
    G4double cost6 = std::cos(t6);
    
    dirPositron.set(sint6*cosp6, sint6*sinp6, cost6);
    
    PhiRotation(dirElectron, phi3);
    PhiRotation(dirPositron, phi3);

  } else if (randNumbs[7] >= W[2]) {
    G4double phi3 = CLHEP::twopi*randNumbs[0];
    G4double phi6 = CLHEP::twopi*randNumbs[1];
    G4double phi5 = CLHEP::twopi*randNumbs[2];
    G4double minmuFinalEnergy = muMass;
    G4double muFinalEnergyinterval = muEnergy - 2.*eMass - minmuFinalEnergy;
    G4double muFEnergy = minmuFinalEnergy + muFinalEnergyinterval*randNumbs[3];
    
    G4double minpEnergy = eMass;
    G4double maxpEnergy = muEnergy - muFEnergy - eMass;
    G4double pEnergyinterval = maxpEnergy - minpEnergy;
    G4double pEnergy = minpEnergy + pEnergyinterval*randNumbs[4];
    
    G4double cosp3 = 1.;
    G4double sinp3 = 0.;
    G4double cosp5 = std::cos(phi5);
    G4double sinp5 = std::sin(phi5);
    G4double cosp6 = std::cos(phi6);
    G4double sinp6 = std::sin(phi6);

    G4double muFMomentum = std::sqrt(muFEnergy*muFEnergy - muMass*muMass);
    G4double pMomentum = std::sqrt(pEnergy*pEnergy - eMass*eMass);
    G4double eEnergy = muEnergy - muFEnergy - pEnergy;
    G4double eMomentum = std::sqrt(eEnergy*eEnergy - eMass*eMass);
    
    G4double A3 = -2.*muMass*muMass + 2.*muEnergy*muFEnergy;
    G4double B3 = -2.*std::sqrt(muEnergy*muEnergy - muMass*muMass)*muFMomentum;
    G4double cost3interval = G4Log((A3 + B3*Cmax)/(A3 + B3*Cmin))/B3;

    G4double expanCost3r6 = G4Exp(B3*cost3interval*randNumbs[5]);
    G4double cost3 = A3*(expanCost3r6 - 1.)/B3 + Cmin*expanCost3r6;
    G4double sint3 = std::sqrt((1. - cost3)*(1. + cost3));

    G4ThreeVector muFinalMomentumVector;
    muFinalMomentumVector.set(muFMomentum*sint3*cosp3, muFMomentum*sint3*sinp3,
			      muFMomentum*cost3);

    G4LorentzVector muFourMomentum(muMomentumVector, muEnergy);
    muFinalFourMomentum.set(muFinalMomentumVector, muFEnergy);
    G4LorentzVector auxVec1 = muFourMomentum - muFinalFourMomentum;
    G4double A6 = auxVec1.mag2() - 2.*pEnergy*(muEnergy - muFEnergy) +
      2.*muMomentumVector[2]*pMomentum - 2.*muFMomentum*pMomentum*cost3;
    G4double B6 = -2.*muFMomentum*pMomentum*(sint3*cosp3*cosp6 + sint3*sinp3*sinp6);
    G4double absA6 = std::abs(A6);
    G4double absB6 = std::abs(B6);
    G4double mint6 = 0.;
    G4double maxt6 = CLHEP::pi;
    G4double t6interval = G4Log((absA6 + absB6*maxt6)/(absA6 + absB6*mint6))/absB6;
    G4double argexp = absB6*t6interval*randNumbs[6] + G4Log(absA6 + absB6*mint6);
    G4double t6 = -absA6/absB6 + G4Exp(argexp)/absB6;
    G4double sint6 = std::sin(t6);
    G4double cost6 = std::cos(t6);

    dirPositron.set(sint6*cosp6, sint6*sinp6, cost6);
    G4ThreeVector pMomentumVector = pMomentum*dirPositron;

    G4ThreeVector auxVec2 = muMomentumVector - muFinalMomentumVector - pMomentumVector;
    G4double p1mp3mp62 = auxVec2.dot(auxVec2);
    G4double Bp = muFMomentum*(sint3*cosp3*cosp5 + sint3*sinp3*sinp5) +
      pMomentum*(sint6*cosp6*cosp5 + sint6*sinp6*sinp5);
    G4double Cp = -muMomentum + muFMomentum*cost3 + pMomentum*cost6;
    G4double A5 = p1mp3mp62 + eMomentum*eMomentum;
    G4double B5 = 2.*eMomentum*Bp;
    G4double C5 = 2.*eMomentum*Cp;
    G4double mint5 = 0.;
    G4double maxt5 = CLHEP::pi;
    G4double absA5C5 = std::abs(A5 + C5);
    G4double absB5 = std::abs(B5);
    G4double t5interval = (1./(absA5C5 + absB5*mint5) - 1./(absA5C5 + absB5*maxt5))/absB5;
    G4double t5 = (-absA5C5 + 1./(1./(absA5C5 + absB5*mint5) - absB5*t5interval*randNumbs[8]))/absB5;
    G4double sint5 = std::sin(t5);
    G4double cost5 = std::cos(t5);

    dirElectron.set(sint5*cosp5, sint5*sinp5, cost5);

    PhiRotation(dirElectron, phi3);
    PhiRotation(dirPositron, phi3);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4RiGeAngularGenerator::PhiRotation(G4ThreeVector& dir, G4double phi)
{
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  G4double newX = dir.x()*cosp + dir.y()*sinp;
  G4double newY = -dir.x()*sinp + dir.y()*cosp;

  dir.setX(newX);
  dir.setY(newY);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LorentzVector G4RiGeAngularGenerator::eDP2(G4double x1, G4double x2,
                                             G4double x3, G4double x4,
                                             G4double x5)
{
  G4double sint = std::sqrt((1.0 - x4)*(1.0 + x4));
  G4double cosp = std::cos(x5);
  G4double sinp = std::sin(x5);

  G4double QJM2 = (x1 + x3 - x2)*(x1 + x3 - x2)/(4.*x1) - x3;

  if (QJM2 < 0.) {
    QJM2 = 1.e-13;
  }

  G4double QJM = std::sqrt(QJM2);
  
  G4LorentzVector x6(QJM*sint*cosp, QJM*sint*sinp, QJM*x4, std::sqrt(x2 + QJM2));

  return x6;
} 
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LorentzVector G4RiGeAngularGenerator::pDP2(G4double x3, const G4LorentzVector& x6)
{
  G4LorentzVector x7(-x6.vect(), std::sqrt(x3 + x6.vect().dot(x6.vect())));
  return x7;
} 
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4RiGeAngularGenerator::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Angular Generator by RiGe algorithm" << G4endl;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
