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

// G4 Low energy model: n-p scattering
// F.W. Jones, L.G. Greeniaus, H.P. Wellisch

// 11-OCT-2007 F.W. Jones: removed erroneous code for identity
//             exchange of particles.
// FWJ 27-AUG-2010: extended to 5 GeV by Tony Kwan TRIUMF
// 06-Aug-15 A.Ribon migrating to G4Pow

#include "G4LEHadronProtonElastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"

#include "G4Pow.hh"


G4LEHadronProtonElastic::G4LEHadronProtonElastic():
 G4HadronElastic("G4LEHadronProtonElastic") 
{
  SetMinEnergy(0.);
  SetMaxEnergy(20.*MeV);
}

G4LEHadronProtonElastic::~G4LEHadronProtonElastic()
{
      theParticleChange.Clear();
}

G4HadFinalState*
G4LEHadronProtonElastic::ApplyYourself(const G4HadProjectile& aTrack, 
                            G4Nucleus& targetNucleus)
{
    theParticleChange.Clear();
    const G4HadProjectile* aParticle = &aTrack;

    G4double P = aParticle->GetTotalMomentum();
    G4double Px = aParticle->Get4Momentum().x();
    G4double Py = aParticle->Get4Momentum().y();
    G4double Pz = aParticle->Get4Momentum().z();
    G4double ek = aParticle->GetKineticEnergy();
    G4ThreeVector theInitial = aParticle->Get4Momentum().vect();

    if (verboseLevel > 1) 
    {
      G4double E = aParticle->GetTotalEnergy();
      G4double E0 = aParticle->GetDefinition()->GetPDGMass();
      G4double Q = aParticle->GetDefinition()->GetPDGCharge();
      G4int A = targetNucleus.GetA_asInt();
      G4int Z = targetNucleus.GetZ_asInt();
      G4cout << "G4LEHadronProtonElastic:ApplyYourself: incident particle: "
             << aParticle->GetDefinition()->GetParticleName() << G4endl;
      G4cout << "P = " << P/GeV << " GeV/c"
             << ", Px = " << Px/GeV << " GeV/c"
             << ", Py = " << Py/GeV << " GeV/c"
             << ", Pz = " << Pz/GeV << " GeV/c" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", kinetic energy = " << ek/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
      G4cout << "G4LEHadronProtonElastic:ApplyYourself: material:" << G4endl;
      G4cout << "A = " << A
             << ", Z = " << Z
             << ", atomic mass " 
             <<  G4Proton::Proton()->GetPDGMass()/GeV << "GeV" 
             << G4endl;
      //
      // GHEISHA ADD operation to get total energy, mass, charge
      //
      E += proton_mass_c2;
      G4double E02 = E*E - P*P;
      E0 = std::sqrt(std::abs(E02));
      if (E02 < 0)E0 *= -1;
      Q += Z;
      G4cout << "G4LEHadronProtonElastic:ApplyYourself: total:" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
    }

    G4double theta = (0.5)*pi/180.;

    // Get the target particle

    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();

    G4double E1 = aParticle->GetTotalEnergy();
    G4double M1 = aParticle->GetDefinition()->GetPDGMass();
    G4double E2 = targetParticle->GetTotalEnergy();
    G4double M2 = targetParticle->GetDefinition()->GetPDGMass();
    G4double totalEnergy = E1 + E2;
    G4double pseudoMass = std::sqrt(totalEnergy*totalEnergy - P*P);

    // Transform into centre of mass system

    G4double px = (M2/pseudoMass)*Px;
    G4double py = (M2/pseudoMass)*Py;
    G4double pz = (M2/pseudoMass)*Pz;
    G4double p = std::sqrt(px*px + py*py + pz*pz);

    if (verboseLevel > 1) {
      G4cout << "  E1, M1 (GeV) " << E1/GeV << " " << M1/GeV << G4endl;
      G4cout << "  E2, M2 (GeV) " << E2/GeV << " " << M2/GeV << G4endl;
      G4cout << "  particle  1 momentum in CM " << px/GeV << " " << py/GeV << " "
           << pz/GeV << " " << p/GeV << G4endl;
    }

    // First scatter w.r.t. Z axis
    G4double phi = G4UniformRand()*twopi;
    G4double pxnew = p*std::sin(theta)*std::cos(phi);
    G4double pynew = p*std::sin(theta)*std::sin(phi);
    G4double pznew = p*std::cos(theta);

    // Rotate according to the direction of the incident particle
    if (px*px + py*py > 0) 
    {
      G4double cost, sint, ph, cosp, sinp;
      cost = pz/p;
      sint = (std::sqrt(std::fabs((1-cost)*(1+cost))) 
            + std::sqrt(px*px+py*py)/p)/2;
      py < 0 ? ph = 3*halfpi : ph = halfpi;
      if (std::abs(px) > 0.000001*GeV) ph = std::atan2(py,px);
      cosp = std::cos(ph);
      sinp = std::sin(ph);
      px = (cost*cosp*pxnew - sinp*pynew + sint*cosp*pznew);
      py = (cost*sinp*pxnew + cosp*pynew + sint*sinp*pznew);
      pz = (-sint*pxnew                  + cost*pznew);
    }
    else {
      px = pxnew;
      py = pynew;
      pz = pznew;
    }

    if (verboseLevel > 1) {
      G4cout << "  AFTER SCATTER..." << G4endl;
      G4cout << "  particle 1 momentum in CM " << px/GeV 
             << " " << py/GeV << " " << pz/GeV << " " << p/GeV 
             << G4endl;
    }

    // Transform to lab system

    G4double E1pM2 = E1 + M2;
    G4double betaCM  = P/E1pM2;
    G4double betaCMx = Px/E1pM2;
    G4double betaCMy = Py/E1pM2;
    G4double betaCMz = Pz/E1pM2;
    G4double gammaCM = E1pM2/std::sqrt(E1pM2*E1pM2 - P*P);

    if (verboseLevel > 1) {
      G4cout << "  betaCM " << betaCMx << " " << betaCMy << " "
             << betaCMz << " " << betaCM << G4endl;
      G4cout << "  gammaCM " << gammaCM << G4endl;
    }

    // Now following GLOREN...

    G4double BETA[5], PA[5], PB[5];
    BETA[1] = -betaCMx;
    BETA[2] = -betaCMy;
    BETA[3] = -betaCMz;
    BETA[4] = gammaCM;

    //The incident particle...

    PA[1] = px;
    PA[2] = py;
    PA[3] = pz;
    PA[4] = std::sqrt(M1*M1 + p*p);

    G4double BETPA  = BETA[1]*PA[1] + BETA[2]*PA[2] + BETA[3]*PA[3];
    G4double BPGAM  = (BETPA * BETA[4]/(BETA[4] + 1.) - PA[4]) * BETA[4];

    PB[1] = PA[1] + BPGAM  * BETA[1];
    PB[2] = PA[2] + BPGAM  * BETA[2];
    PB[3] = PA[3] + BPGAM  * BETA[3];
    PB[4] = (PA[4] - BETPA) * BETA[4];

    G4DynamicParticle* newP = new G4DynamicParticle;
    newP->SetDefinition(aParticle->GetDefinition());
    newP->SetMomentum(G4ThreeVector(PB[1], PB[2], PB[3]));

    //The target particle...

    PA[1] = -px;
    PA[2] = -py;
    PA[3] = -pz;
    PA[4] = std::sqrt(M2*M2 + p*p);

    BETPA  = BETA[1]*PA[1] + BETA[2]*PA[2] + BETA[3]*PA[3];
    BPGAM  = (BETPA * BETA[4]/(BETA[4] + 1.) - PA[4]) * BETA[4];

    PB[1] = PA[1] + BPGAM  * BETA[1];
    PB[2] = PA[2] + BPGAM  * BETA[2];
    PB[3] = PA[3] + BPGAM  * BETA[3];
    PB[4] = (PA[4] - BETPA) * BETA[4];

    targetParticle->SetMomentum(G4ThreeVector(PB[1], PB[2], PB[3]));

    if (verboseLevel > 1) {
      G4cout << "  particle 1 momentum in LAB " 
           << newP->GetMomentum()*(1./GeV) 
           << " " << newP->GetTotalMomentum()/GeV << G4endl;
      G4cout << "  particle 2 momentum in LAB " 
           << targetParticle->GetMomentum()*(1./GeV) 
           << " " << targetParticle->GetTotalMomentum()/GeV << G4endl;
      G4cout << "  TOTAL momentum in LAB " 
           << (newP->GetMomentum()+targetParticle->GetMomentum())*(1./GeV) 
           << " " 
           << (newP->GetMomentum()+targetParticle->GetMomentum()).mag()/GeV
           << G4endl;
    }

    theParticleChange.SetMomentumChange(newP->GetMomentumDirection());
    theParticleChange.SetEnergyChange(newP->GetKineticEnergy());
    delete newP;
    theParticleChange.AddSecondary(targetParticle);    

    return &theParticleChange;
}

////////////////////////////////////////////////////////////////////
//
// sample momentum transfer using Lab. momentum

G4double 
G4LEHadronProtonElastic::SampleInvariantT(const G4ParticleDefinition* p, 
					  G4double plab, G4int , G4int )
{
  G4double hMass = p->GetPDGMass();
  G4double pCMS = 0.5*plab;
  // pCMS *= 50;
  G4double hEcms = std::sqrt(pCMS*pCMS+hMass*hMass);
  // G4double gamma = hEcms/hMass;
  // gamma  *= 15;
  G4double beta =  pCMS/hEcms; // std::sqrt(1-1./gamma/gamma); //
  // beta /= 0.8; // 0.95; // 1.0; // 1.1 // 0.5*pi; // pi; twopi;
  G4double cosDipole =  RandCosThetaDipPen();
  G4double cosTheta = cosDipole + beta;
  cosTheta /= 1. + cosDipole*beta;
  G4double t = 2.*pCMS*pCMS*(1.-cosTheta);
  return t;

}

///////////////////////////////////////////////////////////////
//
// 1 + cos^2(theta) random distribution in the projectile rest frame, Penelope algorithm

G4double G4LEHadronProtonElastic::RandCosThetaDipPen()
{
  G4double x, cosTheta, signX, modX, power = 1./3.;

  if( G4UniformRand() > 0.25) 
  {
    cosTheta = 2.*G4UniformRand()-1.;
  }
  else
  {
    x     = 2.*G4UniformRand()-1.;

    if ( x < 0. ) 
    {
      modX = -x;
      signX = -1.;
    }
    else
    {
      modX = x;
      signX = 1.;
    }
    cosTheta = signX*G4Pow::GetInstance()->powA(modX,power);
  }
  return cosTheta; 
}

 // end of file
