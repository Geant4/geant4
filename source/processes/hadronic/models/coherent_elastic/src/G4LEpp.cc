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

 // G4 Low energy model: n-n or p-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch

// FWJ 27-AUG-2010: extended Coulomb-suppressed data to 5 GeV

#include "G4LEpp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"

// Initialization of static data arrays:
#include "G4LEppData.hh"

G4LEpp::G4LEpp():G4HadronElastic("G4LEpp")
{
  SetMinEnergy(0.);
  SetMaxEnergy(5.*GeV);
}

G4LEpp::~G4LEpp()
{}

G4HadFinalState*
G4LEpp::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
    theParticleChange.Clear();
    const G4HadProjectile* aParticle = &aTrack;

    G4double P = aParticle->GetTotalMomentum();
    G4double Px = aParticle->Get4Momentum().x();
    G4double Py = aParticle->Get4Momentum().y();
    G4double Pz = aParticle->Get4Momentum().z();
    G4double E  = aParticle->GetTotalEnergy();
    G4ThreeVector theInitial = aParticle->Get4Momentum().vect().unit();

    if (verboseLevel > 1) {
      G4double ek = aParticle->GetKineticEnergy();
      G4double E0 = aParticle->GetDefinition()->GetPDGMass();
      G4double Q = aParticle->GetDefinition()->GetPDGCharge();
      G4int A = targetNucleus.GetA_asInt();
      G4int Z = targetNucleus.GetZ_asInt();
      G4cout << "G4LEpp:ApplyYourself: incident particle: "
             << aParticle->GetDefinition()->GetParticleName() << G4endl;
      G4cout << "P = " << P/GeV << " GeV/c"
             << ", Px = " << Px/GeV << " GeV/c"
             << ", Py = " << Py/GeV << " GeV/c"
             << ", Pz = " << Pz/GeV << " GeV/c" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", kinetic energy = " << ek/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
      G4cout << "G4LEpp:ApplyYourself: material:" << G4endl;
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
      E0 = std::sqrt(std::fabs(E02));
      if (E02 < 0)E0 *= -1;
      Q += Z;
      G4cout << "G4LEpp:ApplyYourself: total:" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
    }
    G4double t = SampleInvariantT(aParticle->GetDefinition(), P, 0, 0);
    G4double cost = 1.0 - 2*t/(P*P);
    if(cost > 1.0) { cost = 1.0; }
    if(cost <-1.0) { cost =-1.0; }
    G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
    G4double phi = twopi*G4UniformRand();
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
      G4cout << "  particle  1 momentum in CM " << px/GeV 
	     << " " << py/GeV << " "
	     << pz/GeV << " " << p/GeV << G4endl;
    }

    // First scatter w.r.t. Z axis
    G4double pxnew = p*sint*std::cos(phi);
    G4double pynew = p*sint*std::sin(phi);
    G4double pznew = p*cost;

    // Rotate according to the direction of the incident particle
    if (px*px + py*py > 0) {
      G4double ph, cosp, sinp;
      cost = pz/p;
      sint = (std::sqrt((1-cost)*(1+cost)) + std::sqrt(px*px+py*py)/p)/2;
      py < 0 ? ph = 3*halfpi : ph = halfpi;
      if (std::fabs(px) > 0.000001*GeV) ph = std::atan2(py,px);
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
      G4cout << "  particle 1 momentum in CM " << px/GeV << " " << py/GeV << " "
           << pz/GeV << " " << p/GeV << G4endl;
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
           << newP->GetMomentum()/GeV 
           << " " << newP->GetTotalMomentum()/GeV << G4endl;
      G4cout << "  particle 2 momentum in LAB " 
           << targetParticle->GetMomentum()/GeV 
           << " " << targetParticle->GetTotalMomentum()/GeV << G4endl;
      G4cout << "  TOTAL momentum in LAB " 
           << (newP->GetMomentum()+targetParticle->GetMomentum())/GeV 
           << " " 
           << (newP->GetMomentum()+targetParticle->GetMomentum()).mag()/GeV
           << G4endl;
    }

    theParticleChange.SetMomentumChange( newP->GetMomentumDirection());
    theParticleChange.SetEnergyChange(newP->GetKineticEnergy());
    delete newP;

    // Recoil particle
    theParticleChange.AddSecondary(targetParticle);    
    return &theParticleChange;
}

////////////////////////////////////////////////////////////////////
//
// sample momentum transfer using Lab. momentum

G4double G4LEpp::SampleInvariantT(const G4ParticleDefinition* p, 
				  G4double plab, G4int , G4int )
{
  G4double nMass = p->GetPDGMass(); // 939.565346*MeV;
  G4double ek = std::sqrt(plab*plab+nMass*nMass) - nMass;

    // Find energy bin

  G4int je1 = 0;
  G4int je2 = NENERGY - 1;
  ek /= GeV;

  do 
  {
    G4int midBin = (je1 + je2)/2;

    if (ek < elab[midBin]) je2 = midBin;
    else                   je1 = midBin;
  } 
  while (je2 - je1 > 1);  /* Loop checking, 10.08.2015, A.Ribon */

  G4double delab = elab[je2] - elab[je1];

    // Sample the angle

  G4double sample = G4UniformRand();
  G4int ke1 = 0;
  G4int ke2 = NANGLE - 1;
  G4double dsig, b, rc;

  dsig = Sig[je2][0] - Sig[je1][0]; 
  rc = dsig/delab;
  b = Sig[je1][0] - rc*elab[je1];

  G4double sigint1 = rc*ek + b;
  G4double sigint2 = 0.;

  do
  {
      G4int midBin = (ke1 + ke2)/2;
      dsig = Sig[je2][midBin] - Sig[je1][midBin]; 
      rc = dsig/delab;
      b = Sig[je1][midBin] - rc*elab[je1];
      G4double sigint = rc*ek + b;

      if (sample < sigint) 
      {
        ke2 = midBin;
        sigint2 = sigint;
      }
      else 
      {
        ke1 = midBin;
        sigint1 = sigint;
      }
  } 
  while (ke2 - ke1 > 1);  /* Loop checking, 10.08.2015, A.Ribon */ 

  dsig = sigint2 - sigint1;
  rc = 1./dsig;
  b = ke1 - rc*sigint1;

  G4double kint = rc*sample + b;
  G4double theta = (0.5 + kint)*pi/180.;
  G4double t = 0.5*plab*plab*(1 - std::cos(theta));

  return t;
}
// end of file
