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

 // G4 Low energy model: n-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch


#include "G4LEnp.hh"
#include "Randomize.hh"
#include "G4ios.hh"

// Initialization of static data arrays:
#include "G4LEnpData.hh"
#include "Randomize.hh"


G4LEnp::G4LEnp() :
  G4HadronicInteraction()
{
  //    theParticleChange.SetNumberOfSecondaries(1);
  
  //    SetMinEnergy(10.*MeV);
  //    SetMaxEnergy(1200.*MeV);
  SetMinEnergy(0.);
  SetMaxEnergy(1200.*GeV);
}

G4LEnp::~G4LEnp()
{
  //    theParticleChange.Clear();
}

G4VParticleChange*
G4LEnp::ApplyYourself(const G4Track& aTrack, G4Nucleus& targetNucleus)
{
    theParticleChange.Initialize(aTrack);

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    G4double P = aParticle->GetTotalMomentum();
    G4double Px = aParticle->GetMomentum().x();
    G4double Py = aParticle->GetMomentum().y();
    G4double Pz = aParticle->GetMomentum().z();
    G4double ek = aParticle->GetKineticEnergy();
    G4ThreeVector theInitial = aParticle->GetMomentum();

    if (verboseLevel > 1) {
      G4double E = aParticle->GetTotalEnergy();
      G4double E0 = aParticle->GetDefinition()->GetPDGMass();
      G4double Q = aParticle->GetDefinition()->GetPDGCharge();
      G4double N = targetNucleus.GetN();
      G4double Z = targetNucleus.GetZ();
      G4cout << "G4LEnp:ApplyYourself: incident particle: "
             << aParticle->GetDefinition()->GetParticleName() << G4endl;
      G4cout << "P = " << P/GeV << " GeV/c"
             << ", Px = " << Px/GeV << " GeV/c"
             << ", Py = " << Py/GeV << " GeV/c"
             << ", Pz = " << Pz/GeV << " GeV/c" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", kinetic energy = " << ek/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
      G4cout << "G4LEnp:ApplyYourself: material:" << G4endl;
      G4cout << "A = " << N
             << ", Z = " << Z
             << ", atomic mass " 
             <<  G4Proton::Proton()->GetPDGMass()/GeV << "GeV" 
             << G4endl;
      //
      // GHEISHA ADD operation to get total energy, mass, charge
      //
      E += G4Proton::Proton()->GetPDGMass();
      G4double E02 = E*E - P*P;
      E0 = sqrt(abs(E02));
      if (E02 < 0)E0 *= -1;
      Q += Z;
      G4cout << "G4LEnp:ApplyYourself: total:" << G4endl;
      G4cout << "E = " << E/GeV << " GeV"
             << ", mass = " << E0/GeV << " GeV"
             << ", charge = " << Q << G4endl;
    }

    // Find energy bin

    G4int je1 = 0;
    G4int je2 = NENERGY - 1;
    ek = ek/GeV;
    do {
      G4int midBin = (je1 + je2)/2;
      if (ek < elab[midBin])
        je2 = midBin;
      else
        je1 = midBin;
    } while (je2 - je1 > 1); 
    //    G4int j;
    //abs(ek-elab[je1]) < abs(ek-elab[je2]) ? j = je1 : j = je2;
    G4double delab = elab[je2] - elab[je1];

    // Sample the angle

    G4float sample = G4UniformRand();
    G4int ke1 = 0;
    G4int ke2 = NANGLE - 1;
    G4double dsig = sig[je2][0] - sig[je1][0];
    G4double rc = dsig/delab;
    G4double b = sig[je1][0] - rc*elab[je1];
    G4double sigint1 = rc*ek + b;
    G4double sigint2 = 0.;

    if (verboseLevel > 1) G4cout << "sample=" << sample << G4endl
                                 << ke1 << " " << ke2 << " " 
                                 << sigint1 << " " << sigint2 << G4endl;

    do {
      G4int midBin = (ke1 + ke2)/2;
      dsig = sig[je2][midBin] - sig[je1][midBin];
      rc = dsig/delab;
      b = sig[je1][midBin] - rc*elab[je1];
      G4double sigint = rc*ek + b;
      if (sample < sigint) {
        ke2 = midBin;
        sigint2 = sigint;
      }
      else {
        ke1 = midBin;
        sigint1 = sigint;
      }
      if (verboseLevel > 1)G4cout << ke1 << " " << ke2 << " " 
                                  << sigint1 << " " << sigint2 << G4endl;
    } while (ke2 - ke1 > 1); 

    // sigint1 and sigint2 should be recoverable from above loop

    //    G4double dsig = sig[je2][ke1] - sig[je1][ke1];
    //    G4double rc = dsig/delab;
    //    G4double b = sig[je1][ke1] - rc*elab[je1];
    //    G4double sigint1 = rc*ek + b;

    //    G4double dsig = sig[je2][ke2] - sig[je1][ke2];
    //    G4double rc = dsig/delab;
    //    G4double b = sig[je1][ke2] - rc*elab[je1];
    //    G4double sigint2 = rc*ek + b;

    dsig = sigint2 - sigint1;
    rc = 1./dsig;
    b = ke1 - rc*sigint1;
    G4double kint = rc*sample + b;
    G4double theta = (0.5 + kint)*pi/180.;

    //    G4int k;
    //abs(sample-sig[j][ke1]) < abs(sample-sig[j][ke2]) ? k = ke1 : k = ke2;
    //    G4double theta = (0.5 + k)*pi/180.;

    if (verboseLevel > 1) {
      G4cout << "   energy bin " << je1 << " energy=" << elab[je1] << G4endl;
      G4cout << "   angle bin " << kint << " angle=" << theta/degree << G4endl;
    }


    // Get the target particle

    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();

    G4double E1 = aParticle->GetTotalEnergy();
    G4double M1 = aParticle->GetDefinition()->GetPDGMass();
    G4double E2 = targetParticle->GetTotalEnergy();
    G4double M2 = targetParticle->GetDefinition()->GetPDGMass();
    G4double totalEnergy = E1 + E2;
    G4double pseudoMass = sqrt(totalEnergy*totalEnergy - P*P);
    // pseudoMass also = sqrt(M1*M1 + M2*M2 + 2*M2*E1)

    // Transform into centre of mass system

    G4double px = (M2/pseudoMass)*Px;
    G4double py = (M2/pseudoMass)*Py;
    G4double pz = (M2/pseudoMass)*Pz;
    G4double p = sqrt(px*px + py*py + pz*pz);

    if (verboseLevel > 1) {
      G4cout << "  E1, M1 (GeV) " << E1/GeV << " " << M1/GeV << G4endl;
      G4cout << "  E2, M2 (GeV) " << E2/GeV << " " << M2/GeV << G4endl;
      G4cout << "  particle  1 momentum in CM " << px/GeV << " " << py/GeV << " "
           << pz/GeV << " " << p/GeV << G4endl;
    }

    // First scatter w.r.t. Z axis
    G4double phi = G4UniformRand()*twopi;
    G4double pxnew = p*sin(theta)*cos(phi);
    G4double pynew = p*sin(theta)*sin(phi);
    G4double pznew = p*cos(theta);

    // Rotate according to the direction of the incident particle
    if (px*px + py*py > 0) {
      G4double cost, sint, ph, cosp, sinp, a, b, c;
      cost = pz/p;
      sint = (sqrt(abs((1-cost)*(1+cost))) + sqrt(px*px+py*py)/p)/2;
      py < 0 ? ph = 3*halfpi : ph = halfpi;
      if (abs(px) > 0.000001*GeV) ph = atan2(py,px);
      cosp = cos(ph);
      sinp = sin(ph);
      px = (cost*cosp*pxnew - sinp*pynew + sint*cosp*pznew);
      py = (cost*sinp*pxnew + cosp*pynew + sint*sinp*pznew);
      pz = (-sint*pxnew                  + cost*pznew);
      //      G4ThreeVector it(a,b,c);
      //      p0->SetMomentum(it);
      //      G4ThreeVector aTargetMom = theInitial - it;
      //      targetParticle->SetMomentum(aTargetMom);
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
    G4double gammaCM = E1pM2/sqrt(E1pM2*E1pM2 - P*P);

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
    PA[4] = sqrt(M1*M1 + p*p);

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
    PA[4] = sqrt(M2*M2 + p*p);

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

    // charge symmetry....
    if(G4UniformRand()<.5)
    {
      theParticleChange.SetNumberOfSecondaries(1);
      theParticleChange.SetMomentumChange(newP->GetMomentumDirection());
      theParticleChange.SetEnergyChange(newP->GetKineticEnergy());
      delete newP;
      G4DynamicParticle* p1 = new G4DynamicParticle;
      p1->SetDefinition(targetParticle->GetDefinition());
      p1->SetMomentum(targetParticle->GetMomentum());
      theParticleChange.AddSecondary(p1);    
    }
    else
    {
      theParticleChange.SetNumberOfSecondaries(2);
      theParticleChange.SetStatusChange(fStopAndKill);
      G4DynamicParticle * pA = new G4DynamicParticle;
      pA->SetDefinition(targetParticle->GetDefinition());
      pA->SetMomentum(newP->GetMomentum());
      G4DynamicParticle * pB = new G4DynamicParticle;
      pB->SetDefinition(newP->GetDefinition());
      pB->SetMomentum(targetParticle->GetMomentum());
      delete newP;
      theParticleChange.AddSecondary(pA);    
      theParticleChange.AddSecondary(pB);    
    }
    return &theParticleChange;
}

 // end of file
