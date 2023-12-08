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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 25-08-06 New Final State type (refFlag==3 , Legendre (Low Energy) + Probability (High Energy) )
//          is added by T. KOI
// 080904 Add Protection for negative energy results in very low energy ( 1E-6 eV ) scattering by T.
// Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPElasticFS.hh"

#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4HadronicParameters.hh"
#include "G4IonTable.hh"
#include "G4LorentzVector.hh"
#include "G4Nucleus.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4Pow.hh"
#include "G4Proton.hh"
#include "G4ReactionProduct.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Triton.hh"

#include "zlib.h"

G4ParticleHPElasticFS::G4ParticleHPElasticFS()
{
  svtEmax = 0.0;
  dbrcEmax = 0.0;
  dbrcEmin = 0.0;
  dbrcAmin = 0.0;
  dbrcUse = false;
  xsForDBRC = nullptr;

  secID = G4PhysicsModelCatalog::GetModelID("model_NeutronHPElastic");

  hasXsec = false;
  theCoefficients = nullptr;
  theProbArray = nullptr;

  repFlag = 0;
  tE_of_repFlag3 = 0.0;
  targetMass = 0.0;
  frameFlag = 0;
}

void G4ParticleHPElasticFS::Init(G4double A, G4double Z, G4int M,
                                 G4String& dirName, G4String&,
                                 G4ParticleDefinition*)
{
  G4String tString = "/FS";
  G4bool dbool = true;
  SetA_Z(A, Z, M);
  G4ParticleHPDataUsed aFile =
    theNames.GetName(theBaseA, theBaseZ, M, dirName, tString, dbool);
  G4String filename = aFile.GetName();
  SetAZMs(aFile);
  if (!dbool) {
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }

  // 130205 For compressed data files
  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename, theData);
  // 130205 END
  theData >> repFlag >> targetMass >> frameFlag;

  if (repFlag == 1) {
    G4int nEnergy;
    theData >> nEnergy;
    theCoefficients = new G4ParticleHPLegendreStore(nEnergy);
    theCoefficients->InitInterpolation(theData);
    G4double temp, energy;
    G4int tempdep, nLegendre;
    G4int i, ii;
    for (i = 0; i < nEnergy; i++) {
      theData >> temp >> energy >> tempdep >> nLegendre;
      energy *= eV;
      theCoefficients->Init(i, energy, nLegendre);
      theCoefficients->SetTemperature(i, temp);
      G4double coeff = 0;
      for (ii = 0; ii < nLegendre; ii++) {
        // load legendre coefficients.
        theData >> coeff;
        theCoefficients->SetCoeff(i, ii + 1, coeff);  // @@@HPW@@@
      }
    }
  }
  else if (repFlag == 2) {
    G4int nEnergy;
    theData >> nEnergy;
    theProbArray = new G4ParticleHPPartial(nEnergy, nEnergy);
    theProbArray->InitInterpolation(theData);
    G4double temp, energy;
    G4int tempdep, nPoints;
    for (G4int i = 0; i < nEnergy; i++) {
      theData >> temp >> energy >> tempdep >> nPoints;
      energy *= eV;
      theProbArray->InitInterpolation(i, theData);
      theProbArray->SetT(i, temp);
      theProbArray->SetX(i, energy);
      G4double prob, costh;
      for (G4int ii = 0; ii < nPoints; ii++) {
        // fill probability arrays.
        theData >> costh >> prob;
        theProbArray->SetX(i, ii, costh);
        theProbArray->SetY(i, ii, prob);
      }
      theProbArray->DoneSetXY(i);
    }
  }
  else if (repFlag == 3) {
    G4int nEnergy_Legendre;
    theData >> nEnergy_Legendre;
    if (nEnergy_Legendre <= 0) {
      std::stringstream iss;
      iss << "G4ParticleHPElasticFS::Init Data Error repFlag is 3 but nEnergy_Legendre <= 0";
      iss << "Z, A and M of problematic file is " << theNDLDataZ << ", " << theNDLDataA << " and "
          << theNDLDataM << " respectively.";
      throw G4HadronicException(__FILE__, __LINE__, iss.str());
    }
    theCoefficients = new G4ParticleHPLegendreStore(nEnergy_Legendre);
    theCoefficients->InitInterpolation(theData);
    G4double temp, energy;
    G4int tempdep, nLegendre;

    for (G4int i = 0; i < nEnergy_Legendre; i++) {
      theData >> temp >> energy >> tempdep >> nLegendre;
      energy *= eV;
      theCoefficients->Init(i, energy, nLegendre);
      theCoefficients->SetTemperature(i, temp);
      G4double coeff = 0;
      for (G4int ii = 0; ii < nLegendre; ii++) {
        // load legendre coefficients.
        theData >> coeff;
        theCoefficients->SetCoeff(i, ii + 1, coeff);  // @@@HPW@@@
      }
    }

    tE_of_repFlag3 = energy;

    G4int nEnergy_Prob;
    theData >> nEnergy_Prob;
    theProbArray = new G4ParticleHPPartial(nEnergy_Prob, nEnergy_Prob);
    theProbArray->InitInterpolation(theData);
    G4int nPoints;
    for (G4int i = 0; i < nEnergy_Prob; i++) {
      theData >> temp >> energy >> tempdep >> nPoints;
      energy *= eV;

      // consistency check
      if (i == 0)
        // if ( energy != tE_of_repFlag3 ) //110620TK This is too tight for 32bit machines
        if (std::abs(energy - tE_of_repFlag3) / tE_of_repFlag3 > 1.0e-15)
          G4cout << "Warning Transition Energy of repFlag3 is not consistent." << G4endl;

      theProbArray->InitInterpolation(i, theData);
      theProbArray->SetT(i, temp);
      theProbArray->SetX(i, energy);
      G4double prob, costh;
      for (G4int ii = 0; ii < nPoints; ii++) {
        // fill probability arrays.
        theData >> costh >> prob;
        theProbArray->SetX(i, ii, costh);
        theProbArray->SetY(i, ii, prob);
      }
      theProbArray->DoneSetXY(i);
    }
  }
  else if (repFlag == 0) {
    theData >> frameFlag;
  }
  else {
    G4cout << "unusable number for repFlag: repFlag=" << repFlag << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
                              "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
  }
  // 130205 For compressed data files(theData changed from ifstream to istringstream)
  // theData.close();
}

G4HadFinalState* G4ParticleHPElasticFS::ApplyYourself(const G4HadProjectile& theTrack)
{
  if (theResult.Get() == nullptr) theResult.Put(new G4HadFinalState);
  theResult.Get()->Clear();
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile* incidentParticle = &theTrack;
  G4ReactionProduct theNeutron(
    const_cast<G4ParticleDefinition*>(incidentParticle->GetDefinition()));
  theNeutron.SetMomentum(incidentParticle->Get4Momentum().vect());
  theNeutron.SetKineticEnergy(eKinetic);

  G4ThreeVector neuVelo =
    (1. / incidentParticle->GetDefinition()->GetPDGMass()) * theNeutron.GetMomentum();
  G4ReactionProduct theTarget =
    GetBiasedThermalNucleus(targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());

  // Neutron and target defined as G4ReactionProducts
  // Prepare Lorentz transformation to lab

  G4ThreeVector the3Neutron = theNeutron.GetMomentum();
  G4double nEnergy = theNeutron.GetTotalEnergy();
  G4ThreeVector the3Target = theTarget.GetMomentum();
  G4double tEnergy = theTarget.GetTotalEnergy();
  G4ReactionProduct theCMS;
  G4double totE = nEnergy + tEnergy;
  G4ThreeVector the3CMS = the3Target + the3Neutron;
  theCMS.SetMomentum(the3CMS);
  G4double cmsMom = std::sqrt(the3CMS * the3CMS);
  G4double sqrts = std::sqrt((totE - cmsMom) * (totE + cmsMom));
  theCMS.SetMass(sqrts);
  theCMS.SetTotalEnergy(totE);

  // Data come as function of n-energy in nuclear rest frame
  G4ReactionProduct boosted;
  boosted.Lorentz(theNeutron, theTarget);
  eKinetic = boosted.GetKineticEnergy();  // get kinetic energy for scattering
  G4double cosTh = -2;

  if (repFlag == 1) {
    cosTh = theCoefficients->SampleElastic(eKinetic);
  }
  else if (repFlag == 2) {
    cosTh = theProbArray->Sample(eKinetic);
  }
  else if (repFlag == 3) {
    if (eKinetic <= tE_of_repFlag3) {
      cosTh = theCoefficients->SampleElastic(eKinetic);
    }
    else {
      cosTh = theProbArray->Sample(eKinetic);
    }
  }
  else if (repFlag == 0) {
    cosTh = 2. * G4UniformRand() - 1.;
  }
  else {
    G4cout << "Unusable number for repFlag: repFlag=" << repFlag << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
                              "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
  }

  if (cosTh < -1.1) {
    return nullptr;
  }

  G4double phi = twopi * G4UniformRand();
  G4double cosPhi = std::cos(phi);
  G4double sinPhi = std::sin(phi);
  G4double theta = std::acos(cosTh);
  G4double sinth = std::sin(theta);

  if (frameFlag == 1) {
    // Projectile scattering values cosTh are in target rest frame
    // In this frame, do relativistic calculation of scattered projectile and
    // target 4-momenta

    theNeutron.Lorentz(theNeutron, theTarget);
    G4double mN = theNeutron.GetMass();
    G4double Pinit = theNeutron.GetTotalMomentum();  // Incident momentum
    G4double Einit = theNeutron.GetTotalEnergy();  // Incident energy
    G4double mT = theTarget.GetMass();

    G4double ratio = mT / mN;
    G4double sqt = std::sqrt(ratio * ratio - 1.0 + cosTh * cosTh);
    G4double beta = Pinit / (mT + Einit);  // CMS beta
    G4double denom = 1. - beta * beta * cosTh * cosTh;
    G4double term1 = cosTh * (Einit * ratio + mN) / (mN * ratio + Einit);
    G4double pN = beta * mN * (term1 + sqt) / denom;

    // Get the scattered momentum and rotate it in theta and phi
    G4ThreeVector pDir = theNeutron.GetMomentum() / Pinit;
    G4double px = pN * pDir.x();
    G4double py = pN * pDir.y();
    G4double pz = pN * pDir.z();

    G4ThreeVector pcmRot;
    pcmRot.setX(px * cosTh * cosPhi - py * sinPhi + pz * sinth * cosPhi);
    pcmRot.setY(px * cosTh * sinPhi + py * cosPhi + pz * sinth * sinPhi);
    pcmRot.setZ(-px * sinth + pz * cosTh);
    theNeutron.SetMomentum(pcmRot);
    G4double eN = std::sqrt(pN * pN + mN * mN);  // Scattered neutron energy
    theNeutron.SetTotalEnergy(eN);

    // Get the scattered target momentum
    G4ReactionProduct toLab(-1. * theTarget);
    theTarget.SetMomentum(pDir * Pinit - pcmRot);
    G4double eT = Einit - eN + mT;
    theTarget.SetTotalEnergy(eT);

    // Now back to lab frame
    theNeutron.Lorentz(theNeutron, toLab);
    theTarget.Lorentz(theTarget, toLab);

    // 111005 Protection for not producing 0 kinetic energy target
    if (theNeutron.GetKineticEnergy() <= 0)
      theNeutron.SetTotalEnergy(theNeutron.GetMass()
                                * (1. + G4Pow::GetInstance()->powA(10, -15.65)));
    if (theTarget.GetKineticEnergy() <= 0)
      theTarget.SetTotalEnergy(theTarget.GetMass() * (1. + G4Pow::GetInstance()->powA(10, -15.65)));
  }
  else if (frameFlag == 2) {
    // Projectile scattering values cosTh taken from center of mass tabulation

    G4LorentzVector proj(nEnergy, the3Neutron);
    G4LorentzVector targ(tEnergy, the3Target);
    G4ThreeVector boostToCM = proj.findBoostToCM(targ);
    proj.boost(boostToCM);
    targ.boost(boostToCM);

    // Rotate projectile and target momenta by CM scattering angle
    // Note: at this point collision axis is not along z axis, due to
    //       momentum given target nucleus by thermal process
    G4double px = proj.px();
    G4double py = proj.py();
    G4double pz = proj.pz();

    G4ThreeVector pcmRot;
    pcmRot.setX(px * cosTh * cosPhi - py * sinPhi + pz * sinth * cosPhi);
    pcmRot.setY(px * cosTh * sinPhi + py * cosPhi + pz * sinth * sinPhi);
    pcmRot.setZ(-px * sinth + pz * cosTh);
    proj.setVect(pcmRot);
    targ.setVect(-pcmRot);

    // Back to lab frame
    proj.boost(-boostToCM);
    targ.boost(-boostToCM);

    theNeutron.SetMomentum(proj.vect());
    theNeutron.SetTotalEnergy(proj.e());

    theTarget.SetMomentum(targ.vect());
    theTarget.SetTotalEnergy(targ.e());

    // 080904 Add Protection for very low energy (1e-6eV) scattering
    if (theNeutron.GetKineticEnergy() <= 0) {
      theNeutron.SetTotalEnergy(theNeutron.GetMass()
                                * (1. + G4Pow::GetInstance()->powA(10, -15.65)));
    }

    // 080904 Add Protection for very low energy (1e-6eV) scattering
    if (theTarget.GetKineticEnergy() <= 0) {
      theTarget.SetTotalEnergy(theTarget.GetMass() * (1. + G4Pow::GetInstance()->powA(10, -15.65)));
    }
  }
  else {
    G4cout << "Value of frameFlag (1=LAB, 2=CMS): " << frameFlag;
    throw G4HadronicException(__FILE__, __LINE__,
                              "G4ParticleHPElasticFS::ApplyYourSelf frameflag incorrect");
  }

  // Everything is now in the lab frame
  // Set energy change and momentum change
  theResult.Get()->SetEnergyChange(theNeutron.GetKineticEnergy());
  theResult.Get()->SetMomentumChange(theNeutron.GetMomentum().unit());

  // Make recoil a G4DynamicParticle
  auto theRecoil = new G4DynamicParticle;
  theRecoil->SetDefinition(G4IonTable::GetIonTable()->GetIon(static_cast<G4int>(theBaseZ),
                                                             static_cast<G4int>(theBaseA), 0));
  theRecoil->SetMomentum(theTarget.GetMomentum());
  theResult.Get()->AddSecondary(theRecoil, secID);

  // Postpone the tracking of the primary neutron
  theResult.Get()->SetStatusChange(suspend);
  return theResult.Get();
}

void G4ParticleHPElasticFS::InitializeScatteringKernelParameters()
{
  // Initialize DBRC variables
  svtEmax = G4HadronicParameters::Instance()->GetNeutronKineticEnergyThresholdForSVT();
  G4ParticleHPManager* manager = G4ParticleHPManager::GetInstance();
  dbrcUse = manager->GetUseDBRC();
  dbrcEmax = manager->GetMaxEnergyDBRC();
  dbrcEmin = manager->GetMinEnergyDBRC();
  dbrcAmin = manager->GetMinADBRC();
}

G4ReactionProduct G4ParticleHPElasticFS::GetBiasedThermalNucleus(const G4double aMass,
                                                                 G4ThreeVector aVelocity,
                                                                 const G4double temp)
{
  // This new method implements the DBRC (Doppler Broadening Rejection Correction) algorithm
  // on top of the SVT (Sampling of the Velocity of the Target nucleus) algorithm.
  // The SVT algorithm was written by Loic Thulliez (CEA-Saclay) on 2021/05/04 in
  // the method G4Nucleus::GetBiasedThermalNucleus; Marek Zmeskal on 2022/11/30
  // implemented the DBRC algorithm on top of the SVT one.
  // While the SVT algorithm is still present also in G4Nucleus::GetBiasedThermalNucleus,
  // the DBRC algorithm on top of the SVT one has been moved in this new method, in
  // order to avoid a cycle dependency between hadronic/util and hadronic/model/particle_hp.

  InitializeScatteringKernelParameters();

  // Set threshold for SVT algorithm
  G4double E_threshold = svtEmax;
  if (svtEmax == -1.) {
    // If E_neutron <= 400*kB*T (400 is a common value encounter in MC neutron transport code)
    // then apply the Sampling ot the Velocity of the Target (SVT) method;
    // else consider the target nucleus being without motion
    E_threshold = 400.0 * 8.617333262E-11 * temp;
  }

  // If DBRC is enabled and the nucleus is heavy enough, then update the energy threshold
  if (dbrcUse && aMass >= dbrcAmin) {
    E_threshold = std::max(svtEmax, dbrcEmax);
  }

  G4double E_neutron = 0.5 * aVelocity.mag2() * G4Neutron::Neutron()->GetPDGMass();  // E=0.5*m*v2

  G4bool dbrcIsOn = dbrcUse && E_neutron >= dbrcEmin && aMass >= dbrcAmin && E_neutron <= dbrcEmax;

  G4Nucleus aNucleus;
  if (E_neutron > E_threshold || !dbrcIsOn) {
    // Apply only the SVT algorithm, not the DBRC one
    return aNucleus.GetBiasedThermalNucleus(targetMass, aVelocity, temp);
  }

  G4ReactionProduct result;
  result.SetMass(aMass * G4Neutron::Neutron()->GetPDGMass());

  // Beta = sqrt(m/2kT)
  G4double beta =
    std::sqrt(result.GetMass()
              / (2. * 8.617333262E-11 * temp));  // kT E-5[eV] mass E-11[MeV] => beta in [m/s]-1

  // Neutron speed vn
  G4double vN_norm = aVelocity.mag();
  G4double vN_norm2 = vN_norm * vN_norm;
  G4double y = beta * vN_norm;

  // Normalize neutron velocity
  aVelocity = (1. / vN_norm) * aVelocity;

  // Variables for sampling of target speed and SVT rejection
  G4double x2;
  G4double randThresholdSVT;
  G4double vT_norm, vT_norm2, mu;
  G4double acceptThresholdSVT;
  G4double vRelativeSpeed;
  G4double cdf0 = 2. / (2. + std::sqrt(CLHEP::pi) * y);

  // DBRC variables
  G4double xsRelative = -99.;
  G4double randThresholdDBRC = 0.;
  // Calculate max cross-section in interval from  v - 4/beta  to  v + 4/beta  for rejection
  G4double eMin =
    0.5 * G4Neutron::Neutron()->GetPDGMass() * (vN_norm - 4. / beta) * (vN_norm - 4. / beta);
  G4double eMax =
    0.5 * G4Neutron::Neutron()->GetPDGMass() * (vN_norm + 4. / beta) * (vN_norm + 4. / beta);
  G4double xsMax = xsForDBRC->GetMaxY(eMin, eMax);

  do {
    do {
      // Sample the target velocity vT in the laboratory frame
      if (G4UniformRand() < cdf0) {
        // Sample in C45 from https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf
        x2 = -std::log(G4UniformRand() * G4UniformRand());
      }
      else {
        // Sample in C61 from https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf
        G4double ampl = std::cos(CLHEP::pi / 2.0 * G4UniformRand());
        x2 = -std::log(G4UniformRand()) - std::log(G4UniformRand()) * ampl * ampl;
      }

      vT_norm = std::sqrt(x2) / beta;
      vT_norm2 = vT_norm * vT_norm;

      // Sample cosine between the incident neutron and the target in the laboratory frame
      mu = 2. * G4UniformRand() - 1.;

      // Define acceptance threshold for SVT
      vRelativeSpeed = std::sqrt(vN_norm2 + vT_norm2 - 2 * vN_norm * vT_norm * mu);
      acceptThresholdSVT = vRelativeSpeed / (vN_norm + vT_norm);
      randThresholdSVT = G4UniformRand();
    } while (randThresholdSVT >= acceptThresholdSVT);

    // Apply DBRC rejection
    xsRelative = xsForDBRC->GetXsec(0.5 * G4Neutron::Neutron()->GetPDGMass() * vRelativeSpeed
                                    * vRelativeSpeed);
    randThresholdDBRC = G4UniformRand();

  } while (randThresholdDBRC >= xsRelative / xsMax);

  aNucleus.DoKinematicsOfThermalNucleus(mu, vT_norm, aVelocity, result);

  return result;
}
