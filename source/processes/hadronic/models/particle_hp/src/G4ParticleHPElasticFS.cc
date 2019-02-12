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
// 080904 Add Protection for negative energy results in very low energy ( 1E-6 eV ) scattering by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPElasticFS.hh"
#include "G4ParticleHPManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4IonTable.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4Pow.hh"
#include "zlib.h"

void G4ParticleHPElasticFS::Init(G4double A, G4double Z, G4int M,
                                 G4String& dirName, G4String&,
                                 G4ParticleDefinition* )
{
  G4String tString = "/FS";
  G4bool dbool;
  G4ParticleHPDataUsed aFile =
    theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, dirName, tString, dbool);
  G4String filename = aFile.GetName();
  SetAZMs( A, Z, M, aFile ); 
    //theBaseA = aFile.GetA();
    //theBaseZ = aFile.GetZ();
  if (!dbool) {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }

  //130205 For compressed data files 
  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);
  //130205 END
  theData >> repFlag >> targetMass >> frameFlag;

  if (repFlag == 1) {
    G4int nEnergy;
    theData >> nEnergy; 
    theCoefficients = new G4ParticleHPLegendreStore(nEnergy);
    theCoefficients->InitInterpolation(theData);
    G4double temp, energy;
    G4int tempdep, nLegendre;
    G4int i, ii;
    for (i=0; i < nEnergy; i++) {
      theData >> temp >> energy >> tempdep >> nLegendre;
      energy *=eV;
      theCoefficients->Init(i, energy, nLegendre);
      theCoefficients->SetTemperature(i, temp);
      G4double coeff = 0;
      for (ii = 0; ii < nLegendre; ii++) {
        // load legendre coefficients.
        theData >> coeff;
        theCoefficients->SetCoeff(i, ii+1, coeff); // @@@HPW@@@
      }
    }

  } else if (repFlag == 2) {
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
      theProbArray->DoneSetXY( i );
    }

  } else if (repFlag == 3) {
    G4int nEnergy_Legendre;
    theData >> nEnergy_Legendre; 
    if (nEnergy_Legendre <= 0 ) {
      std::stringstream iss;
      iss << "G4ParticleHPElasticFS::Init Data Error repFlag is 3 but nEnergy_Legendre <= 0";
      iss << "Z, A and M of problematic file is " << theNDLDataZ << ", " 
          << theNDLDataA << " and " << theNDLDataM << " respectively.";
      throw G4HadronicException(__FILE__, __LINE__, iss.str() );
    }
    theCoefficients = new G4ParticleHPLegendreStore( nEnergy_Legendre );
    theCoefficients->InitInterpolation( theData );
    G4double temp, energy;
    G4int tempdep, nLegendre;

    for (G4int i = 0; i < nEnergy_Legendre; i++) {
      theData >> temp >> energy >> tempdep >> nLegendre;
      energy *=eV;
      theCoefficients->Init( i , energy , nLegendre );
      theCoefficients->SetTemperature( i , temp );
      G4double coeff = 0;
      for (G4int ii = 0; ii < nLegendre; ii++) {
        // load legendre coefficients.
        theData >> coeff;
        theCoefficients->SetCoeff(i, ii+1, coeff); // @@@HPW@@@
      }
    } 

    tE_of_repFlag3 = energy; 

    G4int nEnergy_Prob;
    theData >> nEnergy_Prob;
    theProbArray = new G4ParticleHPPartial( nEnergy_Prob , nEnergy_Prob );
    theProbArray->InitInterpolation( theData );
    G4int nPoints;
    for (G4int i = 0; i < nEnergy_Prob; i++) {
      theData >> temp >> energy >> tempdep >> nPoints;
      energy *= eV;

      // consistency check
      if (i == 0)
        //if ( energy != tE_of_repFlag3 ) //110620TK This is too tight for 32bit machines 
        if (std::abs(energy - tE_of_repFlag3) / tE_of_repFlag3 > 1.0e-15)
          G4cout << "Warning Transition Energy of repFlag3 is not consistent." << G4endl; 

      theProbArray->InitInterpolation( i , theData );
      theProbArray->SetT( i , temp );
      theProbArray->SetX( i , energy );
      G4double prob, costh;
      for (G4int ii = 0; ii < nPoints; ii++) {
        // fill probability arrays.
        theData >> costh >> prob;
        theProbArray->SetX( i , ii , costh );
        theProbArray->SetY( i , ii , prob );
      }
      theProbArray->DoneSetXY( i );
    }

  } else if (repFlag==0) {
      theData >> frameFlag;

  } else {
      G4cout << "unusable number for repFlag: repFlag="<<repFlag<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
  }
   //130205 For compressed data files(theData changed from ifstream to istringstream)
   //theData.close();
}


G4HadFinalState*
G4ParticleHPElasticFS::ApplyYourself(const G4HadProjectile& theTrack)
{  
  if (theResult.Get() == NULL) theResult.Put(new G4HadFinalState);
  theResult.Get()->Clear();
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile *incidentParticle = &theTrack;
  G4ReactionProduct theNeutron(const_cast<G4ParticleDefinition*>(incidentParticle->GetDefinition() ));
  theNeutron.SetMomentum(incidentParticle->Get4Momentum().vect() );
  theNeutron.SetKineticEnergy(eKinetic);

  G4ReactionProduct theTarget; 
  G4Nucleus aNucleus;
  G4ThreeVector neuVelo =
    (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
  theTarget =
    aNucleus.GetBiasedThermalNucleus(targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());

  // Neutron and target defined as G4ReactionProducts
  // Prepare Lorentz transformation to lab

  G4ThreeVector the3Neutron = theNeutron.GetMomentum();
  G4double nEnergy = theNeutron.GetTotalEnergy();
  G4ThreeVector the3Target = theTarget.GetMomentum();
  G4double tEnergy = theTarget.GetTotalEnergy();
  G4ReactionProduct theCMS;
  G4double totE = nEnergy+tEnergy;
  G4ThreeVector the3CMS = the3Target+the3Neutron;
  theCMS.SetMomentum(the3CMS);
  G4double cmsMom = std::sqrt(the3CMS*the3CMS);
  G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
  theCMS.SetMass(sqrts);
  theCMS.SetTotalEnergy(totE);
    
  // Data come as function of n-energy in nuclear rest frame
  G4ReactionProduct boosted;
  boosted.Lorentz(theNeutron, theTarget);
  eKinetic = boosted.GetKineticEnergy(); // get kinetic energy for scattering
  G4double cosTh = -2;

  if (repFlag == 1) {
    cosTh = theCoefficients->SampleElastic(eKinetic);

  } else if (repFlag == 2) {
    cosTh = theProbArray->Sample(eKinetic);

  } else if (repFlag == 3) {
    if (eKinetic <= tE_of_repFlag3) {
      cosTh = theCoefficients->SampleElastic(eKinetic);
    } else {
      cosTh = theProbArray->Sample(eKinetic);
    }

  } else if (repFlag == 0) {
    cosTh = 2.*G4UniformRand() - 1.;

  } else {
    G4cout << "Unusable number for repFlag: repFlag=" << repFlag << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
                              "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
  }

  if (cosTh < -1.1) { return 0; }

  G4double phi = twopi*G4UniformRand();
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
    G4double Einit = theNeutron.GetTotalEnergy();    // Incident energy
    G4double mT = theTarget.GetMass();

    G4double ratio = mT/mN;
    G4double sqt = std::sqrt(ratio*ratio - 1.0 + cosTh*cosTh);
    G4double beta = Pinit/(mT + Einit);              // CMS beta
    G4double denom = 1. - beta*beta*cosTh*cosTh;
    G4double term1 = cosTh*(Einit*ratio + mN)/(mN*ratio + Einit);
    G4double pN = beta*mN*(term1 + sqt)/denom;

    // Get the scattered momentum and rotate it in theta and phi
    G4ThreeVector pDir = theNeutron.GetMomentum()/Pinit;
    G4double px = pN*pDir.x();
    G4double py = pN*pDir.y();
    G4double pz = pN*pDir.z();

    G4ThreeVector pcmRot;
    pcmRot.setX(px*cosTh*cosPhi - py*sinPhi + pz*sinth*cosPhi);
    pcmRot.setY(px*cosTh*sinPhi + py*cosPhi + pz*sinth*sinPhi);
    pcmRot.setZ(-px*sinth + pz*cosTh);
    theNeutron.SetMomentum(pcmRot);
    G4double eN = std::sqrt(pN*pN + mN*mN);    // Scattered neutron energy
    theNeutron.SetTotalEnergy(eN);

    // Get the scattered target momentum 
    G4ReactionProduct toLab(-1.*theTarget);
    theTarget.SetMomentum(pDir*Pinit - pcmRot);
    G4double eT = Einit - eN + mT;
    theTarget.SetTotalEnergy(eT);

    // Now back to lab frame     
    theNeutron.Lorentz(theNeutron, toLab);
    theTarget.Lorentz(theTarget, toLab);

    //111005 Protection for not producing 0 kinetic energy target
    if (theNeutron.GetKineticEnergy() <= 0)
      theNeutron.SetTotalEnergy(theNeutron.GetMass()*(1. + G4Pow::GetInstance()->powA(10, -15.65) ) );
    if (theTarget.GetKineticEnergy() <= 0) 
      theTarget.SetTotalEnergy(theTarget.GetMass()*(1. + G4Pow::GetInstance()->powA(10, -15.65) ) );

  } else if (frameFlag == 2) {
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
    pcmRot.setX(px*cosTh*cosPhi - py*sinPhi + pz*sinth*cosPhi);
    pcmRot.setY(px*cosTh*sinPhi + py*cosPhi + pz*sinth*sinPhi);
    pcmRot.setZ(-px*sinth + pz*cosTh);
    proj.setVect(pcmRot);
    targ.setVect(-pcmRot);

    // Back to lab frame
    proj.boost(-boostToCM);
    targ.boost(-boostToCM);

    theNeutron.SetMomentum(proj.vect() );
    theNeutron.SetTotalEnergy(proj.e() );

    theTarget.SetMomentum(targ.vect() );
    theTarget.SetTotalEnergy(targ.e() );

    //080904 Add Protection for very low energy (1e-6eV) scattering 
    if (theNeutron.GetKineticEnergy() <= 0) {
      theNeutron.SetTotalEnergy(theNeutron.GetMass()*(1. + G4Pow::GetInstance()->powA(10, -15.65) ) );
    }

    //080904 Add Protection for very low energy (1e-6eV) scattering 
    if (theTarget.GetKineticEnergy() <= 0) {
      theTarget.SetTotalEnergy(theTarget.GetMass()*(1. + G4Pow::GetInstance()->powA(10, -15.65) ) );
    }

  } else {
    G4cout << "Value of frameFlag (1=LAB, 2=CMS): " << frameFlag;
    throw G4HadronicException(__FILE__, __LINE__,
                   "G4ParticleHPElasticFS::ApplyYourSelf frameflag incorrect");
  }

  // Everything is now in the lab frame
  // Set energy change and momentum change
  theResult.Get()->SetEnergyChange(theNeutron.GetKineticEnergy());
  theResult.Get()->SetMomentumChange(theNeutron.GetMomentum().unit());

  // Make recoil a G4DynamicParticle
  G4DynamicParticle* theRecoil = new G4DynamicParticle;
  theRecoil->SetDefinition(G4IonTable::GetIonTable()->GetIon(static_cast<G4int>(theBaseZ),
                           static_cast<G4int>(theBaseA), 0) );
  theRecoil->SetMomentum(theTarget.GetMomentum());
  theResult.Get()->AddSecondary(theRecoil);

  // Postpone the tracking of the primary neutron
  theResult.Get()->SetStatusChange(suspend);
  return theResult.Get();
}

