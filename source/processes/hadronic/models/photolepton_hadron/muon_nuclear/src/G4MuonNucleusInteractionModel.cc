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
// $Id: G4MuonNucleusInteractionModel.cc,v 1.6 2006-06-29 20:57:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
// G4MuonNucleusInteractionModel.cc
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#include "G4MuonNucleusInteractionModel.hh"


//-----------------------------------------------------------------------------
  G4MuonNucleusInteractionModel::G4MuonNucleusInteractionModel()
    : G4LeptonHadronInteractionModel()
//-----------------------------------------------------------------------------
  {
    // build the physics vector
    Nbin = 90;
    kEmin = 1.0e-5*GeV;
    kEmax = 1.0e+4*GeV;
    cascadeModelMarginalEnergy = 25.0*GeV;
    theCoefficientVector = new G4PhysicsLogVector(kEmin, kEmax, Nbin);
    makePhysicsVector();

    // construct variables
    LEPionMinusInelastic  = new G4LEPionMinusInelastic;
    LEPionPlusInelastic   = new G4LEPionPlusInelastic;
    HEPionMinusInelastic  = new G4HEPionMinusInelastic;
    HEPionPlusInelastic   = new G4HEPionPlusInelastic;
  }


//-----------------------------------------------------------------------------
  G4MuonNucleusInteractionModel::~G4MuonNucleusInteractionModel()
//-----------------------------------------------------------------------------
  {
    delete LEPionMinusInelastic;
    delete LEPionPlusInelastic;
    delete HEPionMinusInelastic;
    delete HEPionPlusInelastic;

    delete theCoefficientVector;
  }


//-----------------------------------------------------------------------------
  G4double G4MuonNucleusInteractionModel::tetal[35] = {
//-----------------------------------------------------------------------------
    1.0000000,  0.9999995,  0.9999990,  0.9999981,  0.9999962,
    0.9999943,  0.9999905,  0.9999847,  0.9999752,  0.9999599,
    0.9999352,  0.9998951,  0.9998302,  0.9997253,  0.9995556,
    0.9992810,  0.9988368,  0.9981183,  0.9969561,  0.9950773,
    0.9920409,  0.9871377,  0.9792297,  0.9665010,  0.9460785,
    0.9134827,  0.8618938,  0.7813507,  0.6583430,  0.4770452,
    0.2247237, -0.0955139, -0.4461272, -0.7495149, -0.9900000
  };


//-----------------------------------------------------------------------------
  G4double G4MuonNucleusInteractionModel::xeml[23] = {
//-----------------------------------------------------------------------------
    1.000,  0.998,  0.997,  0.996,  0.995,  0.994, 0.992, 0.990,
    0.970,  0.950,  0.920,  0.890,  0.850,  0.800, 0.750, 0.700,
    0.600,  0.500,  0.400,  0.300,  0.200,  0.100, 0.050
  };


//-----------------------------------------------------------------------------
  G4double G4MuonNucleusInteractionModel::computeMicroscopicCrossSection
    (const G4Track &muonTrack)
//-----------------------------------------------------------------------------
  {
    const G4DynamicParticle *muonDynamics = muonTrack.GetDynamicParticle();
    G4double kineticEnergy = muonDynamics->GetKineticEnergy();
    G4double muonMass      = muonDynamics->GetDefinition()->GetPDGMass();

    G4double totalEnergy = kineticEnergy + muonMass;

    G4double microscopicCrossSection;
    if(totalEnergy <= 30.0*GeV) {
      microscopicCrossSection
        = 0.0003*millibarn;
    } else {
      microscopicCrossSection
        = 0.0003*std::pow((totalEnergy/(30.0*GeV)), 0.25)*millibarn;
    }

    return microscopicCrossSection;
  }


//-----------------------------------------------------------------------------
  void G4MuonNucleusInteractionModel::makePhysicsVector()
//-----------------------------------------------------------------------------
  {
    G4double Ei, Ef;  // initial and final energy of incident muon;
    G4double muonMass = G4MuonMinus::MuonMinus()->GetPDGMass();

    for (G4int i=0; i<=(Nbin-1); i++)
    {
      G4double totalCrossSection = 0.0;
      Ei = theCoefficientVector->GetLowEdgeEnergy(i) + muonMass;
      for (G4int j=1; j<=34; j++)
      {
        cosTheta = 0.5 * (tetal[j] + tetal[j-1]);
        for (G4int k=1; k<=22; k++)
        {
          Ef = 0.5 * Ei * (xeml[k]+xeml[k-1]);
          G4double dsigma = computeDifferentialCrossSection(Ei,Ef,cosTheta);
          totalCrossSection = totalCrossSection 
            + Ei * (tetal[j-1]-tetal[j])*(xeml[k-1]-xeml[k]) * dsigma;
        }
      }
      theCoefficientVector->PutValue(i, totalCrossSection);
    }
  }


//-----------------------------------------------------------------------------
  G4VParticleChange* G4MuonNucleusInteractionModel::applyInteractionModel
    (const G4Track &muonTrack, G4Nucleus &targetNucleus )
//-----------------------------------------------------------------------------
  {
    G4int icos=0, ie1=0;
    G4double E1=0., P1=0.;
    G4double rndm[3];
    G4bool isOutRange;

    // Initialization
    aParticleChange.Initialize(muonTrack);

    const G4DynamicParticle *muonDynamics = muonTrack.GetDynamicParticle();
    G4double kineticEnergy = muonDynamics->GetKineticEnergy();
    G4double totalMomentum = muonDynamics->GetTotalMomentum();
    G4double totalEnergy   = muonDynamics->GetTotalEnergy();
    G4double muonMass      = muonDynamics->GetDefinition()->GetPDGMass();


    G4double W2 = 0.0;  G4int W2try = 0;
    while (W2 <= 0.0)
    {
      G4double totalCrossSection = 0.0;
      G4bool interpolated = false;
      G4double fRndm = G4UniformRand();
      G4double Hmax 
        = theCoefficientVector->GetValue(kineticEnergy, isOutRange);
      for (G4int i=1; i<=34; i++)
      {
        cosTheta = 0.5 * (tetal[i] + tetal[i-1]);
        for (G4int j=1; j<=22; j++)
        {
          E1 = 0.5 * totalEnergy * (xeml[j]+xeml[j-1]);
          G4double dsigma 
            = computeDifferentialCrossSection(totalEnergy, E1, cosTheta);
          totalCrossSection = totalCrossSection 
             + totalEnergy*(tetal[i-1]-tetal[i])*(xeml[j-1]-xeml[j])*dsigma; 

          if((fRndm*Hmax)<totalCrossSection) {
            interpolated = true;
            icos = i;  ie1 = j;
            break;
          }
        }

        if(interpolated) {
          // calculate energy, momentum and angle of outgoing muon
          CLHEP::RandFlat::shootArray(3, rndm);
          G4double theta = std::acos(tetal[icos-1]) 
                          + rndm[0]*(std::acos(tetal[icos])-std::acos(tetal[icos-1]));
          cosTheta = std::cos(theta);
          E1 = (xeml[ie1] + rndm[1]*(xeml[ie1-1]-xeml[ie1])) * totalEnergy;
            if(E1<muonMass) E1 = muonMass + 0.0001*GeV;
          P1 = std::sqrt(std::abs(E1*E1-muonMass*muonMass));

          // invariant mass of final hadron state must be greater than zero
          W2 = proton_mass_c2*proton_mass_c2
              +2.0*proton_mass_c2*(totalEnergy-E1)
              -2.0*(totalEnergy*E1-totalMomentum*P1*cosTheta-muonMass*muonMass);

          break;
        }
      }

      W2try++;
      if (W2try>100) return &aParticleChange;
    }

    // calculate momentum of outgoing muon / pion(photon)
    G4double sinTheta = std::sqrt(std::abs(1.0 - cosTheta*cosTheta));
    G4double phi = rndm[2]*twopi;
    G4ThreeVector muonDirection(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
    G4ThreeVector muonDirectionInit = muonTrack.GetMomentumDirection();
    muonDirection.rotateUz(muonDirectionInit);

    G4ParticleMomentum pionMomentum 
      = muonDynamics->GetMomentum() - P1*muonDirection;
    G4double muonKineticEnergy 
      = std::sqrt(P1*P1 + muonMass*muonMass) - muonMass;
    aParticleChange.ProposeMomentumDirection(muonDirection);
    aParticleChange.ProposeEnergy(muonKineticEnergy);
    aParticleChange.ProposeTrackStatus(fAlive);


    // virtual photon is exchanged with a pion of same Q2
    // select pi+/pi- randomly and generate pion track
    G4ParticleDefinition* pdPion;
    if(CLHEP::RandBit::shootBit())
      pdPion = G4PionMinus::PionMinusDefinition();
    else
      pdPion = G4PionPlus::PionPlusDefinition();

    G4DynamicParticle* pionDynamics
        = new G4DynamicParticle(pdPion, pionMomentum);

    G4Track* pionTrack = new G4Track(pionDynamics,
                                     muonTrack.GetGlobalTime(), 
                                     muonTrack.GetPosition() );
    pionTrack->SetStep(muonTrack.GetStep());

    // Invoke pion-nucleus inelastic process
    invokePionNucleus(*pionTrack, targetNucleus);

    // Termination
    delete pionTrack;

    return &aParticleChange;
  }


//-----------------------------------------------------------------------------
  void G4MuonNucleusInteractionModel::invokePionNucleus
    (const G4Track &pionTrack, G4Nucleus &targetNucleus )
//-----------------------------------------------------------------------------
  {
    // force interaction of pion with target nucleus
    G4double pionKineticEnergy = pionTrack.GetKineticEnergy();
    if(pionTrack.GetDefinition()->GetParticleName() == "pi-") {
      if(pionKineticEnergy <= cascadeModelMarginalEnergy)
        pionChange 
          = LEPionMinusInelastic->ApplyYourself(pionTrack, targetNucleus);
      else
        pionChange 
          = HEPionMinusInelastic->ApplyYourself(pionTrack, targetNucleus);
    } else if(pionTrack.GetDefinition()->GetParticleName() == "pi+") {
      if(pionKineticEnergy <= cascadeModelMarginalEnergy)
        pionChange 
          = LEPionPlusInelastic->ApplyYourself(pionTrack, targetNucleus);
      else
        pionChange 
          = HEPionPlusInelastic->ApplyYourself(pionTrack, targetNucleus);
    }


    // add local energy deposit
    G4double localEnergyDeposited = 0.0;
    localEnergyDeposited = pionChange->GetLocalEnergyDeposit();
    aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposited);


    // register secondary particles
    G4int numSecondaries = pionChange->GetNumberOfSecondaries();
    aParticleChange.SetNumberOfSecondaries(numSecondaries);

    G4ParticleMomentum secondaryMomentum = G4ThreeVector(0.,0.,0.);
    for(G4int iS=0; iS<=(numSecondaries-1); iS++) {
      secondaryMomentum 
        = secondaryMomentum + pionChange->GetSecondary(iS)->GetParticle()->GetMomentum();
      aParticleChange.AddSecondary(pionChange->GetSecondary(iS)->GetParticle());
    }
    pionChange->Clear();

    return;
  }


//-----------------------------------------------------------------------------
  G4double G4MuonNucleusInteractionModel::computeDifferentialCrossSection
    (G4double initialEnergy, G4double finalEnergy, G4double aCosTheta)
//-----------------------------------------------------------------------------
  {
    G4double muonMass = G4MuonMinus::MuonMinus()->GetPDGMass();

    if(finalEnergy < muonMass) return(0.0);
    if(aCosTheta >= 1.0) return DBL_MAX;

    G4double initialMomentum 
      = std::sqrt(initialEnergy*initialEnergy - muonMass*muonMass);
    G4double finalMomentum 
      = std::sqrt(finalEnergy*finalEnergy - muonMass*muonMass);


    // calculate momentum transfer (Q2)
    // and invariant mass of final state of hadrons (W2)
    G4double Q2 
      = 2.0*(initialEnergy*finalEnergy
             -initialMomentum*finalMomentum*aCosTheta-muonMass*muonMass);
    if(Q2 < 0.0) return(0.0);

    G4double W2
      = proton_mass_c2*proton_mass_c2
       +2.0*proton_mass_c2*(initialEnergy-finalEnergy)-Q2;
    if(W2 < 0.0) return(0.0);


    // calculate factors
    //   Nu      : energy transfer
    //   K       : incident flux of photon
    //   Epsilon : virtual photon polarization
    G4double fNu = initialEnergy-finalEnergy;
    G4double fK = fNu+Q2/(2.0*fNu);
    G4double fEpsilon 
      = 1.0/(1.0+2.0*((1.0-aCosTheta)/(1.0+aCosTheta))*(Q2+fNu*fNu)/Q2);
    if(fEpsilon > 1.0) return DBL_MAX;


    // calculate photoabsorption cross sections
    //   fGamma  : flux of transverse photons
    //   sigma_t : for transverse photons
    //   sigma_l : for longitudinal photons
    G4double fGamma 
      = fine_structure_const*fK*finalEnergy
        / (pi*Q2*initialEnergy*(1.0 - fEpsilon));
    G4double sigma_t = 0.12*millibarn;
    G4double sigma_l = 0.3*(1.0-Q2/(1.868*GeV*fNu))*sigma_t;
    if(sigma_l < 0.) sigma_l = 0.;

    return fGamma*(sigma_t+fEpsilon*sigma_l);
  }
