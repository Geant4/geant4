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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 070606 bug fix and migrate to enable to Partial cases by T. Koi
// 080603 bug fix for Hadron Hyper News #932 by T. Koi
// 080612 bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #4,6
// 080717 bug fix of calculation of residual momentum by T. Koi
// 080801 protect negative available energy by T. Koi
//        introduce theNDLDataA,Z which has A and Z of NDL data by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
// 090514 Fix bug in IC electron emission case
//        Contribution from Chao Zhang (Chao.Zhang@usd.edu) and Dongming Mei(Dongming.Mei@usd.edu)
// 100406 "nothingWasKnownOnHadron=1" then sample mu isotropic in CM
//        add two_body_reaction
// 100909 add safty
// 101111 add safty for _nat_ data case in Binary reaction, but break conservation
// 110430 add Reaction Q value and break up flag (MF3::QI and LR)
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// June-2019 - E. Mendoza - re-build "two_body_reaction", to be used by incident charged particles
// (now isotropic emission in the CMS). Also restrict nresp use below 20 MeV (for future
// developments). Add photon emission when no data available.
//
// nresp71_m03.hh and nresp71_m02.hh are alike. The only difference between m02 and m03
// is in the total carbon cross section that is properly included in the latter.
// These data are not used in nresp71_m0*.hh.
//

#include "G4ParticleHPInelasticCompFS.hh"

#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4He3.hh"
#include "G4IonTable.hh"
#include "G4NRESP71M03.hh"
#include "G4NucleiProperties.hh"
#include "G4Nucleus.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4ParticleHPManager.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

G4ParticleHPInelasticCompFS::G4ParticleHPInelasticCompFS()
  : G4ParticleHPFinalState()
{
  QI.resize(51);
  LR.resize(51);
  for (G4int i = 0; i < 51; ++i) {
    hasXsec = true;
    theXsection[i] = nullptr;
    theEnergyDistribution[i] = nullptr;
    theAngularDistribution[i] = nullptr;
    theEnergyAngData[i] = nullptr;
    theFinalStatePhotons[i] = nullptr;
    QI[i] = 0.0;
    LR[i] = 0;
  }
}

G4ParticleHPInelasticCompFS::~G4ParticleHPInelasticCompFS()
{
  for (G4int i = 0; i < 51; ++i) {
    if (theXsection[i] != nullptr) delete theXsection[i];
    if (theEnergyDistribution[i] != nullptr) delete theEnergyDistribution[i];
    if (theAngularDistribution[i] != nullptr) delete theAngularDistribution[i];
    if (theEnergyAngData[i] != nullptr) delete theEnergyAngData[i];
    if (theFinalStatePhotons[i] != nullptr) delete theFinalStatePhotons[i];
  }
}

void G4ParticleHPInelasticCompFS::InitDistributionInitialState(G4ReactionProduct& inPart,
                                  G4ReactionProduct& aTarget, G4int it)
{
  if (theAngularDistribution[it] != nullptr) {
    theAngularDistribution[it]->SetTarget(aTarget);
    theAngularDistribution[it]->SetProjectileRP(inPart);
  }

  if (theEnergyAngData[it] != nullptr) {
    theEnergyAngData[it]->SetTarget(aTarget);
    theEnergyAngData[it]->SetProjectileRP(inPart);
  }
}

void G4ParticleHPInelasticCompFS::InitGammas(G4double AR, G4double ZR)
{
  G4int Z = G4lrint(ZR);
  G4int A = G4lrint(AR);
  std::ostringstream ost;
  ost << gammaPath << "z" << Z << ".a" << A;
  G4String aName = ost.str();
  std::ifstream from(aName, std::ios::in);

  if (!from) return;  // no data found for this isotope
  std::ifstream theGammaData(aName, std::ios::in);

  theGammas.Init(theGammaData);
}

void G4ParticleHPInelasticCompFS::Init(G4double A, G4double Z, G4int M, G4String& dirName,
                                       G4String& aFSType, G4ParticleDefinition*)
{
  gammaPath = fManager->GetNeutronHPPath() + "/Inelastic/Gammas/";
  G4String tString = dirName;
  SetA_Z(A, Z, M);
  G4bool dbool;
  G4ParticleHPDataUsed aFile =
    theNames.GetName(theBaseA, theBaseZ, M, tString, aFSType, dbool);
  SetAZMs(aFile);
  G4String filename = aFile.GetName();
#ifdef G4VERBOSE
  if (fManager->GetDEBUG())
    G4cout << " G4ParticleHPInelasticCompFS::Init FILE " << filename << G4endl;
#endif

  SetAZMs(A, Z, M, aFile);

  if (!dbool || (theBaseZ <= 2 && (theNDLDataZ != theBaseZ || theNDLDataA != theBaseA)))
  {
#ifdef G4VERBOSE
    if (fManager->GetDEBUG())
      G4cout << "Skipped = " << filename << " " << A << " " << Z << G4endl;
#endif
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }
  std::istringstream theData(std::ios::in);
  fManager->GetDataStream(filename, theData);
  if (!theData)  //"!" is a operator of ios
  {
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }
  // here we go
  G4int infoType, dataType, dummy;
  G4int sfType, it;
  hasFSData = false;
  while (theData >> infoType)  // Loop checking, 11.05.2015, T. Koi
  {
    hasFSData = true;
    theData >> dataType;
    theData >> sfType >> dummy;
    it = 50;
    if (sfType >= 600 || (sfType < 100 && sfType >= 50))
      it = sfType % 50;
    if (dataType == 3)
    {
      G4double dqi;
      G4int ilr;
      theData >> dqi >> ilr;

      QI[it] = dqi * CLHEP::eV;
      LR[it] = ilr;
      theXsection[it] = new G4ParticleHPVector;
      G4int total;
      theData >> total;
      theXsection[it]->Init(theData, total, CLHEP::eV);
    }
    else if (dataType == 4) {
      theAngularDistribution[it] = new G4ParticleHPAngular;
      theAngularDistribution[it]->Init(theData);
    }
    else if (dataType == 5) {
      theEnergyDistribution[it] = new G4ParticleHPEnergyDistribution;
      theEnergyDistribution[it]->Init(theData);
    }
    else if (dataType == 6) {
      theEnergyAngData[it] = new G4ParticleHPEnAngCorrelation(theProjectile);
      //      G4cout << this << " CompFS theEnergyAngData " << it << theEnergyAngData[it] << G4endl;
      theEnergyAngData[it]->Init(theData);
    }
    else if (dataType == 12) {
      theFinalStatePhotons[it] = new G4ParticleHPPhotonDist;
      theFinalStatePhotons[it]->InitMean(theData);
    }
    else if (dataType == 13) {
      theFinalStatePhotons[it] = new G4ParticleHPPhotonDist;
      theFinalStatePhotons[it]->InitPartials(theData, theXsection[50]);
    }
    else if (dataType == 14) {
      theFinalStatePhotons[it]->InitAngular(theData);
    }
    else if (dataType == 15) {
      theFinalStatePhotons[it]->InitEnergies(theData);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Z=" << theBaseZ << " A=" << theBaseA << " dataType=" << dataType
	 << " projectile: " << theProjectile->GetParticleName();
      G4Exception("G4ParticleHPInelasticCompFS::Init", "hadr01", JustWarning,
		  ed, "Data-type unknown");
    }
  }
}

G4int G4ParticleHPInelasticCompFS::SelectExitChannel(G4double eKinetic)
{
  G4double running[50];
  running[0] = 0;
  G4int i;
  for (i = 0; i < 50; ++i) {
    if (i != 0) running[i] = running[i - 1];
    if (theXsection[i] != nullptr) {
      running[i] += std::max(0., theXsection[i]->GetXsec(eKinetic));
    }
  }
  G4double random = G4UniformRand();
  G4double sum = running[49];
  G4int it = 50;
  if (0 != sum) {
    G4int i0;
    for (i0 = 0; i0 < 50; ++i0) {
      it = i0;
      if (random < running[i0] / sum) break;
    }
  }
  return it;
}

// n,p,d,t,he3,a
void G4ParticleHPInelasticCompFS::CompositeApply(const G4HadProjectile& theTrack,
                                                 G4ParticleDefinition* aDefinition)
{
  // prepare neutron
  if (theResult.Get() == nullptr) theResult.Put(new G4HadFinalState);
  theResult.Get()->Clear();
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile* hadProjectile = &theTrack;
  G4ReactionProduct incidReactionProduct(hadProjectile->GetDefinition());
  incidReactionProduct.SetMomentum(hadProjectile->Get4Momentum().vect());
  incidReactionProduct.SetKineticEnergy(eKinetic);

  // prepare target
  for (G4int i = 0; i < 50; ++i) {
    if (theXsection[i] != nullptr) {
      break;
    }
  }

  G4double targetMass = G4NucleiProperties::GetNuclearMass(theBaseA, theBaseZ);
#ifdef G4VERBOSE
  if (fManager->GetDEBUG())
    G4cout << "G4ParticleHPInelasticCompFS::CompositeApply A=" << theBaseA << " Z="
	   << theBaseZ << " incident " << hadProjectile->GetDefinition()->GetParticleName()
	   << G4endl;
#endif
  G4ReactionProduct theTarget;
  G4Nucleus aNucleus;
  // G4ThreeVector neuVelo =
  // (1./hadProjectile->GetDefinition()->GetPDGMass())*incidReactionProduct.GetMomentum(); theTarget
  // = aNucleus.GetBiasedThermalNucleus( targetMass/hadProjectile->GetDefinition()->GetPDGMass() ,
  // neuVelo, theTrack.GetMaterial()->GetTemperature()); G4Nucleus::GetBiasedThermalNucleus requests
  // normalization of mass and velocity in neutron mass
  G4ThreeVector neuVelo = incidReactionProduct.GetMomentum() / CLHEP::neutron_mass_c2;
  theTarget = aNucleus.GetBiasedThermalNucleus(targetMass / CLHEP::neutron_mass_c2,
                                               neuVelo, theTrack.GetMaterial()->GetTemperature());

  theTarget.SetDefinition(G4IonTable::GetIonTable()->GetIon(theBaseZ, theBaseA, 0.0));

  // prepare the residual mass
  G4double residualMass = 0;
  G4int residualZ = theBaseZ +
    G4lrint((theProjectile->GetPDGCharge() - aDefinition->GetPDGCharge())/CLHEP::eplus);
  G4int residualA = theBaseA + theProjectile->GetBaryonNumber() - aDefinition->GetBaryonNumber();
  residualMass = G4NucleiProperties::GetNuclearMass(residualA, residualZ);

  // prepare energy in target rest frame
  G4ReactionProduct boosted;
  boosted.Lorentz(incidReactionProduct, theTarget);
  eKinetic = boosted.GetKineticEnergy();

  // select exit channel for composite FS class.
  G4int it = SelectExitChannel(eKinetic);

  // E. Mendoza (2018) -- to use JENDL/AN-2005
  if (theEnergyDistribution[it] == nullptr && theAngularDistribution[it] == nullptr
      && theEnergyAngData[it] == nullptr)
  {
    if (theEnergyDistribution[50] != nullptr || theAngularDistribution[50] != nullptr
        || theEnergyAngData[50] != nullptr)
    {
      it = 50;
    }
  }

  // set target and neutron in the relevant exit channel
  InitDistributionInitialState(incidReactionProduct, theTarget, it);

  //---------------------------------------------------------------------//
  // Hook for NRESP71MODEL
  if (fManager->GetUseNRESP71Model() && eKinetic < 20 * CLHEP::MeV) {
    if (theBaseZ == 6)  // If the reaction is with Carbon...
    {
      if (theProjectile == G4Neutron::Definition()) {
        if (use_nresp71_model(aDefinition, it, theTarget, boosted)) return;
      }
    }
  }
  //---------------------------------------------------------------------//

  G4ReactionProductVector* thePhotons = nullptr;
  G4ReactionProductVector* theParticles = nullptr;
  G4ReactionProduct aHadron;
  aHadron.SetDefinition(aDefinition);  // what if only cross-sections exist ==> Na 23 11 @@@@
  G4double availableEnergy = incidReactionProduct.GetKineticEnergy()
                             + incidReactionProduct.GetMass() - aHadron.GetMass()
                             + (targetMass - residualMass);

  if (availableEnergy < 0) {
    availableEnergy = 0;
  }
  G4int nothingWasKnownOnHadron = 0;
  G4double eGamm = 0;
  G4int iLevel = -1;
  // max gamma energy and index
  G4int imaxEx = theGammas.GetNumberOfLevels() - 1;

  // without photon has it = 0
  if (50 == it) {
    // Excitation level is not determined
    aHadron.SetKineticEnergy(availableEnergy * residualMass / (aHadron.GetMass() + residualMass));

    // TK add safty 100909
    G4double p2 =
      (aHadron.GetTotalEnergy() * aHadron.GetTotalEnergy() - aHadron.GetMass() * aHadron.GetMass());
    G4double p = (p2 > 0.0) ? std::sqrt(p2) : 0.0;
    aHadron.SetMomentum(p * incidReactionProduct.GetMomentum() /
			incidReactionProduct.GetTotalMomentum());
  }
  else {
    iLevel = imaxEx;
  }

  if (theAngularDistribution[it] != nullptr)  // MF4
  {
    if (theEnergyDistribution[it] != nullptr)  // MF5
    {
      //************************************************************
      /*
            aHadron.SetKineticEnergy(theEnergyDistribution[it]->Sample(eKinetic, dummy));
            G4double eSecN = aHadron.GetKineticEnergy();
      */
      //************************************************************
      // EMendoza --> maximum allowable energy should be taken into account.
      G4double dqi = 0.0;
      if (QI[it] < 0 || 849 < QI[it])
        dqi = QI[it];  // For backword compatibility QI introduced since G4NDL3.15
      G4double MaxEne = eKinetic + dqi;
      G4double eSecN = 0.;

      G4int icounter = 0;
      G4int icounter_max = 1024;
      G4int dummy = 0;
      do {
        ++icounter;
        if (icounter > icounter_max) {
          G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of "
                 << __FILE__ << "." << G4endl;
          break;
        }
        eSecN = theEnergyDistribution[it]->Sample(eKinetic, dummy);
      } while (eSecN > MaxEne);  // Loop checking, 11.05.2015, T. Koi
      aHadron.SetKineticEnergy(eSecN);
      //************************************************************
      eGamm = eKinetic - eSecN;
      for (iLevel = imaxEx; iLevel >= 0; --iLevel) {
        if (theGammas.GetLevelEnergy(iLevel) < eGamm) break;
      }
      if (iLevel < imaxEx && iLevel >= 0) {
        if (G4UniformRand() > 0.5) {
          ++iLevel;
        }
      }
    }
    else {
      G4double eExcitation = 0;
      for (iLevel = imaxEx; iLevel >= 0; --iLevel) {
        if (theGammas.GetLevelEnergy(iLevel) < eKinetic) break;
      }

      // Use QI value for calculating excitation energy of residual.
      G4bool useQI = false;
      G4double dqi = QI[it];
      if (dqi < 0 || 849 < dqi) useQI = true;  // Former libraries do not have values in this range

      if (useQI) {
        eExcitation = std::max(0., QI[0] - QI[it]);  // Bug fix 2333

        // Re-evaluate iLevel based on this eExcitation
        iLevel = 0;
        G4bool find = false;
        const G4double level_tolerance = 1.0 * CLHEP::keV;

        // VI: the first level is ground
        if (0 < imaxEx) {
          for (iLevel = 1; iLevel <= imaxEx; ++iLevel) {
            G4double elevel = theGammas.GetLevelEnergy(iLevel);
            if (std::abs(eExcitation - elevel) < level_tolerance) {
              find = true;
              break;
            }
            if (eExcitation < elevel) {
              find = true;
              iLevel = std::max(iLevel - 1, 0);
              break;
            }
          }

          // If proper level cannot be found, use the maximum level
          if (!find) iLevel = imaxEx;
        }
      }

      if (fManager->GetDEBUG() && eKinetic - eExcitation < 0) {
        throw G4HadronicException(
          __FILE__, __LINE__,
          "SEVERE: InelasticCompFS: Consistency of data not good enough, please file report");
      }
      if (eKinetic - eExcitation < 0) eExcitation = 0;
      if (iLevel != -1) aHadron.SetKineticEnergy(eKinetic - eExcitation);
    }
    theAngularDistribution[it]->SampleAndUpdate(aHadron);

    if (theFinalStatePhotons[it] == nullptr) {
      thePhotons = theGammas.GetDecayGammas(iLevel);
      eGamm -= theGammas.GetLevelEnergy(iLevel);
    }
  }
  else if (theEnergyAngData[it] != nullptr)  // MF6
  {
    theParticles = theEnergyAngData[it]->Sample(eKinetic);

    // Adjust A and Z in the case of miss much between selected data and target nucleus
    if (theParticles != nullptr) {
      G4int sumA = 0;
      G4int sumZ = 0;
      G4int maxA = 0;
      G4int jAtMaxA = 0;
      for (G4int j = 0; j != (G4int)theParticles->size(); ++j) {
	auto ptr = theParticles->at(j);
	G4int barnum = ptr->GetDefinition()->GetBaryonNumber();
        if (barnum > maxA) {
          maxA = barnum;
          jAtMaxA = j;
        }
        sumA += barnum;
        sumZ += G4lrint(ptr->GetDefinition()->GetPDGCharge()/CLHEP::eplus);
      }
      G4int dA = theBaseA + hadProjectile->GetDefinition()->GetBaryonNumber() - sumA;
      G4int dZ = theBaseZ +
	G4lrint(hadProjectile->GetDefinition()->GetPDGCharge()/CLHEP::eplus) - sumZ;
      if (dA < 0 || dZ < 0) {
        G4int newA = theParticles->at(jAtMaxA)->GetDefinition()->GetBaryonNumber() + dA;
        G4int newZ = 
	  G4lrint(theParticles->at(jAtMaxA)->GetDefinition()->GetPDGCharge()/CLHEP::eplus) + dZ;
        G4ParticleDefinition* pd = ionTable->GetIon(newZ, newA);
        theParticles->at(jAtMaxA)->SetDefinition(pd);
      }
    }
  }
  else {
    // @@@ what to do, if we have photon data, but no info on the hadron itself
    nothingWasKnownOnHadron = 1;
  }

  if (theFinalStatePhotons[it] != nullptr) {
    // the photon distributions are in the Nucleus rest frame.
    // TK residual rest frame
    G4ReactionProduct boosted_tmp;
    boosted_tmp.Lorentz(incidReactionProduct, theTarget);
    G4double anEnergy = boosted_tmp.GetKineticEnergy();
    thePhotons = theFinalStatePhotons[it]->GetPhotons(anEnergy);
    G4double aBaseEnergy = theFinalStatePhotons[it]->GetLevelEnergy();
    G4double testEnergy = 0;
    if (thePhotons != nullptr && !thePhotons->empty()) {
      aBaseEnergy -= (*thePhotons)[0]->GetTotalEnergy();
    }
    if (theFinalStatePhotons[it]->NeedsCascade()) {
      while (aBaseEnergy > 0.01 * CLHEP::keV)  // Loop checking, 11.05.2015, T. Koi
      {
        // cascade down the levels
        G4bool foundMatchingLevel = false;
        G4int closest = 2;
        G4double deltaEold = -1;
        for (G4int j = 1; j < it; ++j) {
          if (theFinalStatePhotons[j] != nullptr) {
            testEnergy = theFinalStatePhotons[j]->GetLevelEnergy();
          }
          else {
            testEnergy = 0;
          }
          G4double deltaE = std::abs(testEnergy - aBaseEnergy);
          if (deltaE < 0.1 * CLHEP::keV) {
            G4ReactionProductVector* theNext = theFinalStatePhotons[j]->GetPhotons(anEnergy);
            if (thePhotons != nullptr) thePhotons->push_back(theNext->operator[](0));
            aBaseEnergy = testEnergy - theNext->operator[](0)->GetTotalEnergy();
            delete theNext;
            foundMatchingLevel = true;
            break;  // ===>
          }
          if (theFinalStatePhotons[j] != nullptr && (deltaE < deltaEold || deltaEold < 0.)) {
            closest = j;
            deltaEold = deltaE;
          }
        }  // <=== the break goes here.
        if (!foundMatchingLevel) {
          G4ReactionProductVector* theNext = theFinalStatePhotons[closest]->GetPhotons(anEnergy);
          if (thePhotons != nullptr) thePhotons->push_back(theNext->operator[](0));
          aBaseEnergy = aBaseEnergy - theNext->operator[](0)->GetTotalEnergy();
          delete theNext;
        }
      }
    }
  }

  if (thePhotons != nullptr) {
    for (auto const & p : *thePhotons) {
      // back to lab
      p->Lorentz(*p, -1. * theTarget);
    }
  }
  if (nothingWasKnownOnHadron != 0) {
    // In this case, hadron should be isotropic in CM
    // Next 12 lines are Emilio's replacement
    // G4double QM=(incidReactionProduct.GetMass()+targetMass)-(aHadron.GetMass()+residualMass);
    // G4double eExcitation = QM-QI[it];
    // G4double eExcitation = QI[0] - QI[it];  // Fix of bug #1838
    // if(eExcitation<20*CLHEP::keV){eExcitation=0;}

    G4double eExcitation = std::max(0., QI[0] - QI[it]);  // Fix of bug #2333

    two_body_reaction(&incidReactionProduct, &theTarget, &aHadron, eExcitation);
    if (thePhotons == nullptr && eExcitation > 0) {
      for (iLevel = imaxEx; iLevel >= 0; --iLevel) {
        if (theGammas.GetLevelEnergy(iLevel) < eExcitation + 5 * keV) break;  // 5 keV tolerance
      }
      thePhotons = theGammas.GetDecayGammas(iLevel);
    }
  }

  // fill the result
  // Beware - the recoil is not necessarily in the particles...
  // Can be calculated from momentum conservation?
  // The idea is that the particles ar emitted forst, and the gammas only once the
  // recoil is on the residual; assumption is that gammas do not contribute to
  // the recoil.
  // This needs more design @@@

  G4bool needsSeparateRecoil = false;
  G4int totalBaryonNumber = 0;
  G4int totalCharge = 0;
  G4ThreeVector totalMomentum(0);
  if (theParticles != nullptr) {
    const G4ParticleDefinition* aDef;
    for (std::size_t ii0 = 0; ii0 < theParticles->size(); ++ii0) {
      aDef = (*theParticles)[ii0]->GetDefinition();
      totalBaryonNumber += aDef->GetBaryonNumber();
      totalCharge += G4lrint(aDef->GetPDGCharge()/CLHEP::eplus);
      totalMomentum += (*theParticles)[ii0]->GetMomentum();
    }
    if (totalBaryonNumber
        != theBaseA +  hadProjectile->GetDefinition()->GetBaryonNumber())
    {
      needsSeparateRecoil = true;
      residualA = theBaseA + hadProjectile->GetDefinition()->GetBaryonNumber()
	- totalBaryonNumber;
      residualZ = theBaseZ +
	G4lrint((hadProjectile->GetDefinition()->GetPDGCharge() - totalCharge)/CLHEP::eplus);
    }
  }

  std::size_t nPhotons = 0;
  if (thePhotons != nullptr) {
    nPhotons = thePhotons->size();
  }

  G4DynamicParticle* theSec;

  if (theParticles == nullptr) {
    theSec = new G4DynamicParticle;
    theSec->SetDefinition(aHadron.GetDefinition());
    theSec->SetMomentum(aHadron.GetMomentum());
    theResult.Get()->AddSecondary(theSec, secID);
#ifdef G4VERBOSE
    if (fManager->GetDEBUG())
      G4cout << " G4ParticleHPInelasticCompFS::BaseApply  add secondary1 "
             << theSec->GetParticleDefinition()->GetParticleName()
             << " E= " << theSec->GetKineticEnergy() << " NSECO "
             << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif

    aHadron.Lorentz(aHadron, theTarget);
    G4ReactionProduct theResidual;
    theResidual.SetDefinition(ionTable->GetIon(residualZ, residualA, 0));
    theResidual.SetKineticEnergy(aHadron.GetKineticEnergy() * aHadron.GetMass()
                                 / theResidual.GetMass());

    // 080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #6
    // theResidual.SetMomentum(-1.*aHadron.GetMomentum());
    G4ThreeVector incidentNeutronMomentum = incidReactionProduct.GetMomentum();
    theResidual.SetMomentum(incidentNeutronMomentum - aHadron.GetMomentum());

    theResidual.Lorentz(theResidual, -1. * theTarget);
    G4ThreeVector totalPhotonMomentum(0, 0, 0);
    if (thePhotons != nullptr) {
      for (std::size_t i = 0; i < nPhotons; ++i) {
        totalPhotonMomentum += (*thePhotons)[i]->GetMomentum();
      }
    }
    theSec = new G4DynamicParticle;
    theSec->SetDefinition(theResidual.GetDefinition());
    theSec->SetMomentum(theResidual.GetMomentum() - totalPhotonMomentum);
    theResult.Get()->AddSecondary(theSec, secID);
#ifdef G4VERBOSE
    if (fManager->GetDEBUG())
      G4cout << this << " G4ParticleHPInelasticCompFS::BaseApply add secondary2 "
             << theSec->GetParticleDefinition()->GetParticleName()
             << " E= " << theSec->GetKineticEnergy() << " NSECO "
             << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
  }
  else {
    for (std::size_t i0 = 0; i0 < theParticles->size(); ++i0) {
      theSec = new G4DynamicParticle;
      theSec->SetDefinition((*theParticles)[i0]->GetDefinition());
      theSec->SetMomentum((*theParticles)[i0]->GetMomentum());
      theResult.Get()->AddSecondary(theSec, secID);
#ifdef G4VERBOSE
      if (fManager->GetDEBUG())
        G4cout << " G4ParticleHPInelasticCompFS::BaseApply add secondary3 "
               << theSec->GetParticleDefinition()->GetParticleName()
               << " E= " << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
      delete (*theParticles)[i0];
    }
    delete theParticles;
    if (needsSeparateRecoil && residualZ != 0) {
      G4ReactionProduct theResidual;
      theResidual.SetDefinition(ionTable->GetIon(residualZ, residualA, 0));
      G4double resiualKineticEnergy = theResidual.GetMass() * theResidual.GetMass();
      resiualKineticEnergy += totalMomentum * totalMomentum;
      resiualKineticEnergy = std::sqrt(resiualKineticEnergy) - theResidual.GetMass();
      theResidual.SetKineticEnergy(resiualKineticEnergy);

      // 080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #4
      // theResidual.SetMomentum(-1.*totalMomentum);
      // G4ThreeVector incidentNeutronMomentum = incidReactionProduct.GetMomentum();
      // theResidual.SetMomentum(incidentNeutronMomentum - aHadron.GetMomentum());
      // 080717 TK Comment still do NOT include photon's mometum which produce by thePhotons
      theResidual.SetMomentum(incidReactionProduct.GetMomentum() + theTarget.GetMomentum()
                              - totalMomentum);

      theSec = new G4DynamicParticle;
      theSec->SetDefinition(theResidual.GetDefinition());
      theSec->SetMomentum(theResidual.GetMomentum());
      theResult.Get()->AddSecondary(theSec, secID);
#ifdef G4VERBOSE
      if (fManager->GetDEBUG())
        G4cout << " G4ParticleHPInelasticCompFS::BaseApply add secondary4 "
               << theSec->GetParticleDefinition()->GetParticleName()
               << " E= " << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
    }
  }
  if (thePhotons != nullptr) {
    for (std::size_t i = 0; i < nPhotons; ++i) {
      theSec = new G4DynamicParticle;
      // Bug reported Chao Zhang (Chao.Zhang@usd.edu), Dongming Mei(Dongming.Mei@usd.edu) Feb. 25,
      // 2009 theSec->SetDefinition(G4Gamma::Gamma());
      theSec->SetDefinition((*thePhotons)[i]->GetDefinition());
      // But never cause real effect at least with G4NDL3.13 TK
      theSec->SetMomentum((*thePhotons)[i]->GetMomentum());
      theResult.Get()->AddSecondary(theSec, secID);
#ifdef G4VERBOSE
      if (fManager->GetDEBUG())
        G4cout << " G4ParticleHPInelasticCompFS::BaseApply add secondary5 "
               << theSec->GetParticleDefinition()->GetParticleName()
               << " E= " << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif

      delete thePhotons->operator[](i);
    }
    // some garbage collection
    delete thePhotons;
  }

  G4ParticleDefinition* targ_pd = ionTable->GetIon(theBaseZ, theBaseA, 0.0);
  G4LorentzVector targ_4p_lab(
    theTarget.GetMomentum(),
    std::sqrt(targ_pd->GetPDGMass() * targ_pd->GetPDGMass() + theTarget.GetMomentum().mag2()));
  G4LorentzVector proj_4p_lab = theTrack.Get4Momentum();
  G4LorentzVector init_4p_lab = proj_4p_lab + targ_4p_lab;
  adjust_final_state(init_4p_lab);

  // clean up the primary neutron
  theResult.Get()->SetStatusChange(stopAndKill);
}

// Re-implemented by E. Mendoza (2019). Isotropic emission in the CMS:
//  proj: projectile in target-rest-frame (input)
//  targ: target in target-rest-frame (input)
//  product: secondary particle in target-rest-frame (output)
//  resExcitationEnergy: excitation energy of the residual nucleus
//
void G4ParticleHPInelasticCompFS::two_body_reaction(G4ReactionProduct* proj,
                                                    G4ReactionProduct* targ,
                                                    G4ReactionProduct* product,
                                                    G4double resExcitationEnergy)
{
  // CMS system:
  G4ReactionProduct theCMS = *proj + *targ;

  // Residual definition:
  G4int resZ = G4lrint((proj->GetDefinition()->GetPDGCharge() + targ->GetDefinition()->GetPDGCharge()
			- product->GetDefinition()->GetPDGCharge())/CLHEP::eplus);
  G4int resA = proj->GetDefinition()->GetBaryonNumber() + targ->GetDefinition()->GetBaryonNumber()
               - product->GetDefinition()->GetBaryonNumber();
  G4ReactionProduct theResidual;
  theResidual.SetDefinition(ionTable->GetIon(resZ, resA, 0.0));

  // CMS system:
  G4ReactionProduct theCMSproj;
  G4ReactionProduct theCMStarg;
  theCMSproj.Lorentz(*proj, theCMS);
  theCMStarg.Lorentz(*targ, theCMS);
  // final Momentum in the CMS:
  G4double totE = std::sqrt(theCMSproj.GetMass() * theCMSproj.GetMass()
                            + theCMSproj.GetTotalMomentum() * theCMSproj.GetTotalMomentum())
                  + std::sqrt(theCMStarg.GetMass() * theCMStarg.GetMass()
                              + theCMStarg.GetTotalMomentum() * theCMStarg.GetTotalMomentum());
  G4double prodmass = product->GetMass();
  G4double resmass = theResidual.GetMass() + resExcitationEnergy;
  G4double fmomsquared = (totE * totE - (prodmass - resmass) * (prodmass - resmass)) *
    (totE * totE - (prodmass + resmass) * (prodmass + resmass)) / (4.*totE*totE);
  G4double fmom = (fmomsquared > 0) ? std::sqrt(fmomsquared) : 0.0;

  // random (isotropic direction):
  product->SetMomentum(fmom * G4RandomDirection());
  product->SetTotalEnergy(std::sqrt(prodmass * prodmass + fmom * fmom));  // CMS
  // Back to the LAB system:
  product->Lorentz(*product, -1. * theCMS);
}

G4bool G4ParticleHPInelasticCompFS::use_nresp71_model(const G4ParticleDefinition* aDefinition,
                                                      const G4int itt,
                                                      const G4ReactionProduct& theTarget,
                                                      G4ReactionProduct& boosted)
{
  if (aDefinition == G4Neutron::Definition())  // If the outgoing particle is a neutron...
  {
    // LR: flag LR in ENDF. It indicates whether there is breakup of the residual nucleus or not.
    // it: exit channel (index of the carbon excited state)

    // Added by A. R. Garcia (CIEMAT) to include the physics of C(N,N'3A) reactions from NRESP71.

    if (LR[itt] > 0) // If there is breakup of the residual nucleus LR(flag LR in ENDF)>0 (i.e. Z=6
                     // MT=52-91 (it=MT-50)).
    {
      // Defining carbon as the target in the reference frame at rest.
      G4ReactionProduct theCarbon(theTarget);

      theCarbon.SetMomentum(G4ThreeVector());
      theCarbon.SetKineticEnergy(0.);

      // Creating four reaction products.
      G4ReactionProduct theProds[4];

      // Applying C(N,N'3A) reaction mechanisms in the target rest frame.
      if (itt == 41) {
        // QI=QM=-7.275 MeV for C-0(N,N')C-C(3A) in ENDF/B-VII.1.
        // This is not the value of the QI of the first step according
        // to the model. So we don't take it. Instead, we set the one
        // we have calculated: QI=(mn+m12C)-(ma+m9Be+Ex9Be)=-8.130 MeV.
        nresp71_model.ApplyMechanismI_NBeA2A(boosted, theCarbon, theProds, -8.130 /*QI[it]*/);
        // N+C --> A[0]+9BE* | 9BE* --> N[1]+8BE | 8BE --> 2*A[2,3].
      }
      else {
        nresp71_model.ApplyMechanismII_ACN2A(boosted, theCarbon, theProds, QI[itt]);
        // N+C --> N'[0]+C* | C* --> A[1]+8BE | 8BE --> 2*A[2,3].
      }

      // Returning to the reference frame where the target was in motion.
      for (auto& theProd : theProds) {
        theProd.Lorentz(theProd, -1. * theTarget);
        theResult.Get()->AddSecondary(
          new G4DynamicParticle(theProd.GetDefinition(), theProd.GetMomentum()), secID);
      }

      // Killing the primary neutron.
      theResult.Get()->SetStatusChange(stopAndKill);

      return true;
    }
  }
  else if (aDefinition == G4Alpha::Definition())  // If the outgoing particle is an alpha, ...
  {
    // Added by A. R. Garcia (CIEMAT) to include the physics of C(N,A)9BE reactions from NRESP71.

    if (LR[itt] == 0) // If Z=6, an alpha particle is emitted and there is no breakup of the
                      // residual nucleus LR(flag LR in ENDF)==0.
    {
      // Defining carbon as the target in the reference frame at rest.
      G4ReactionProduct theCarbon(theTarget);
      theCarbon.SetMomentum(G4ThreeVector());
      theCarbon.SetKineticEnergy(0.);

      // Creating four reaction products.
      G4ReactionProduct theProds[2];

      // Applying C(N,A)9BE reaction mechanism.
      nresp71_model.ApplyMechanismABE(boosted, theCarbon, theProds);
      // N+C --> A[0]+9BE[1].

      for (auto& theProd : theProds) {
        // Returning to the system of reference where the target was in motion.
        theProd.Lorentz(theProd, -1. * theTarget);
        theResult.Get()->AddSecondary(
          new G4DynamicParticle(theProd.GetDefinition(), theProd.GetMomentum()), secID);
      }

      // Killing the primary neutron.
      theResult.Get()->SetStatusChange(stopAndKill);
      return true;
    }
    G4Exception("G4ParticleHPInelasticCompFS::CompositeApply()", "G4ParticleInelasticCompFS.cc",
                FatalException, "Alpha production with LR!=0.");
  }
  return false;
}
