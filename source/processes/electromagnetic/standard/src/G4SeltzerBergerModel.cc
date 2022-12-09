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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SeltzerBergerModel
//
// Author:        Vladimir Ivanchenko use inheritance from Andreas Schaelicke
//                base class implementing ultra relativistic bremsstrahlung
//                model
//
// Creation date: 04.10.2011
//
// Modifications:
//
// 24.07.2018 Introduced possibility to use sampling tables to sample the
//            emitted photon energy (instead of using rejectio) from the 
//            Seltzer-Berger scalled DCS for bremsstrahlung photon emission. 
//            Using these sampling tables option gives faster(30-70%) final 
//            state generation than the original rejection but takes some 
//            extra memory (+ ~6MB in the case of the full CMS detector). 
//            (M Novak)
//
// -------------------------------------------------------------------
//

#include "G4SeltzerBergerModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4SBBremTable.hh"
#include "G4ModifiedTsai.hh"

#include "G4EmParameters.hh"
#include "G4ProductionCutsTable.hh"

#include "G4Physics2DVector.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4AutoLock.hh"

#include "G4ios.hh"

#include <fstream>
#include <iomanip>
#include <sstream>

G4Physics2DVector* G4SeltzerBergerModel::gSBDCSData[]     = { nullptr };
G4SBBremTable*     G4SeltzerBergerModel::gSBSamplingTable =   nullptr;
G4double           G4SeltzerBergerModel::gYLimitData[]    = { 0.0 };
G4String           G4SeltzerBergerModel::gDataDirectory   = "";

namespace
{
  G4Mutex theSBMutex = G4MUTEX_INITIALIZER;
}

static const G4double kMC2   = CLHEP::electron_mass_c2;
static const G4double kAlpha = CLHEP::twopi*CLHEP::fine_structure_const;

G4SeltzerBergerModel::G4SeltzerBergerModel(const G4ParticleDefinition* p,
                                           const G4String& nam)
: G4eBremsstrahlungRelModel(p,nam), fIsUseBicubicInterpolation(false),
  fIsUseSamplingTables(true), fNumWarnings(0), fIndx(0), fIndy(0)
{
  fLowestKinEnergy = 1.0*keV;
  SetLowEnergyLimit(fLowestKinEnergy);
  SetLPMFlag(false);
  SetAngularDistribution(new G4ModifiedTsai());
}

G4SeltzerBergerModel::~G4SeltzerBergerModel()
{
  // delete SB-DCS data per Z
  if (IsMaster()) {
    for (std::size_t iz = 0; iz < gMaxZet; ++iz) {
      if (gSBDCSData[iz]) {
        delete gSBDCSData[iz];
        gSBDCSData[iz] = nullptr;
      }
    }
    if (gSBSamplingTable) {
      delete gSBSamplingTable;
      gSBSamplingTable = nullptr;
    }
  }
}

void G4SeltzerBergerModel::Initialise(const G4ParticleDefinition* p,
                                      const G4DataVector& cuts)
{
  if (p) {
    SetParticle(p);
  }
  fIsUseSamplingTables = G4EmParameters::Instance()->EnableSamplingTable();
  // Access to elements
  if (IsMaster()) {

    auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
    for(G4int j=0; j<numOfCouples; ++j) {
      auto mat = theCoupleTable->GetMaterialCutsCouple(j)->GetMaterial();
      auto elmVec = mat->GetElementVector();
      for (auto & elm : *elmVec) {
	G4int Z = std::max(1,std::min(elm->GetZasInt(), gMaxZet-1));
	// load SB-DCS data for this atomic number if it has not been loaded yet
        if (gSBDCSData[Z] == nullptr) ReadData(Z);
      }
    }
    // elem.selectr. only for master: base class init-local will set for workers
    if (LowEnergyLimit() < HighEnergyLimit()) {
      InitialiseElementSelectors(p,cuts);
    }
    // init sampling tables if it was requested
    if (fIsUseSamplingTables) {
      if (!gSBSamplingTable) {
        gSBSamplingTable = new G4SBBremTable();
      }
      gSBSamplingTable->Initialize(std::max(fLowestKinEnergy,LowEnergyLimit()),
                                   HighEnergyLimit());
    }
  }
  //
  if (!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }
  if (GetTripletModel()) {
    GetTripletModel()->Initialise(p, cuts);
    fIsScatOffElectron = true;
  }
}

const G4String& G4SeltzerBergerModel::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4LEDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/brem_SB/br";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4SeltzerBergerModel::FindDirectoryPath()","em0006",
                  FatalException,
                  "Environment variable G4LEDATA not defined");
    }
  }
  return gDataDirectory;
}

void G4SeltzerBergerModel::ReadData(G4int Z) {
  // return if it has been already loaded
  if (gSBDCSData[Z] != nullptr) return;

  G4AutoLock l(&theSBMutex);
  if (gSBDCSData[Z] == nullptr) {
    std::ostringstream ost;
    ost << FindDirectoryPath() << Z;
    std::ifstream fin(ost.str().c_str());
    if (!fin.is_open()) {
      G4ExceptionDescription ed;
      ed << "Bremsstrahlung data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4SeltzerBergerModel::ReadData()","em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.23 or later.");
      return;
    }
    //G4cout << "G4SeltzerBergerModel read from <" << ost.str().c_str()
    //         << ">" << G4endl;
    auto v = new G4Physics2DVector();
    if (v->Retrieve(fin)) {
      v->SetBicubicInterpolation(fIsUseBicubicInterpolation);
      static const G4double emaxlog = 4*G4Log(10.);
      gYLimitData[Z] = v->Value(0.97, emaxlog, fIndx, fIndy);
      gSBDCSData[Z]  = v;
    } else {
      G4ExceptionDescription ed;
      ed << "Bremsstrahlung data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4SeltzerBergerModel::ReadData()","em0005",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.23 or later.");
      delete v;
    }
  }
  l.unlock();
}

G4double G4SeltzerBergerModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
{
  G4double dxsec = 0.0;
  if (gammaEnergy < 0.0 || fPrimaryKinEnergy <= 0.0) {
    return dxsec;
  }
  // reduced photon energy
  const G4double x = gammaEnergy/fPrimaryKinEnergy;
  // l-kinetic energy of the e-/e+
  const G4double y = G4Log(fPrimaryKinEnergy/CLHEP::MeV);
  // make sure that the Z-related SB-DCS are loaded
  // NOTE: fCurrentIZ should have been set before.
  fCurrentIZ = std::max(std::min(fCurrentIZ, gMaxZet-1), 1);
  if (nullptr == gSBDCSData[fCurrentIZ]) {
    ReadData(fCurrentIZ);
  }
  // NOTE: SetupForMaterial should have been called before!
  const G4double pt2   = fPrimaryKinEnergy*(fPrimaryKinEnergy+2.*kMC2);
  const G4double invb2 = fPrimaryTotalEnergy*fPrimaryTotalEnergy/pt2;
  G4double val = gSBDCSData[fCurrentIZ]->Value(x,y,fIndx,fIndy);
  dxsec = val*invb2*CLHEP::millibarn/gBremFactor;
  // e+ correction
  if (!fIsElectron) {
    const G4double invbeta1 = std::sqrt(invb2);
    const G4double e2       = fPrimaryKinEnergy-gammaEnergy;
    if (e2 > 0.0) {
      const G4double invbeta2 = (e2+kMC2)/std::sqrt(e2*(e2+2.0*kMC2));
      const G4double dum0     = kAlpha*fCurrentIZ*(invbeta1-invbeta2);
      if (dum0 < gExpNumLimit) {
        dxsec = 0.0;
      } else {
        dxsec *= G4Exp(dum0);
      }
    } else {
      dxsec = 0.0;
    }
  }
  return dxsec;
}

void
G4SeltzerBergerModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                         const G4MaterialCutsCouple* couple,
                                         const G4DynamicParticle* dp,
                                         G4double cutEnergy,
                                         G4double maxEnergy)
{
  const G4double kinEnergy    = dp->GetKineticEnergy();
  const G4double logKinEnergy = dp->GetLogKineticEnergy();
  const G4double tmin = std::min(cutEnergy, kinEnergy);
  const G4double tmax = std::min(maxEnergy, kinEnergy);
  if (tmin >= tmax) {
    return;
  }
  // set local variables and select target element
  SetupForMaterial(fPrimaryParticle, couple->GetMaterial(), kinEnergy);
  const G4Element* elm = SelectTargetAtom(couple, fPrimaryParticle, kinEnergy,
                                          logKinEnergy, tmin, tmax);
  fCurrentIZ = std::max(std::min(elm->GetZasInt(),gMaxZet-1), 1);
  //
  const G4double totMomentum = std::sqrt(kinEnergy*(fPrimaryTotalEnergy+kMC2));
  /*
  G4cout << "G4SeltzerBergerModel::SampleSecondaries E(MeV)= "
         << kinEnergy/MeV
         << " Z= " << fCurrentIZ << " cut(MeV)= " << tmin/MeV
         << " emax(MeV)= " << tmax/MeV << " corr= " << fDensityCorr << G4endl;
  */
  // sample emitted photon energy either by rejection or from samplign tables
  const G4double gammaEnergy = !fIsUseSamplingTables
        ? SampleEnergyTransfer(kinEnergy, logKinEnergy, tmin, tmax)
        : gSBSamplingTable->SampleEnergyTransfer(kinEnergy, logKinEnergy, tmin, 
                     fDensityCorr, fCurrentIZ, couple->GetIndex(), fIsElectron);
  // should never happen under normal conditions but protect it
  if (gammaEnergy <= 0.) {
    return;
  }
  //
  // angles of the emitted gamma. ( Z - axis along the parent particle) use
  // general interface
  G4ThreeVector gamDir = GetAngularDistribution()->SampleDirection(dp,
            fPrimaryTotalEnergy-gammaEnergy, fCurrentIZ, couple->GetMaterial());
  // create G4DynamicParticle object for the emitted Gamma
  auto gamma = new G4DynamicParticle(fGammaParticle, gamDir, gammaEnergy);
  vdp->push_back(gamma);
  //
  // compute post-interaction kinematics of the primary e-/e+
  G4ThreeVector dir = 
    (totMomentum*dp->GetMomentumDirection()-gammaEnergy*gamDir).unit();
  const G4double finalE = kinEnergy - gammaEnergy;
  /*
  G4cout << "### G4SBModel: v= "
         << " Eg(MeV)= " << gammaEnergy
         << " Ee(MeV)= " << kineticEnergy
         << " DirE " << direction << " DirG " << gammaDirection
         << G4endl;
  */
  // if secondary gamma energy is higher than threshold(very high by default)
  // then stop tracking the primary particle and create new secondary e-/e+
  // instead of the primary
  if (gammaEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    auto el = new G4DynamicParticle(
             const_cast<G4ParticleDefinition*>(fPrimaryParticle), dir, finalE);
    vdp->push_back(el);
  } else { // continue tracking the primary e-/e+ otherwise
    fParticleChange->SetProposedMomentumDirection(dir);
    fParticleChange->SetProposedKineticEnergy(finalE);
  }
}

// sample emitted photon energy by usign rejection
G4double
G4SeltzerBergerModel::SampleEnergyTransfer(const G4double kinEnergy,
                                           const G4double logKinEnergy,
                                           const G4double tmin,
                                           const G4double tmax)
{
  // min max of the transformed variable: x(k) = ln(k^2+k_p^2) that is in
  // [ln(k_c^2+k_p^2), ln(E_k^2+k_p^2)]
  const G4double xmin   = G4Log(tmin*tmin+fDensityCorr);
  const G4double xrange = G4Log(tmax*tmax+fDensityCorr)-xmin;
  const G4double y      = logKinEnergy;
  // majoranta
  const G4double x0 = tmin/kinEnergy;
  G4double vmax;
  if (nullptr == gSBDCSData[fCurrentIZ]) {
    ReadData(fCurrentIZ);
  }
  vmax = gSBDCSData[fCurrentIZ]->Value(x0, y, fIndx, fIndy)*1.02;
  //
  static const G4double kEPeakLim = 300.*CLHEP::MeV;
  static const G4double kELowLim  =  20.*CLHEP::keV;
  // majoranta corrected for e-
  if (fIsElectron && x0 < 0.97 && 
      ((kinEnergy>kEPeakLim) || (kinEnergy<kELowLim))) {
    G4double ylim = std::min(gYLimitData[fCurrentIZ],
                         1.1*gSBDCSData[fCurrentIZ]->Value(0.97,y,fIndx,fIndy));
    vmax = std::max(vmax, ylim);
  }
  if (x0 < 0.05) {
    vmax *= 1.2;
  }
  //G4cout<<"y= "<<y<<" xmin= "<<xmin<<" xmax= "<<xmax
  //<<" vmax= "<<vmax<<G4endl;
  static const G4int kNCountMax = 100;
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  G4double rndm[2];
  G4double gammaEnergy, v;
  for (G4int nn = 0; nn < kNCountMax; ++nn) {
    rndmEngine->flatArray(2, rndm);
    gammaEnergy = 
      std::sqrt(std::max(G4Exp(xmin + rndm[0]*xrange)-fDensityCorr,0.));
    v = gSBDCSData[fCurrentIZ]->Value(gammaEnergy/kinEnergy, y, fIndx, fIndy);
    // e+ correction
    if (!fIsElectron) {
      const G4double       e1 = kinEnergy - tmin;
      const G4double invbeta1 = (e1+kMC2)/std::sqrt(e1*(e1+2.*kMC2));
      const G4double       e2 = kinEnergy-gammaEnergy;
      const G4double invbeta2 = (e2+kMC2)/std::sqrt(e2*(e2+2.*kMC2));
      const G4double     dum0 = kAlpha*fCurrentIZ*(invbeta1-invbeta2);
      if (dum0 < gExpNumLimit) {
        v = 0.0;
      } else {
        v *= G4Exp(dum0);
      }
    }
    if (v > 1.05*vmax && fNumWarnings < 11) {
      ++fNumWarnings;
      G4ExceptionDescription ed;
      ed << "### G4SeltzerBergerModel Warning: Majoranta exceeded! "
         << v << " > " << vmax << " by " << v/vmax
         << " Niter= " << nn
         << " Egamma(MeV)= " << gammaEnergy
         << " Ee(MeV)= " << kinEnergy
         << " Z= " << fCurrentIZ << "  " << fPrimaryParticle->GetParticleName();
      //
      if (10 == fNumWarnings) {
        ed << "\n ### G4SeltzerBergerModel Warnings stopped";
      }
      G4Exception("G4SeltzerBergerModel::SampleScattering","em0044",
                  JustWarning, ed,"");
    }
    if (v >= vmax*rndm[1]) {
      break;
    }
  }
  return gammaEnergy;
}

void G4SeltzerBergerModel::SetupForMaterial(const G4ParticleDefinition*,
                                            const G4Material* mat,
                                            G4double kineticEnergy)
{
  fDensityFactor      = gMigdalConstant*mat->GetElectronDensity();
  // calculate threshold for density effect: gamma*k_p = sqrt(fDensityCorr)
  fPrimaryKinEnergy   = kineticEnergy;
  fPrimaryTotalEnergy = kineticEnergy+CLHEP::electron_mass_c2;
  fDensityCorr        = fDensityFactor*fPrimaryTotalEnergy*fPrimaryTotalEnergy;
  fIsLPMActive        = LPMFlag();
}

