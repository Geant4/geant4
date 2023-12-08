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
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4Physics2DVector.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4AutoLock.hh"

#include "G4ios.hh"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <thread>

G4double G4SeltzerBergerModel::gYLimitData[] = { 0.0 };
G4Physics2DVector* G4SeltzerBergerModel::gSBDCSData[] = { nullptr };
G4SBBremTable* G4SeltzerBergerModel::gSBSamplingTable = nullptr;

// constant DCS factor: 16\alpha r_0^2/3
const G4double G4SeltzerBergerModel::gBremFactor
  = 16. * CLHEP::fine_structure_const * CLHEP::classic_electr_radius
    * CLHEP::classic_electr_radius/3.;

// Migdal's constant: 4\pi r_0*electron_reduced_compton_wavelength^2
const G4double G4SeltzerBergerModel::gMigdalConstant
  = 4. * CLHEP::pi * CLHEP::classic_electr_radius
    * CLHEP::electron_Compton_length * CLHEP::electron_Compton_length;

static constexpr G4double twoMass = 2* CLHEP::electron_mass_c2;
static constexpr G4double kAlpha = CLHEP::twopi*CLHEP::fine_structure_const;
static std::once_flag applyOnce;

namespace
{
  G4Mutex theSBMutex = G4MUTEX_INITIALIZER;

  // for numerical integration on [0,1]
  const G4double gXGL[8] = {
    1.98550718e-02, 1.01666761e-01, 2.37233795e-01, 4.08282679e-01,
    5.91717321e-01, 7.62766205e-01, 8.98333239e-01, 9.80144928e-01
  };
  const G4double gWGL[8] = {
    5.06142681e-02, 1.11190517e-01, 1.56853323e-01, 1.81341892e-01,
    1.81341892e-01, 1.56853323e-01, 1.11190517e-01, 5.06142681e-02
  };
}

G4SeltzerBergerModel::G4SeltzerBergerModel(const G4ParticleDefinition* p,
                                           const G4String& nam)
  : G4VEmModel(nam), 
    fGammaParticle(G4Gamma::Gamma()),
    fLowestKinEnergy(1.0*CLHEP::keV)
{
  SetLowEnergyLimit(fLowestKinEnergy);
  SetAngularDistribution(new G4ModifiedTsai());
  if (fPrimaryParticle != p) { SetParticle(p); }
}

G4SeltzerBergerModel::~G4SeltzerBergerModel()
{
  // delete SB-DCS data per Z
  if (isInitializer) {
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
  // parameters in each thread
  if (fPrimaryParticle != p) {
    SetParticle(p);
  }
  fIsUseSamplingTables = G4EmParameters::Instance()->EnableSamplingTable();
  fCurrentIZ = 0;

  // initialise static tables for the Seltzer-Berger model
  std::call_once(applyOnce, [this]() { isInitializer = true; });

  if (isInitializer) {
    G4AutoLock l(&theSBMutex);

    // initialisation per element is done only once
    auto elemTable = G4Element::GetElementTable();
    for (auto const & elm : *elemTable) {
      G4int Z = std::max(1,std::min(elm->GetZasInt(), gMaxZet-1));
      // load SB-DCS data for this atomic number if it has not been loaded yet
      if (gSBDCSData[Z] == nullptr) ReadData(Z);
    }

    // init sampling tables if it was requested
    if (fIsUseSamplingTables) {
      if (nullptr == gSBSamplingTable) {
        gSBSamplingTable = new G4SBBremTable();
      }
      gSBSamplingTable->Initialize(std::max(fLowestKinEnergy, LowEnergyLimit()),
                                   HighEnergyLimit());
    }
    l.unlock();
  }
  // element selectors are initialized in the master thread
  if (IsMaster()) {
    InitialiseElementSelectors(p, cuts);
  }
  // initialisation in all threads
  if (nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForLoss(); 
  }
  auto trmodel = GetTripletModel();
  if (nullptr != trmodel) {
    trmodel->Initialise(p, cuts);
    fIsScatOffElectron = true;
  }
}

void G4SeltzerBergerModel::InitialiseLocal(const G4ParticleDefinition*,
                                                G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

void G4SeltzerBergerModel::SetParticle(const G4ParticleDefinition* p)
{
  fPrimaryParticle = p;
  fIsElectron = (p == G4Electron::Electron());
}

void G4SeltzerBergerModel::ReadData(G4int Z) {
  // return if it has been already loaded
  if (gSBDCSData[Z] != nullptr) return;

  if (gSBDCSData[Z] == nullptr) {
    std::ostringstream ost;
    ost << G4EmParameters::Instance()->GetDirLEDATA() << "/brem_SB/br" << Z;
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
      gSBDCSData[Z] = v;
    } else {
      G4ExceptionDescription ed;
      ed << "Bremsstrahlung data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4SeltzerBergerModel::ReadData()","em0005",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.23 or later.");
      delete v;
    }
  }
}

// minimum primary (e-/e+) energy at which discrete interaction is possible
G4double G4SeltzerBergerModel::MinPrimaryEnergy(const G4Material*,
                                                const G4ParticleDefinition*,
                                                G4double cut)
{
  return std::max(fLowestKinEnergy, cut);
}

// Sets kinematical variables like E_kin, E_t and some material dependent
// for characteristic photon energy k_p (more exactly
// k_p^2) for the Ter-Mikaelian suppression effect.
void G4SeltzerBergerModel::SetupForMaterial(const G4ParticleDefinition*,
                                            const G4Material* mat,
	                                    G4double kinEnergy)
{
  fDensityFactor = gMigdalConstant*mat->GetElectronDensity();
  // calculate threshold for density effect: k_p = sqrt(fDensityCorr)
  fPrimaryKinEnergy = kinEnergy;
  fPrimaryTotalEnergy = kinEnergy + CLHEP::electron_mass_c2;
  fDensityCorr = fDensityFactor*fPrimaryTotalEnergy*fPrimaryTotalEnergy;
}

// Computes the restricted dE/dx as the appropriate weight of the individual
// element contributions that are computed by numerically integrating the DCS.
G4double
G4SeltzerBergerModel::ComputeDEDXPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                           G4double kineticEnergy,
                                           G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (nullptr == fPrimaryParticle) {
    SetParticle(p);
  }
  if (kineticEnergy <= fLowestKinEnergy) {
    return dedx;
  }
  // maximum value of the dE/dx integral (the minimum is 0 of course)
  G4double tmax = std::min(cutEnergy, kineticEnergy);
  if (tmax == 0.0) {
    return dedx;
  }
  // sets kinematical and material related variables
  SetupForMaterial(fPrimaryParticle, material, kineticEnergy);
  // get element compositions of the material
  const G4ElementVector* theElemVector = material->GetElementVector();
  const G4double* theAtomNumDensVector = material->GetAtomicNumDensityVector();
  const std::size_t numberOfElements = theElemVector->size();
  // loop over the elements of the material and compute their contributions to
  // the restricted dE/dx by numerical integration of the dependent part of DCS
  for (std::size_t ie = 0; ie < numberOfElements; ++ie) {
    G4VEmModel::SetCurrentElement((*theElemVector)[ie]);
    G4int Z = (*theElemVector)[ie]->GetZasInt();
    fCurrentIZ = std::min(Z, gMaxZet);
    dedx += (Z*Z)*theAtomNumDensVector[ie]*ComputeBremLoss(tmax);
  }
  // apply the constant factor C/Z = 16\alpha r_0^2/3
  dedx *= gBremFactor;
  return std::max(dedx, 0.);
}

// Computes the integral part of the restricted dE/dx contribution from a given
// element (Z) by numerically integrating the k dependent DCS between
// k_min=0 and k_max = tmax = min[gamma-cut, electron-kinetic-energy].
// The numerical integration is done by dividing the integration range into 'n'
// subintervals and an 8 pint GL integral (on [0,1]) is performed on each sub-
// inteval by tranforming k to alpha=k/E_t (E_t is the total energy of the e-)
// and each sub-interavl is transformed to [0,1]. So the integrastion is done
// in xi(alpha) = xi(k) = [k/E_t-alpha_i]/delta where alpha_i=(i-1)*delta for
// the i = 1,2,..,n-th sub-interval so xi(k) in [0,1] on each sub-intevals.
// This transformation from 'k' to 'xi(k)' results in a multiplicative factor
// of E_t*delta at each step.
// The restricted dE/dx = N int_{0}^{k_max} k*ds/dk dk. In this case not
// the ds/dk(Z,k) but ds/dk(Z,k)*[F*k/C] is computed since:
// (i)    what we need here is ds/dk*k and not k so this multiplication is done
// (ii)   the Ter-Mikaelian suppression i.e. F related factor is done here
// (iii)  the constant factor C (includes Z^2 as well)is accounted in the caller
G4double G4SeltzerBergerModel::ComputeBremLoss(G4double tmax)
{
  // number of intervals and integration step
  const G4double alphaMax = tmax/fPrimaryTotalEnergy;
  const G4int nSub = (G4int)(20*alphaMax)+3;
  const G4double delta = alphaMax/((G4double)nSub);
  // set minimum value of the first sub-inteval
  G4double alpha_i = 0.0;
  G4double dedxInteg = 0.0;
  for (G4int l = 0; l < nSub; ++l) {
    for (G4int igl = 0; igl < 8; ++igl) {
      // compute the emitted photon energy k
      const G4double k   = (alpha_i+gXGL[igl]*delta)*fPrimaryTotalEnergy;
      // compute the DCS value at k (without the constant, the 1/k, 1/F factors)
      const G4double dcs = ComputeDXSectionPerAtom(k);
      // account Ter-Mikaelian suppression: times 1/F with F = 1+(k_p/k)^2
      dedxInteg += gWGL[igl]*dcs/(1.0+fDensityCorr/(k*k));
    }
    // update sub-interval minimum value
    alpha_i += delta;
  }
  // apply corrections due to variable transformation i.e. E_t*delta
  dedxInteg *= delta*fPrimaryTotalEnergy;
  return std::max(dedxInteg,0.);
}

// Computes restrected atomic cross section by numerically integrating the
// DCS between the proper kinematical limits accounting the gamma production cut
G4double
G4SeltzerBergerModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* p,
						 G4double kineticEnergy,
						 G4double Z,
						 G4double,
						 G4double cut,
						 G4double maxEnergy)
{
  G4double crossSection = 0.0;
  if (nullptr == fPrimaryParticle) {
    SetParticle(p);
  }
  if (kineticEnergy <= fLowestKinEnergy) {
    return crossSection;
  }
  // min/max kinetic energy limits of the DCS integration:
  const G4double tmin = std::min(cut, kineticEnergy);
  const G4double tmax = std::min(maxEnergy, kineticEnergy);
  // zero restricted x-section if e- kinetic energy is below gamma cut
  if (tmin >= tmax) {
    return crossSection;
  }
  fCurrentIZ = std::min(G4lrint(Z), gMaxZet);
  // integrate numerically (dependent part of) the DCS between the kin. limits:
  // a. integrate between tmin and kineticEnergy of the e-
  crossSection = ComputeXSectionPerAtom(tmin);
  // allow partial integration: only if maxEnergy < kineticEnergy
  // b. integrate between tmax and kineticEnergy (tmax=maxEnergy in this case)
  // (so the result in this case is the integral of DCS between tmin and
  // maxEnergy)
  if (tmax < kineticEnergy) {
    crossSection -= ComputeXSectionPerAtom(tmax);
  }
  // multiply with the constant factors: 16\alpha r_0^2/3 Z^2
  crossSection *= Z*Z*gBremFactor;
  return std::max(crossSection, 0.);
}

// Numerical integral of the (k dependent part of) DCS between k_min=tmin and
// k_max = E_k (where E_k is the kinetic energy of the e- and tmin is the
// minimum of energy of the  emitted photon). The integration is done in the
// transformed alpha(k) = ln(k/E_t) variable (with E_t being the total energy of
// the primary e-). The integration range is divided into n sub-intervals with
// delta = [ln(k_min/E_t)-ln(k_max/E_t)]/n width each. An 8 point GL integral
// on [0,1] is applied on each sub-inteval so alpha is transformed to
// xi(alpha) = xi(k) = [ln(k/E_t)-alpha_i]/delta where alpha_i = ln(k_min/E_t) +
// (i-1)*delta for the i = 1,2,..,n-th sub-interval and xi(k) in [0,1] on each
// sub-intevals. From the transformed xi, k(xi) = E_t exp[xi*delta+alpha_i].
// Since the integration is done in variable xi instead of k this
// transformation results in a multiplicative factor of k*delta at each step.
// However, DCS differential in k is ~1/k so the multiplicative factor is simple
// becomes delta and the 1/k factor is dropped from the DCS computation.
// Ter-Mikaelian suppression is always accounted
G4double G4SeltzerBergerModel::ComputeXSectionPerAtom(G4double tmin)
{
  G4double xSection = 0.0;
  const G4double alphaMin = G4Log(tmin/fPrimaryTotalEnergy);
  const G4double alphaMax = G4Log(fPrimaryKinEnergy/fPrimaryTotalEnergy);
  const G4int    nSub = (G4int)(0.45*(alphaMax-alphaMin))+4;
  const G4double delta = (alphaMax-alphaMin)/((G4double)nSub);
  // set minimum value of the first sub-inteval
  G4double alpha_i = alphaMin;
  for (G4int l = 0; l < nSub; ++l) {
    for (G4int igl = 0; igl < 8; ++igl) {
      // compute the emitted photon energy k
      const G4double k = G4Exp(alpha_i+gXGL[igl]*delta)*fPrimaryTotalEnergy;
      // compute the DCS value at k (without the constant, the 1/k, 1/F factors)
      const G4double dcs = ComputeDXSectionPerAtom(k);
      // account Ter-Mikaelian suppression: times 1/F with F = 1+(k_p/k)^2
      xSection += gWGL[igl]*dcs/(1.0+fDensityCorr/(k*k));
    }
    // update sub-interval minimum value
    alpha_i += delta;
  }
  // apply corrections due to variable transformation
  xSection *= delta;
  // final check
  return std::max(xSection, 0.);
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
    G4AutoLock l(&theSBMutex);
    ReadData(fCurrentIZ);
    l.unlock();
  }
  // NOTE: SetupForMaterial should have been called before!
  const G4double pt2 = fPrimaryKinEnergy*(fPrimaryKinEnergy + twoMass);
  const G4double invb2 = fPrimaryTotalEnergy*fPrimaryTotalEnergy/pt2;
  G4double val = gSBDCSData[fCurrentIZ]->Value(x,y,fIndx,fIndy);
  dxsec = val*invb2*CLHEP::millibarn/gBremFactor;
  // e+ correction
  if (!fIsElectron) {
    const G4double invbeta1 = std::sqrt(invb2);
    const G4double e2 = fPrimaryKinEnergy - gammaEnergy;
    if (e2 > 0.0) {
      const G4double invbeta2 =
	(e2 + CLHEP::electron_mass_c2)/std::sqrt(e2*(e2 + twoMass));
      const G4double dum0 = kAlpha*fCurrentIZ*(invbeta1-invbeta2);
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
  const G4double kinEnergy = dp->GetKineticEnergy();
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
  fCurrentIZ = std::max(std::min(elm->GetZasInt(), gMaxZet-1), 1);
  //
  const G4double totMomentum = std::sqrt(kinEnergy*(kinEnergy + twoMass));
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
    (totMomentum*dp->GetMomentumDirection() - gammaEnergy*gamDir).unit();
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
      const G4double e1 = kinEnergy - tmin;
      const G4double invbeta1 =
	(e1 + CLHEP::electron_mass_c2)/std::sqrt(e1*(e1 + twoMass));
      const G4double       e2 = kinEnergy-gammaEnergy;
      const G4double invbeta2 =
	(e2 + CLHEP::electron_mass_c2)/std::sqrt(e2*(e2 + twoMass));
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
