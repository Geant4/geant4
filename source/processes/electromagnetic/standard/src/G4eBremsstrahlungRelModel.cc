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
// File name:     G4eBremsstrahlungRelModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 12.08.2008
//
// Modifications:
//
// 13.11.08    add SetLPMflag and SetLPMconstant methods
// 13.11.08    change default LPMconstant value
// 13.10.10    add angular distributon interface (VI)
// 31.05.16    change LPMconstant such that it gives suppression variable 's'
//             that consistent to Migdal's one; fix a small bug in 'logTS1'
//             computation; better agreement with exp.(M.Novak)
// 15.07.18    improved LPM suppression function approximation (no artificial
//             steps), code cleanup and optimizations,more implementation and
//             model related comments, consistent variable naming (M.Novak)
//
// Main References:
//  Y.-S.Tsai, Rev. Mod. Phys. 46 (1974) 815; Rev. Mod. Phys. 49 (1977) 421.
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media,
//  Wiley, 1972.
//
// -------------------------------------------------------------------
//

#include "G4eBremsstrahlungRelModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ModifiedTsai.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

const G4int G4eBremsstrahlungRelModel::gMaxZet = 120;

// constant DCS factor: 16\alpha r_0^2/3
const G4double G4eBremsstrahlungRelModel::gBremFactor
  = 16. * CLHEP::fine_structure_const * CLHEP::classic_electr_radius
    * CLHEP::classic_electr_radius/3.;

// Migdal's constant: 4\pi r_0*electron_reduced_compton_wavelength^2
const G4double G4eBremsstrahlungRelModel::gMigdalConstant
  = 4. * CLHEP::pi * CLHEP::classic_electr_radius
    * CLHEP::electron_Compton_length * CLHEP::electron_Compton_length;

// LPM constant: \alpha(mc^2)^2/(4\pi*\hbar c)
const G4double G4eBremsstrahlungRelModel::gLPMconstant
  = CLHEP::fine_structure_const * CLHEP::electron_mass_c2
    * CLHEP::electron_mass_c2 / (4. * CLHEP::pi * CLHEP::hbarc);

// abscissas and weights of an 8 point Gauss-Legendre quadrature
// for numerical integration on [0,1]
const G4double G4eBremsstrahlungRelModel::gXGL[] = {
  1.98550718e-02, 1.01666761e-01, 2.37233795e-01, 4.08282679e-01,
  5.91717321e-01, 7.62766205e-01, 8.98333239e-01, 9.80144928e-01
};
const G4double G4eBremsstrahlungRelModel::gWGL[] = {
  5.06142681e-02, 1.11190517e-01, 1.56853323e-01, 1.81341892e-01,
  1.81341892e-01, 1.56853323e-01, 1.11190517e-01, 5.06142681e-02
};

// elastic and inelatic radiation logarithms for light elements (where the
// Thomas-Fermi model doesn't work): computed by using Dirac-Fock model of atom.
const G4double G4eBremsstrahlungRelModel::gFelLowZet  [] = {
  0.0, 5.3104, 4.7935, 4.7402, 4.7112, 4.6694, 4.6134, 4.5520
};
const G4double G4eBremsstrahlungRelModel::gFinelLowZet[] = {
  0.0, 5.9173, 5.6125, 5.5377, 5.4728, 5.4174, 5.3688, 5.3236
};

// LPM supression functions evaluated at initialisation time
G4eBremsstrahlungRelModel::LPMFuncs  G4eBremsstrahlungRelModel::gLPMFuncs;

// special data structure per element i.e. per Z
std::vector<G4eBremsstrahlungRelModel::ElementData*> G4eBremsstrahlungRelModel::gElementData;

G4eBremsstrahlungRelModel::G4eBremsstrahlungRelModel(const G4ParticleDefinition* p,
                                                     const G4String& nam)
: G4VEmModel(nam)
{
  fGammaParticle       = G4Gamma::Gamma();
  //
  fLowestKinEnergy     = 1.0*MeV;
  SetLowEnergyLimit(fLowestKinEnergy);
  //
  fLPMEnergyThreshold  = 1.e+39;
  fLPMEnergy           = 0.;

  SetLPMFlag(true);
  //
  SetAngularDistribution(new G4ModifiedTsai());
  //
  if (nullptr != p) {
    SetParticle(p);
  }
}

G4eBremsstrahlungRelModel::~G4eBremsstrahlungRelModel()
{
  if (IsMaster()) {
    // clear ElementData container
    for (std::size_t iz = 0; iz < gElementData.size(); ++iz) {
      if (nullptr != gElementData[iz]) {
        delete gElementData[iz];
      }
    }
    gElementData.clear();
    // clear LPMFunctions (if any)
    if (LPMFlag()) {
      gLPMFuncs.fLPMFuncG.clear();
      gLPMFuncs.fLPMFuncPhi.clear();
      gLPMFuncs.fIsInitialized = false;
    }
  }
}

void G4eBremsstrahlungRelModel::Initialise(const G4ParticleDefinition* p,
                                           const G4DataVector& cuts)
{
  if (nullptr != p) {
    SetParticle(p);
  }
  fCurrentIZ = 0;
  // init element data and precompute LPM functions (only if lpmflag is true)
  if (IsMaster()) {
    InitialiseElementData();
    if (LPMFlag()) { InitLPMFunctions(); }
    if (LowEnergyLimit() < HighEnergyLimit()) {
      InitialiseElementSelectors(p, cuts);
    }
  }
  if (nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForLoss(); 
  }
  if (GetTripletModel()) {
    GetTripletModel()->Initialise(p, cuts);
    fIsScatOffElectron = true;
  }
}

void G4eBremsstrahlungRelModel::InitialiseLocal(const G4ParticleDefinition*,
                                                G4VEmModel* masterModel)
{
  if (LowEnergyLimit() < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

void G4eBremsstrahlungRelModel::SetParticle(const G4ParticleDefinition* p)
{
  fPrimaryParticle     = p;
  fPrimaryParticleMass = p->GetPDGMass();
  fIsElectron          = (p==G4Electron::Electron());
}

// Sets kinematical variables like E_kin, E_t and some material dependent
// variables like LPM energy and characteristic photon energy k_p (more exactly
// k_p^2) for the Ter-Mikaelian suppression effect.
void G4eBremsstrahlungRelModel::SetupForMaterial(const G4ParticleDefinition*,
                                                 const G4Material* mat,
	                                               G4double kineticEnergy)
{
  fDensityFactor = gMigdalConstant*mat->GetElectronDensity();
  fLPMEnergy     = gLPMconstant*mat->GetRadlen();
  // threshold for LPM effect (i.e. below which LPM hidden by density effect)
  if (LPMFlag()) {
    fLPMEnergyThreshold = std::sqrt(fDensityFactor)*fLPMEnergy;
  } else {
    fLPMEnergyThreshold = 1.e+39;   // i.e. do not use LPM effect
  }
  // calculate threshold for density effect: k_p = sqrt(fDensityCorr)
  fPrimaryKinEnergy   = kineticEnergy;
  fPrimaryTotalEnergy = kineticEnergy+fPrimaryParticleMass;
  fDensityCorr        = fDensityFactor*fPrimaryTotalEnergy*fPrimaryTotalEnergy;
  // set activation flag for LPM effects in the DCS
  fIsLPMActive        = (fPrimaryTotalEnergy>fLPMEnergyThreshold);
}

// minimum primary (e-/e+) energy at which discrete interaction is possible
G4double G4eBremsstrahlungRelModel::MinPrimaryEnergy(const G4Material*,
                                                     const G4ParticleDefinition*,
                                                     G4double cut)
{
  return std::max(fLowestKinEnergy, cut);
}

// Computes the restricted dE/dx as the appropriate weight of the individual
// element contributions that are computed by numerically integrating the DCS.
G4double
G4eBremsstrahlungRelModel::ComputeDEDXPerVolume(const G4Material* material,
                                                const G4ParticleDefinition* p,
                                                G4double kineticEnergy,
                                                G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (nullptr == fPrimaryParticle) {
    SetParticle(p);
  }
  if (kineticEnergy < LowEnergyLimit()) {
    return dedx;
  }
  // maximum value of the dE/dx integral (the minimum is 0 of course)
  G4double tmax = std::min(cutEnergy, kineticEnergy);
  if (tmax == 0.0) {
    return dedx;
  }
  // sets kinematical and material related variables
  SetupForMaterial(fPrimaryParticle, material,kineticEnergy);
  // get element compositions of the material
  const G4ElementVector* theElemVector = material->GetElementVector();
  const G4double* theAtomNumDensVector = material->GetAtomicNumDensityVector();
  const std::size_t numberOfElements = theElemVector->size();
  // loop over the elements of the material and compute their contributions to
  // the restricted dE/dx by numerical integration of the dependent part of DCS
  for (std::size_t ie = 0; ie < numberOfElements; ++ie) {
    G4VEmModel::SetCurrentElement((*theElemVector)[ie]);
    G4int zet = (*theElemVector)[ie]->GetZasInt();
    fCurrentIZ = std::min(zet, gMaxZet);
    dedx              += (zet*zet)*theAtomNumDensVector[ie]*ComputeBremLoss(tmax);
  }
  // apply the constant factor C/Z = 16\alpha r_0^2/3
  dedx *= gBremFactor;
  return std::max(dedx,0.);
}

// Computes the integral part of the restricted dE/dx contribution from a given
// element (Z) by numerically integrating the k dependent part of the DCS between
// k_min=0 and k_max = tmax = min[gamma-cut, electron-kinetic-eenrgy].
// The numerical integration is done by dividing the integration range into 'n'
// subintervals and an 8 pint GL integral (on [0,1]) is performed on each sub-
// inteval by tranforming k to alpha=k/E_t (E_t is the total energy of the e-)
// and each sub-interavl is transformed to [0,1]. So the integrastion is done
// in xi(alpha) = xi(k) = [k/E_t-alpha_i]/delta where alpha_i=(i-1)*delta for
// the i = 1,2,..,n-th sub-interval so xi(k) in [0,1] on each sub-intevals.
// This transformation from 'k' to 'xi(k)' results in a multiplicative factor
// of E_t*delta at each step.
// The restricted dE/dx = N int_{0}^{k_max} k*ds/dk dk. There are 2 DCS model
// one with LPM and one without LPM effects (see them below). In both case not
// the ds/dk(Z,k) but ds/dk(Z,k)*[F*k/C] is computed since:
// (i)    what we need here is ds/dk*k and not k so this multiplication is done
// (ii)   the Ter-Mikaelian suppression i.e. F related factor is done here
// (iii)  the constant factor C (includes Z^2 as well)is accounted in the caller
G4double G4eBremsstrahlungRelModel::ComputeBremLoss(G4double tmax)
{
  // number of intervals and integration step
  const G4double alphaMax = tmax/fPrimaryTotalEnergy;
  const G4int        nSub = (G4int)(20*alphaMax)+3;
  const G4double    delta = alphaMax/((G4double)nSub);
  // set minimum value of the first sub-inteval
  G4double alpha_i        = 0.0;
  G4double dedxInteg      = 0.0;
  for (G4int l = 0; l < nSub; ++l) {
    for (G4int igl = 0; igl < 8; ++igl) {
      // compute the emitted photon energy k
      const G4double k   = (alpha_i+gXGL[igl]*delta)*fPrimaryTotalEnergy;
      // compute the DCS value at k (without the constant, the 1/k, 1/F factors)
      const G4double dcs = fIsLPMActive
                          ? ComputeRelDXSectionPerAtom(k)  // DCS WITHOUT LPM
                          : ComputeDXSectionPerAtom(k);    // DCS WITH    LPM
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
G4double G4eBremsstrahlungRelModel::ComputeCrossSectionPerAtom(
                                                  const G4ParticleDefinition* p,
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
  if (kineticEnergy < LowEnergyLimit()) {
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
// NOTE:
//   - LPM suppression is accounted above threshold e- energy (corresponidng
//     flag is set in SetUpForMaterial() => 2 DCS with/without LPM
//   - Ter-Mikaelian suppression is always accounted
G4double G4eBremsstrahlungRelModel::ComputeXSectionPerAtom(G4double tmin)
{
  G4double xSection = 0.0;
  const G4double alphaMin = G4Log(tmin/fPrimaryTotalEnergy);
  const G4double alphaMax = G4Log(fPrimaryKinEnergy/fPrimaryTotalEnergy);
  const G4int    nSub     = (G4int)(0.45*(alphaMax-alphaMin))+4;
  const G4double delta    = (alphaMax-alphaMin)/((G4double)nSub);
  // set minimum value of the first sub-inteval
  G4double alpha_i        = alphaMin;
  for (G4int l = 0; l < nSub; ++l) {
    for (G4int igl = 0; igl < 8; ++igl) {
      // compute the emitted photon energy k
      const G4double k   = G4Exp(alpha_i+gXGL[igl]*delta)*fPrimaryTotalEnergy;
      // compute the DCS value at k (without the constant, the 1/k, 1/F factors)
      const G4double dcs = fIsLPMActive
                          ? ComputeRelDXSectionPerAtom(k) // DCS WITHOUT LPM
                          : ComputeDXSectionPerAtom(k);   // DCS WITH    LPM
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

// DCS WITH LPM EFFECT: complete screening aprx. and includes LPM suppression
// ds/dk(Z,k) = C/[F*k]*{ Xi(s*F)*[y^2*G/4 +(1-y+y^2/3)Phi]*[L_el-f_c+L_inel/Z]
//                        +(1-y)*[1+1/Z]/12}  with C = 16\alpha r_0^2/3 Z^2 and
// Xi(s),G(s), Phi(s) are LPM suppression functions:
//
// LPM SUPPRESSION: The 's' is the suppression variable and F = F(k,k_p) =
// 1+(k_p/k)^2 with k_p = hbar*w_p*E/(m*c^2) is a material (e- density)
// dependent constant. F accounts the Ter-Mikaelian suppression with a smooth
// transition in the emitted photon energy. Also, the LPM suppression functions
// goes to 0 when s goes to 0 and goes to 1 when s is increasing (=1 at s=~2)
// So evaluating the LPM suppression functions at 'sF' instead of 's' ensures a
// smooth transition depending on the emitted photon energy 'k': LPM effect is
// smoothly turned off i.e. Xi(sF)=G(sF)=Phi(sF)=1 when k << k_p because F >> 1
// and sF ~ s when k >> k_p since F ~ 1 in that case.
// HERE, ds/dk(Z,k)*[F*k/C] is computed since:
//  (i)   DCS ~ 1/k factor will disappear due to the variable transformation
//        v(k)=ln(k/E_t) -> dk/dv=E_t*e^v=k -> ds/dv= ds/dk*dk/dv=ds/dk*k so it
//        would cnacell out the 1/k factor => 1/k don't included here
//  (ii)  the constant factor C and Z don't depend on 'k' => not included here
//  (iii) the 1/F(k) factor is accounted in the callers: explicitly (cross sec-
//        tion computation) or implicitly through further variable transformaton
//        (in the final state sampling algorithm)
// COMPLETE SCREENING: see more at the DCS without LPM effect below.
G4double
G4eBremsstrahlungRelModel::ComputeRelDXSectionPerAtom(G4double gammaEnergy)
{
  G4double dxsec = 0.0;
  if (gammaEnergy < 0.) {
    return dxsec;
  }
  const G4double y     = gammaEnergy/fPrimaryTotalEnergy;
  const G4double onemy = 1.-y;
  const G4double dum0  = 0.25*y*y;
  // evaluate LPM functions (combined with the Ter-Mikaelian effect)
  G4double funcGS, funcPhiS, funcXiS;
  ComputeLPMfunctions(funcXiS, funcGS, funcPhiS, gammaEnergy);
  const ElementData* elDat = gElementData[fCurrentIZ];
  const G4double term1     = funcXiS*(dum0*funcGS+(onemy+2.0*dum0)*funcPhiS);
  dxsec = term1*elDat->fZFactor1+onemy*elDat->fZFactor2;
  //
  if (fIsScatOffElectron) {
    fSumTerm = dxsec;
    fNucTerm = term1*elDat->fZFactor11 + onemy/12.;
  }
  return std::max(dxsec,0.0);
}

// DCS WITHOUT LPM EFFECT: DCS with sceening (Z>5) and Coulomb cor. no LPM
// ds/dk(Z,k)=C/[F*k]*{(1-y+3*y^2/4)*[(0.25*phi1(g)-ln(Z)/3-f_c)+(0.25*psi1(e)
// -2*ln(Z)/3)/Z]+ (1-y)*[(phi1(g)-phi2(g))+(psi1(e)-psi2(e))/Z]/8}
// where f_c(Z) is the Coulomb correction factor and phi1(g),phi2(g) and psi1(e),
// psi2(e) are coherent and incoherent screening functions. In the Thomas-Fermi
// model of the atom, the screening functions will have a form that do not
// depend on Z (not explicitly). These numerical screening functions can be
// approximated as Tsai Eqs. [3.38-3.41] with the variables g=gamma and
// e=epsilon given by Tsai Eqs. [3.30 and 3.31] (see more details at the method
// ComputeScreeningFunctions()). Note, that in case of complete screening i.e.
// g = e = 0 => 0.25*phi1(0)-ln(Z)/3 = ln(184.149/Z^(1/3)) = L_el and
// 0.25*psi1(0)-2*ln(Z)/3=ln(1193.923/Z^(2/3))=L_inel and phi1(0)-phi2(0) =
// psi1(0)-psi2(0) = 2/3 so the DCS in complete screening =>
// COMPLETE SCREENING:
// ds/dk(Z,k)=C/k*{(1-y+3*y^2/4)*[L_el-f_c+L_inel/Z] + (1-y)*[1+1/Z]/12} that is
// used in case of DCS with LPM above (if all the suprression functions are
// absent i.e. their value = 1).
// Since the Thomas-Fermi model of the atom is not accurate at low Z, the DCS in
// complete screening is used here at low Z(<5) with L_el(Z), L_inel(Z) values
// computed by using the Dirac-Fock model of the atom.
// NOTE: that the Ter-Mikaelian suppression is accounted in the DCS through the
// 1/F factor but it is included in the caller and not considered here.
// HERE, ds/dk(Z,k)*[F*k/C] is computed exactly like in the DCS with LPM case.
G4double
G4eBremsstrahlungRelModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
{
  G4double dxsec = 0.0;
  if (gammaEnergy < 0.) {
    return dxsec;
  }
  const G4double y         = gammaEnergy/fPrimaryTotalEnergy;
  const G4double onemy     = 1.-y;
  const G4double dum0      = onemy+0.75*y*y;
  const ElementData* elDat = gElementData[fCurrentIZ];
  // use complete screening and L_el, L_inel from Dirac-Fock model instead of TF
  if (fCurrentIZ < 5 || fIsUseCompleteScreening) {
    dxsec  = dum0*elDat->fZFactor1;
    dxsec += onemy*elDat->fZFactor2;
    if (fIsScatOffElectron) {
      fSumTerm = dxsec;
      fNucTerm = dum0*elDat->fZFactor11+onemy/12.;
    }
  } else {
    // use Tsai's analytical approx. (Tsai Eqs. [3.38-3.41]) to the 'universal'
    // numerical screening functions computed by using the TF model of the atom
    const G4double invZ    = 1./(G4double)fCurrentIZ;
    const G4double Fz      = elDat->fFz;
    const G4double logZ    = elDat->fLogZ;
    const G4double dum1    = y/(fPrimaryTotalEnergy-gammaEnergy);
    const G4double gamma   = dum1*elDat->fGammaFactor;
    const G4double epsilon = dum1*elDat->fEpsilonFactor;
    // evaluate the screening functions
    G4double phi1, phi1m2, psi1, psi1m2;
    ComputeScreeningFunctions(phi1, phi1m2, psi1, psi1m2, gamma, epsilon);
    dxsec  = dum0*((0.25*phi1-Fz) + (0.25*psi1-2.*logZ/3.)*invZ);
    dxsec += 0.125*onemy*(phi1m2 + psi1m2*invZ);
    if (fIsScatOffElectron) {
      fSumTerm = dxsec;
      fNucTerm = dum0*(0.25*phi1-Fz) + 0.125*onemy*phi1m2;
    }
  }
  return std::max(dxsec,0.0);
}

// Coherent and incoherent screening function approximations (see Tsai
// Eqs.[3.38-3.41]). Tsai's analytical approximations to the numerical screening
// functions computed by using the Thomas-Fermi model of atom (Moliere's appro-
// ximation to the numerical TF screening function). In the TF-model, these
// screening functions can be expressed in a 'universal' i.e. Z (directly) inde-
// pendent variable (see Tsai Eqs. Eqs. [3.30 and 3.31]).
void G4eBremsstrahlungRelModel::ComputeScreeningFunctions(G4double& phi1,
                                                          G4double& phi1m2,
                                                          G4double& psi1,
                                                          G4double& psi1m2,
                                                          const G4double gam,
                                                          const G4double eps)
{
  const G4double gam2 = gam*gam;
  phi1   = 16.863-2.0*G4Log(1.0+0.311877*gam2)+2.4*G4Exp(-0.9*gam)
          +1.6*G4Exp(-1.5*gam);
  phi1m2 = 2.0/(3.0+19.5*gam+18.0*gam2);    // phi1-phi2
  const G4double eps2 = eps*eps;
  psi1   = 24.34-2.0*G4Log(1.0+13.111641*eps2)+2.8*G4Exp(-8.0*eps)
          +1.2*G4Exp(-29.2*eps);
  psi1m2 = 2.0/(3.0+120.0*eps+1200.0*eps2); //psi1-psi2
}

void
G4eBremsstrahlungRelModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                             const G4MaterialCutsCouple* couple,
                                             const G4DynamicParticle* dp,
                                             G4double cutEnergy,
                                             G4double maxEnergy)
{
  const G4double kineticEnergy    = dp->GetKineticEnergy();
//  const G4double logKineticEnergy = dp->GetLogKineticEnergy();
  if (kineticEnergy < LowEnergyLimit()) {
    return;
  }
  // min, max kinetic energy limits
  const G4double tmin = std::min(cutEnergy, kineticEnergy);
  const G4double tmax = std::min(maxEnergy, kineticEnergy);
  if (tmin >= tmax) {
    return;
  }
  //
  SetupForMaterial(fPrimaryParticle, couple->GetMaterial(), kineticEnergy);
  const G4Element* elm = SelectTargetAtom(couple,fPrimaryParticle,kineticEnergy,
                                          dp->GetLogKineticEnergy(),tmin,tmax);
  //
  fCurrentIZ = elm->GetZasInt();
  const ElementData* elDat = gElementData[fCurrentIZ];
  const G4double funcMax = elDat->fZFactor1+elDat->fZFactor2;
  // get the random engine
  G4double rndm[2];
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  // min max of the transformed variable: x(k) = ln(k^2+k_p^2) that is in [ln(k_c^2+k_p^2), ln(E_k^2+k_p^2)]
  const G4double xmin   = G4Log(tmin*tmin+fDensityCorr);
  const G4double xrange = G4Log(tmax*tmax+fDensityCorr)-xmin;
  G4double gammaEnergy, funcVal;
  do {
    rndmEngine->flatArray(2, rndm);
    gammaEnergy = std::sqrt(std::max(G4Exp(xmin+rndm[0]*xrange)-fDensityCorr, 0.0));
    funcVal     = fIsLPMActive
                 ? ComputeRelDXSectionPerAtom(gammaEnergy)
                 : ComputeDXSectionPerAtom(gammaEnergy);
    // cross-check of proper function maximum in the rejection
//    if (funcVal > funcMax) {
//      G4cout << "### G4eBremsstrahlungRelModel Warning: Majoranta exceeded! "
//	     << funcVal << " > " << funcMax
//	     << " Egamma(MeV)= " << gammaEnergy
//	     << " Ee(MeV)= " << kineticEnergy
//	     << "  " << GetName()
//	     << G4endl;
//    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (funcVal < funcMax*rndm[1]);
  //
  // scattering off nucleus or off e- by triplet model
  if (fIsScatOffElectron && rndmEngine->flat()*fSumTerm>fNucTerm) {
    GetTripletModel()->SampleSecondaries(vdp, couple, dp, cutEnergy, maxEnergy);
    return;
  }
  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  G4ThreeVector gamDir =
    GetAngularDistribution()->SampleDirection(dp,fPrimaryTotalEnergy-gammaEnergy,
                                              fCurrentIZ, couple->GetMaterial());
  // create G4DynamicParticle object for the Gamma
  auto gamma = new G4DynamicParticle(fGammaParticle, gamDir, gammaEnergy);
  vdp->push_back(gamma);
  // compute post-interaction kinematics of primary e-/e+ based on
  // energy-momentum conservation
  const G4double totMomentum = std::sqrt(kineticEnergy*(
                               fPrimaryTotalEnergy + CLHEP::electron_mass_c2));
  G4ThreeVector dir =
             (totMomentum*dp->GetMomentumDirection()-gammaEnergy*gamDir).unit();
  const G4double finalE   = kineticEnergy-gammaEnergy;
  // if secondary gamma energy is higher than threshold(very high by default)
  // then stop tracking the primary particle and create new secondary e-/e+
  // instead of the primary one
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

void G4eBremsstrahlungRelModel::InitialiseElementData()
{
  const G4int size = (G4int)gElementData.size();
  if (size < gMaxZet+1) {
    gElementData.resize(gMaxZet+1, nullptr);
  }
  // create for all elements that are in the detector
  const G4ElementTable* elemTable = G4Element::GetElementTable();
  std::size_t numElems = (*elemTable).size();
  for (std::size_t ielem=0; ielem<numElems; ++ielem) {
    const G4Element* elem = (*elemTable)[ielem];
    const G4double    zet = elem->GetZ();
    const G4int      izet = std::min(G4lrint(zet),gMaxZet);
    if (!gElementData[izet]) {
      auto elemData  = new ElementData();
      const G4double fc = elem->GetfCoulomb();
      G4double Fel      = 1.;
      G4double Finel    = 1.;
      elemData->fLogZ   = G4Log(zet);
      elemData->fFz     = elemData->fLogZ/3.+fc;
      if (izet < 5) {
        Fel   = gFelLowZet[izet];
        Finel = gFinelLowZet[izet];
      } else {
        Fel   = G4Log(184.15) -    elemData->fLogZ/3.;
        Finel = G4Log(1194)   - 2.*elemData->fLogZ/3.;
      }
      const G4double z23       = std::pow(zet,2./3.);
      const G4double z13       = std::pow(zet,1./3.);
      elemData->fZFactor1      = (Fel-fc)+Finel/zet;
      elemData->fZFactor11     = (Fel-fc); // used only for the triplet
      elemData->fZFactor2      = (1.+1./zet)/12.;
      elemData->fVarS1         = z23/(184.15*184.15);
      elemData->fILVarS1Cond   = 1./(G4Log(std::sqrt(2.0)*elemData->fVarS1));
      elemData->fILVarS1       = 1./G4Log(elemData->fVarS1);
      elemData->fGammaFactor   = 100.0*electron_mass_c2/z13;
      elemData->fEpsilonFactor = 100.0*electron_mass_c2/z23;
      gElementData[izet] = elemData;
    }
  }
}

void G4eBremsstrahlungRelModel::ComputeLPMfunctions(G4double& funcXiS,
                                                    G4double& funcGS,
                                                    G4double& funcPhiS,
	                                                  const G4double egamma)
{
  static const G4double sqrt2 = std::sqrt(2.);
  const G4double    redegamma = egamma/fPrimaryTotalEnergy;
  const G4double    varSprime = std::sqrt(0.125*redegamma*fLPMEnergy/
                                ((1.0-redegamma)*fPrimaryTotalEnergy));
  const ElementData* elDat    = gElementData[fCurrentIZ];
  const G4double varS1        = elDat->fVarS1;
  const G4double condition    = sqrt2*varS1;
  G4double funcXiSprime = 2.0;
  if (varSprime > 1.0) {
    funcXiSprime = 1.0;
  } else if (varSprime > condition) {
    const G4double ilVarS1Cond = elDat->fILVarS1Cond;
    const G4double funcHSprime = G4Log(varSprime)*ilVarS1Cond;
    funcXiSprime = 1.0 + funcHSprime - 0.08*(1.0-funcHSprime)*funcHSprime
                                      *(2.0-funcHSprime)*ilVarS1Cond;
  }
  const G4double varS    = varSprime/std::sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal
  const G4double varShat = varS*(1.0+fDensityCorr/(egamma*egamma));
  funcXiS = 2.0;
  if (varShat > 1.0) {
    funcXiS = 1.0;
  } else if (varShat > varS1) {
    funcXiS = 1.0+G4Log(varShat)*elDat->fILVarS1;
  }
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  //ComputeLPMGsPhis(funcGS, funcPhiS, varShat);
  //
  //MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  if (funcXiS*funcPhiS > 1. || varShat > 0.57) {
    funcXiS=1./funcPhiS;
  }
}

void G4eBremsstrahlungRelModel::ComputeLPMGsPhis(G4double& funcGS,
                                                 G4double& funcPhiS,
                                                 const G4double varShat)
{
  if (varShat < 0.01) {
    funcPhiS = 6.0*varShat*(1.0-CLHEP::pi*varShat);
    funcGS   = 12.0*varShat-2.0*funcPhiS;
  } else {
    const G4double varShat2 = varShat*varShat;
    const G4double varShat3 = varShat*varShat2;
    const G4double varShat4 = varShat2*varShat2;
    // use Stanev approximation: for \psi(s) and compute G(s)
    if (varShat < 0.415827) {
      funcPhiS = 1.0-G4Exp(-6.0*varShat*(1.0+varShat*(3.0-CLHEP::pi))
                + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      // 1-\exp \left\{-4s-\frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right\}
      const G4double funcPsiS = 1.0 - G4Exp(-4.0*varShat
                               - 8.0*varShat2/(1.0+3.936*varShat+4.97*varShat2
                               - 0.05*varShat3 + 7.5*varShat4));
      // G(s) = 3 \psi(s) - 2 \phi(s)
      funcGS = 3.0*funcPsiS - 2.0*funcPhiS;
    } else if (varShat<1.55) {
      funcPhiS = 1.0-G4Exp(-6.0*varShat*(1.0+varShat*(3.0-CLHEP::pi))
                + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      const G4double dum0  = -0.160723          + 3.755030*varShat
                             -1.798138*varShat2 + 0.672827*varShat3
                             -0.120772*varShat4;
      funcGS = std::tanh(dum0);
    } else {
      funcPhiS = 1.0-0.011905/varShat4;
      if (varShat<1.9156) {
        const G4double dum0 = -0.160723          + 3.755030*varShat
                              -1.798138*varShat2 + 0.672827*varShat3
                              -0.120772*varShat4;
        funcGS = std::tanh(dum0);
      } else {
        funcGS   = 1.0-0.023065/varShat4;
      }
    }
  }
}

// s goes up to 2 with ds = 0.01 to be the default bining
void G4eBremsstrahlungRelModel::InitLPMFunctions()
{
  if (!gLPMFuncs.fIsInitialized) {
    const G4int num = gLPMFuncs.fSLimit*gLPMFuncs.fISDelta+1;
    gLPMFuncs.fLPMFuncG.resize(num);
    gLPMFuncs.fLPMFuncPhi.resize(num);
    for (G4int i = 0; i < num; ++i) {
      const G4double sval=i/gLPMFuncs.fISDelta;
      ComputeLPMGsPhis(gLPMFuncs.fLPMFuncG[i],gLPMFuncs.fLPMFuncPhi[i],sval);
    }
    gLPMFuncs.fIsInitialized = true;
  }
}

void G4eBremsstrahlungRelModel::GetLPMFunctions(G4double& lpmGs,
                                                G4double& lpmPhis,
                                                const G4double sval)
{
  if (sval < gLPMFuncs.fSLimit) {
    G4double     val = sval*gLPMFuncs.fISDelta;
    const G4int ilow = (G4int)val;
    val    -= ilow;
    lpmGs   = (gLPMFuncs.fLPMFuncG[ilow+1]-gLPMFuncs.fLPMFuncG[ilow])*val
              + gLPMFuncs.fLPMFuncG[ilow];
    lpmPhis = (gLPMFuncs.fLPMFuncPhi[ilow+1]-gLPMFuncs.fLPMFuncPhi[ilow])*val
              + gLPMFuncs.fLPMFuncPhi[ilow];
  } else {
    G4double ss = sval*sval;
    ss *= ss;
    lpmPhis = 1.0-0.01190476/ss;
    lpmGs   = 1.0-0.0230655/ss;
  }
}

