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
// File name:     G4PairProductionRelModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 02.04.2009
//
// Modifications:
// 20.03.17 Change LPMconstant such that it gives suppression variable 's'
//          that consistent to Migdal's one; fix a small bug in 'logTS1'
//          computation; suppression is consistent now with the one in the
//          brem. model (F.Hariri)
// 28-05-18 New version with improved screening function approximation, improved
//          LPM function approximation, efficiency, documentation and cleanup. 
//          Corrected call to selecting target atom in the final state sampling. 
//          (M. Novak)
//
// Class Description:
//
// Main References:
//  J.W.Motz et.al., Rev. Mod. Phys. 41 (1969) 581.
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media,
//                     Wiley, 1972.
//
// -------------------------------------------------------------------

#include "G4PairProductionRelModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4ModifiedTsai.hh"


const G4int G4PairProductionRelModel::gMaxZet = 120; 

// LPM constant: \alpha(mc^2)^2/(4\pi*\hbar c)
const G4double G4PairProductionRelModel::gLPMconstant = 
  CLHEP::fine_structure_const*CLHEP::electron_mass_c2*CLHEP::electron_mass_c2
  /(4.*CLHEP::pi*CLHEP::hbarc);

// abscissas and weights of an 8 point Gauss-Legendre quadrature 
// for numerical integration on [0,1]
const G4double G4PairProductionRelModel::gXGL[] = { 
  1.98550718e-02, 1.01666761e-01, 2.37233795e-01, 4.08282679e-01,
  5.91717321e-01, 7.62766205e-01, 8.98333239e-01, 9.80144928e-01 
};
const G4double G4PairProductionRelModel::gWGL[] = {
  5.06142681e-02, 1.11190517e-01, 1.56853323e-01, 1.81341892e-01,
  1.81341892e-01, 1.56853323e-01, 1.11190517e-01, 5.06142681e-02 
};

// elastic and inelatic radiation logarithms for light elements (where the 
// Thomas-Fermi model doesn't work): computed by using Dirac-Fock model of atom. 
const G4double G4PairProductionRelModel::gFelLowZet  [] = {
  0.0, 5.3104, 4.7935, 4.7402, 4.7112, 4.6694, 4.6134, 4.5520
};
const G4double G4PairProductionRelModel::gFinelLowZet[] = {
  0.0, 5.9173, 5.6125, 5.5377, 5.4728, 5.4174, 5.3688, 5.3236
};

// constant cross section factor
const G4double G4PairProductionRelModel::gXSecFactor = 
  4.*CLHEP::fine_structure_const*CLHEP::classic_electr_radius
  *CLHEP::classic_electr_radius;

// gamma energy limit above which LPM suppression will be applied (if the 
// fIsUseLPMCorrection flag is true)
const G4double G4PairProductionRelModel::gEgLPMActivation = 100.*CLHEP::GeV;

// special data structure per element i.e. per Z 
std::vector<G4PairProductionRelModel::ElementData*> G4PairProductionRelModel::gElementData;

// LPM supression functions evaluated at initialisation time
G4PairProductionRelModel::LPMFuncs G4PairProductionRelModel::gLPMFuncs;

// CTR
G4PairProductionRelModel::G4PairProductionRelModel(const G4ParticleDefinition*,
                                                   const G4String& nam)
  : G4VEmModel(nam), fIsUseLPMCorrection(true), fIsUseCompleteScreening(false),
  fLPMEnergy(0.), fG4Calc(G4Pow::GetInstance()), fTheGamma(G4Gamma::Gamma()),
  fTheElectron(G4Electron::Electron()), fThePositron(G4Positron::Positron()),
  fParticleChange(nullptr) 
{
  SetAngularDistribution(new G4ModifiedTsai());
}

// DTR
G4PairProductionRelModel::~G4PairProductionRelModel()
{
  if (IsMaster()) {
    // clear ElementData container
    for (size_t iz = 0; iz < gElementData.size(); ++iz) {
      if (gElementData[iz]) delete gElementData[iz];
    }
    gElementData.clear(); 
    // clear LPMFunctions (if any)
    if (fIsUseLPMCorrection) {
      gLPMFuncs.fLPMFuncG.clear();
      gLPMFuncs.fLPMFuncPhi.clear();
      gLPMFuncs.fIsInitialized = false;
    }
  }
}

void G4PairProductionRelModel::Initialise(const G4ParticleDefinition* p,
                                          const G4DataVector& cuts)
{
  if (IsMaster()) {
    // init element data and LPM funcs
    if (IsMaster()) {
      InitialiseElementData();
      if (fIsUseLPMCorrection) {
        InitLPMFunctions();
      }
    }
  }
  if(!fParticleChange) { fParticleChange = GetParticleChangeForGamma(); }
  if(IsMaster() && LowEnergyLimit() < HighEnergyLimit()) {
    InitialiseElementSelectors(p, cuts);
  }
}

void G4PairProductionRelModel::InitialiseLocal(const G4ParticleDefinition*,
                                               G4VEmModel* masterModel)
{
  if(LowEnergyLimit() < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

G4double G4PairProductionRelModel::ComputeXSectionPerAtom(G4double gammaEnergy, 
                                                          G4double Z)
{
  G4double xSection = 0.0;
  // check if LPM suppression needs to be used
  const G4bool   isLPM  = (fIsUseLPMCorrection && gammaEnergy>gEgLPMActivation); 
  // determine the kinematical limits (taken into account the correction due to
  // the way in which the Coulomb correction is applied i.e. avoid negative DCS)
  const G4int    iz     = std::min(gMaxZet, G4lrint(Z)); 
  const G4double eps0   = CLHEP::electron_mass_c2/gammaEnergy;
  const G4double dmax   = gElementData[iz]->fDeltaMax;
  const G4double dmin   = 4.*eps0*gElementData[iz]->fDeltaFactor;
  const G4double eps1   = 0.5 - 0.5*std::sqrt(1.-dmin/dmax);
  const G4double epsMin = std::max(eps0, eps1);
  const G4double epsMax = 0.5; // DCS is symmetric around eps=0.5
  // let Et be the total energy transferred to the e- or to the e+
  // the [Et-min, Et-max] interval will be divided into i=1,2,..,n subintervals
  // with width of dInterv = (Et-max - Et-min)/n and numerical integration will 
  // be done in each sub-inteval using the xi = (Et - Et_i-min)/dInterv variable
  // that is in [0,1]. The 8-point GL q. is used for the integration on [0,1]. 
  const G4int    numSub = 2;
  const G4double dInterv= (epsMax - epsMin)*gammaEnergy/G4double(numSub); 
  G4double minEti       = epsMin*gammaEnergy;  // Et-min i.e. Et_0-min 
  for (G4int i = 0; i < numSub; ++i) {
    for (G4int ngl = 0; ngl < 8; ++ngl) {
      const G4double Et = (minEti + gXGL[ngl]*dInterv);
      const G4double xs = isLPM ? ComputeRelDXSectionPerAtom(Et, gammaEnergy, Z)
                                : ComputeDXSectionPerAtom(Et, gammaEnergy, Z);
      xSection += gWGL[ngl]*xs;
    }
    // update minimum Et of the sub-inteval
    minEti += dInterv;
  }
  // apply corrections of variable transformation and half interval integration 
  xSection = std::max(2.*xSection*dInterv, 0.);
  return xSection;
}

// DCS WITHOUT LPM SUPPRESSION
// Computes DCS value for a given target element (Z), initial gamma energy (Eg), 
// total energy transferred to one of the e-/e+ pair(Et) WITHOUT LPM suppression
// The constant factor 4 \alpha r_0^2 Z (Z +\eta(Z)) is not included here and 
// the returned value will be differential in total energy transfer instead of 
// the eps=Et/Eg. The computed part of the DCS
// NORMAL CASE: DEFAULT STTING (i.e. fIsUseCompleteScreening = FALSE)
// ds/deps(Et,Eg,Z) = ds/deps(eps,Z) = (eps^2+(1-eps)^2)*[phi1(d)/4-ln(Z)/3-fc]
// + 2*eps(1-eps)*[phi2(d)/4-ln(Z)/3-fc]/3 where the universal (in the TF model)
// screening variable d=d(eps)=136Z^(-1/3)eps0/[eps*(1-eps)] with eps0=mc^2/Eg.  
// COMPLETE SCREENING (when d(eps) approx-equal-to 0) : NEED TO BE SET BY USER 
// ds/deps(Et,Eg,Z) = ds/deps(eps,Z) = (eps^2+(1-eps)^2+eps*(1-eps)/3)*[Lel-fc]
// -eps(1-eps)/9 where Lel=phi1(0)/4-ln(Z)/3 is the elastic(coherent) radiation
// logarithm, fc is the Coulomb correction and the relation phi2(0)/4-ln(Z)/3 =
// phi1(0)/4-1/6-ln(Z)/3 = Lel-1/6 (due to phi2(0)=phi1(0)-2/3) was used.
G4double G4PairProductionRelModel::ComputeDXSectionPerAtom(G4double pEnergy,
                                                  G4double gammaEnergy,
                                                  G4double Z)
{
  G4double xSection = 0.;
  const G4int    iz   = std::min(gMaxZet, G4lrint(Z)); 
  const G4double eps  = pEnergy/gammaEnergy;
  const G4double epsm = 1.-eps;
  const G4double dum  = eps*epsm;
  if (fIsUseCompleteScreening) {
    // complete screening: 
    const G4double Lel   = gElementData[iz]->fLradEl;
    const G4double fc    = gElementData[iz]->fCoulomb;
    xSection = (eps*eps + epsm*epsm + 2.*dum/3.)*(Lel-fc) - dum/9.;  
  } else {
    // normal case:
    const G4double eps0  = CLHEP::electron_mass_c2/gammaEnergy;
    const G4double fc    = gElementData[iz]->fCoulomb;
    const G4double lnZ13 = gElementData[iz]->fLogZ13;
    const G4double delta = gElementData[iz]->fDeltaFactor*eps0/dum;
    G4double phi1, phi2;
    ComputePhi12(delta, phi1, phi2);
    xSection =  (eps*eps + epsm*epsm)*(0.25*phi1-lnZ13-fc)
               + 2.*dum*(0.25*phi2-lnZ13-fc)/3.;     
  }
  // non-const. part of the DCS differential in total energy transfer not in eps
  // ds/dEt=ds/deps deps/dEt with deps/dEt=1/Eg 
  return std::max(xSection, 0.0)/gammaEnergy;
}

// DCS WITH POSSIBLE LPM SUPPRESSION
// Computes DCS value for a given target element (Z), initial gamma energy (Eg), 
// total energy transferred to one of the e-/e+ pair(Et) WITH LPM suppression. 
// For a given Z, the LPM suppression will depend on the material through the 
// LMP-Energy. This will determine the suppression variable s and the LPM sup-
// pression functions xi(s), fi(s) and G(s). 
// The constant factor 4 \alpha r_0^2 Z (Z +\eta(Z)) is not included here and 
// the returned value will be differential in total energy transfer instead of 
// the eps=Et/Eg. The computed part of the DCS
// NORMAL CASE: DEFAULT STTING (i.e. fIsUseCompleteScreening = FALSE)
// ds/deps(Et,Eg,Z)=ds/deps(eps,Z) = xi(s)*{ (eps^2+(1-eps)^2)*[2fi(s)/3+G(s)/3]
// *[phi1(d)/4-ln(Z)/3-fc] + 2*eps(1-eps)*G(s)*[phi2(d)/4-ln(Z)/3-fc]/3 } where 
// the universal (in the TF model) screening variable d=d(eps)=136Z^(-1/3)eps0
// /[eps*(1-eps)] with eps0=mc^2/Eg.  
// COMPLETE SCREENING (when d(eps) approx-equal-to 0) : NEED TO BE SET BY USER 
// ds/deps(Et,Eg,Z) = ds/deps(eps,Z) = xi(s)*{ [Lel-fc]*[ (eps^2+(1-eps)^2+eps
// *(1-eps)/3)*2fi(s)/3 + G(s)/3] - eps(1-eps)*G(s)/9 } 
// Note, that when the LPM suppression is absent i.e. xi(s)=fi(s)=G(s)=1, both
// the normal and the complete screening DCS give back the NO-LMP case above. 
G4double G4PairProductionRelModel::ComputeRelDXSectionPerAtom(G4double pEnergy,
                                                           G4double gammaEnergy,
                                                           G4double Z)
{
  G4double xSection = 0.;
  const G4int    iz   = std::min(gMaxZet, G4lrint(Z)); 
  const G4double eps  = pEnergy/gammaEnergy;
  const G4double epsm = 1.-eps;
  const G4double dum  = eps*epsm;
  // evaluate LPM suppression functions
  G4double fXiS, fGS, fPhiS;
  ComputeLPMfunctions(fXiS, fGS, fPhiS, eps, gammaEnergy, iz);   
  if (fIsUseCompleteScreening) {
    // complete screening: 
    const G4double Lel   = gElementData[iz]->fLradEl;
    const G4double fc    = gElementData[iz]->fCoulomb;
    xSection = (Lel-fc)*((eps*eps+epsm*epsm)*2.*fPhiS + fGS)/3. - dum*fGS/9.;  
  } else {
    // normal case:
    const G4double eps0  = CLHEP::electron_mass_c2/gammaEnergy;
    const G4double fc    = gElementData[iz]->fCoulomb;
    const G4double lnZ13 = gElementData[iz]->fLogZ13;
    const G4double delta = gElementData[iz]->fDeltaFactor*eps0/dum;
    G4double phi1, phi2;
    ComputePhi12(delta, phi1, phi2);
    xSection =  (eps*eps + epsm*epsm)*(2.*fPhiS+fGS)*(0.25*phi1-lnZ13-fc)/3.
               + 2.*dum*fGS*(0.25*phi2-lnZ13-fc)/3.;     
  }
  // non-const. part of the DCS differential in total energy transfer not in eps
  // ds/dEt=ds/deps deps/dEt with deps/dEt=1/Eg 
  return std::max(fXiS*xSection, 0.0)/gammaEnergy;
}

G4double
G4PairProductionRelModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
      G4double gammaEnergy, G4double Z, G4double, G4double, G4double)
{
  G4double crossSection = 0.0 ;
  // check kinematical limit
  if ( gammaEnergy <= 2.0*electron_mass_c2 ) { return crossSection; }
  // Computes the cross section with or without LPM suppression depending on 
  // settings (by default with if the gamma energy is above a given threshold) 
  // and using or not using complete sreening approximation (by default not).
  // Only the dependent part is computed in the numerical integration of the DCS
  // i.e. the result must be multiplied here with 4 \alpha r_0^2 Z(Z+\eta(Z))
  crossSection = ComputeXSectionPerAtom(gammaEnergy, Z);
  // apply the constant factors: 
  // - eta(Z) is a correction to account interaction in the field of e-
  // - gXSecFactor = 4 \alpha r_0^2
  const G4int iz     = std::min(gMaxZet, G4lrint(Z));
  const G4double eta = gElementData[iz]->fEtaValue; 
  crossSection      *= gXSecFactor*Z*(Z+eta);
  // final protection
  return std::max(crossSection, 0.);
}

void G4PairProductionRelModel::SetupForMaterial(const G4ParticleDefinition*,
                                                const G4Material* mat, G4double)
{
  fLPMEnergy = mat->GetRadlen()*gLPMconstant;
}

void
G4PairProductionRelModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* aDynamicGamma,
                                            G4double,
                                            G4double)
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.
{
  const G4Material* mat         = couple->GetMaterial();
  const G4double    gammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4double    eps0        = CLHEP::electron_mass_c2/gammaEnergy ;
  //
  // check kinematical limit: gamma energy(Eg) must be at least 2 e- rest mass
  // (but the model should be used at higher energies above 100 MeV)
  if (eps0 > 0.5) { return; }
  // 
  // select target atom of the material
  const G4Element* anElement = SelectRandomAtom(couple, fTheGamma, gammaEnergy);
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  //
  // 'eps' is the total energy transferred to one of the e-/e+ pair in initial
  // gamma energy units Eg. Since the corresponding DCS is symmetric on eps=0.5,
  // the kinematical limits for eps0=mc^2/Eg <= eps <= 0.5 
  //
  // The Coulomb factor for the target element (Z) (Eg>50 MeV is assumed)
  // F(Z) = 8*ln(Z)/3 + 8*fc(Z)
  //
  // The screening variable 'delta(eps)' = 136*Z^{-1/3}*eps0/[eps(1-eps)]
  // Due to the Coulomb correction, the DCS can go below zero even at 
  // kinematicaly allowed eps > eps0 values. In order to exclude this eps 
  // range with negative DCS, the minimum eps value will be set to eps_min = 
  // max[eps0, epsp] with epsp is the solution of SF(delta(epsp)) - F(Z)/2 = 0 
  // with SF being the screening function (SF1=SF2 at high value of delta). 
  // The solution is epsp = 0.5 - 0.5*sqrt[ 1 - 4*136*Z^{-1/3}eps0/deltap] 
  // with deltap = Exp[(42.038-F(Z))/8.29]-0.958. So the limits are:
  // - when eps=eps_max = 0.5            => delta_min = 136*Z^{-1/3}*eps0/4
  // - epsp = 0.5 - 0.5*sqrt[ 1 - delta_min/deltap]
  // - and eps_min = max[eps0, epsp]     
  const G4int    iZet        = std::min(gMaxZet, anElement->GetZasInt());   
  const G4double deltaFactor = gElementData[iZet]->fDeltaFactor*eps0;
  const G4double deltaMin    = 4.*deltaFactor;
  const G4double deltaMax    = gElementData[iZet]->fDeltaMax;
  // compute the limits of eps
  const G4double epsp        = 0.5 - 0.5*std::sqrt(1. - deltaMin/deltaMax) ;
  const G4double epsMin      = std::max(eps0,epsp);
  const G4double epsRange    = 0.5 - epsMin;
  const G4double   FZ        = 8.*(gElementData[iZet]->fLogZ13 + 
                                   gElementData[iZet]->fCoulomb); 
  //
  // sample the energy rate (eps) of the created electron (or positron)
  G4double F10, F20;
  ScreenFunction12(deltaMin, F10, F20); 
  F10 -= FZ;
  F20 -= FZ; 
  const G4double NormF1   = std::max(F10 * epsRange * epsRange, 0.); 
  const G4double NormF2   = std::max(1.5 * F20                , 0.);
  const G4double NormCond = NormF1/(NormF1 + NormF2); 
  // check if LPM correction is active
  const G4bool isLPM = (fIsUseLPMCorrection && gammaEnergy>gEgLPMActivation);
  fLPMEnergy = mat->GetRadlen()*gLPMconstant;
  // we will need 3 uniform random number for each trial of sampling 
  G4double rndmv[3];
  G4double greject = 0.;
  G4double eps;
  do {
    rndmEngine->flatArray(3, rndmv);
    if (NormCond > rndmv[0]) {
      eps = 0.5 - epsRange * fG4Calc->A13(rndmv[1]);
      const G4double delta = deltaFactor/(eps*(1.-eps));
      if (isLPM) {
        G4double lpmXiS, lpmGS, lpmPhiS, phi1, phi2;
        ComputePhi12(delta, phi1, phi2);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, eps, gammaEnergy, iZet); 
        greject = lpmXiS*((2.*lpmPhiS+lpmGS)*phi1-lpmGS*phi2-lpmPhiS*FZ)/F10;
      } else {
        greject = (ScreenFunction1(delta)-FZ)/F10;
      }
    } else {
      eps = epsMin + epsRange*rndmv[1];
      const G4double delta = deltaFactor/(eps*(1.-eps));
      if (isLPM) {
        G4double lpmXiS, lpmGS, lpmPhiS, phi1, phi2;
        ComputePhi12(delta, phi1, phi2);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, eps, gammaEnergy, iZet); 
        greject = lpmXiS*( (lpmPhiS+0.5*lpmGS)*phi1 + 0.5*lpmGS*phi2 
                           -0.5*(lpmGS+lpmPhiS)*FZ )/F20;
      } else {
        greject = (ScreenFunction2(delta)-FZ)/F20;
      }
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (greject < rndmv[2]);
  //  end of eps sampling
  //
  // select charges randomly
  G4double eTotEnergy, pTotEnergy;
  if (rndmEngine->flat() > 0.5) {
    eTotEnergy = (1.-eps)*gammaEnergy;
    pTotEnergy = eps*gammaEnergy;
  } else {
    pTotEnergy = (1.-eps)*gammaEnergy;
    eTotEnergy = eps*gammaEnergy;
  }
  //
  // sample pair kinematics
  // 
  const G4double eKinEnergy = std::max(0.,eTotEnergy - CLHEP::electron_mass_c2);
  const G4double pKinEnergy = std::max(0.,pTotEnergy - CLHEP::electron_mass_c2);
  //
  G4ThreeVector eDirection, pDirection;
  //
  GetAngularDistribution()->SamplePairDirections(aDynamicGamma, 
						 eKinEnergy, pKinEnergy,
                                                 eDirection, pDirection);
  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1= new G4DynamicParticle(
                     fTheElectron,eDirection,eKinEnergy);

  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2= new G4DynamicParticle(
                     fThePositron,pDirection,pKinEnergy);
  // Fill output vector
  fvect->push_back(aParticle1);
  fvect->push_back(aParticle2);
  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

// should be called only by the master and at initialisation
void G4PairProductionRelModel::InitialiseElementData() 
{
  G4int size = gElementData.size();
  if (size < gMaxZet+1) {
    gElementData.resize(gMaxZet+1, nullptr);
  }
  // create for all elements that are in the detector
  const G4ElementTable* elemTable = G4Element::GetElementTable();
  size_t numElems = (*elemTable).size();
  for (size_t ie = 0; ie < numElems; ++ie) {
    const G4Element* elem = (*elemTable)[ie];
    const G4int        iz = std::min(gMaxZet, elem->GetZasInt());
    if (!gElementData[iz]) { // create it if doesn't exist yet
      const G4double logZ13 = elem->GetIonisation()->GetlogZ3();
      const G4double Z13    = elem->GetIonisation()->GetZ3();
      const G4double fc     = elem->GetfCoulomb(); 
      const G4double FZ     = 8.*(logZ13 + fc);
      G4double       Fel;
      G4double       Finel;
      if (iz<5) {  // use data from Dirac-Fock atomic model
        Fel   = gFelLowZet[iz];
        Finel = gFinelLowZet[iz];
      } else {     // use the results of the Thomas-Fermi-Moliere model 
        Fel   = G4Log(184.)  -    logZ13;
        Finel = G4Log(1194.) - 2.*logZ13;
      }
      ElementData* elD     = new ElementData(); 
      elD->fLogZ13         = logZ13;
      elD->fCoulomb        = fc;
      elD->fLradEl         = Fel;
      elD->fDeltaFactor    = 136./Z13;
      elD->fDeltaMax       = G4Exp((42.038 - FZ)/8.29) - 0.958;
      elD->fEtaValue       = Finel/(Fel-fc);
      elD->fLPMVarS1Cond   = std::sqrt(2.)*Z13*Z13/(184.*184.);
      elD->fLPMILVarS1Cond = 1./G4Log(elD->fLPMVarS1Cond);
      gElementData[iz]   = elD;
    }
  }
}

// s goes up to 2 with ds = 0.01 be default 
void G4PairProductionRelModel::InitLPMFunctions() {
  if (!gLPMFuncs.fIsInitialized) {
    const G4int num = gLPMFuncs.fSLimit*gLPMFuncs.fISDelta+1;
    gLPMFuncs.fLPMFuncG.resize(num);
    gLPMFuncs.fLPMFuncPhi.resize(num);
    for (G4int i=0; i<num; ++i) {
      const G4double sval = i/gLPMFuncs.fISDelta;
      ComputeLPMGsPhis(gLPMFuncs.fLPMFuncG[i],gLPMFuncs.fLPMFuncPhi[i],sval);
    }
    gLPMFuncs.fIsInitialized = true;
  }
}

// used only at initialisation time
void G4PairProductionRelModel::ComputeLPMGsPhis(G4double &funcGS, G4double &funcPhiS, const G4double varShat) {
  if (varShat < 0.01) {
    funcPhiS = 6.0*varShat*(1.0-CLHEP::pi*varShat);
    funcGS   = 12.0*varShat-2.0*funcPhiS;
  } else {
    const G4double varShat2 = varShat*varShat;
    const G4double varShat3 = varShat*varShat2;
    const G4double varShat4 = varShat2*varShat2;
    if (varShat < 0.415827397755) { // Stanev ap.: for \psi(s) and compute G(s)
      funcPhiS = 1.0-G4Exp( -6.0*varShat*(1.0+varShat*(3.0-CLHEP::pi)) 
                           + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      // 1-\exp \left\{-4s-\frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right\}
      const G4double funcPsiS = 1.0-G4Exp( -4.0*varShat - 8.0*varShat2/(1.0 
                     + 3.936*varShat+4.97*varShat2-0.05*varShat3+7.5*varShat4));
      // G(s) = 3 \psi(s) - 2 \phi(s)
      funcGS = 3.0*funcPsiS - 2.0*funcPhiS;
    } else if (varShat < 1.55) {
      funcPhiS = 1.0-G4Exp( -6.0*varShat*(1.0+varShat*(3.0-CLHEP::pi)) 
                           + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      const G4double dum0  = -0.16072300849123999+3.7550300067531581*varShat
                             -1.7981383069010097 *varShat2
                             +0.67282686077812381*varShat3
                             -0.1207722909879257 *varShat4;
      funcGS = std::tanh(dum0);
    } else {
      funcPhiS = 1.0-0.01190476/varShat4;
      if (varShat < 1.9156) {
        const G4double dum0 = -0.16072300849123999+3.7550300067531581*varShat
                              -1.7981383069010097 *varShat2
                              +0.67282686077812381*varShat3
                              -0.1207722909879257 *varShat4;
        funcGS = std::tanh(dum0);
      } else {
        funcGS   = 1.0-0.0230655/varShat4;
      }
    }
  }
}

// used at run-time to get some pre-computed LPM function values
void G4PairProductionRelModel::GetLPMFunctions(G4double &lpmGs, 
                                               G4double &lpmPhis, 
                                               const G4double sval) {
  if (sval < gLPMFuncs.fSLimit) {
    G4double     val = sval*gLPMFuncs.fISDelta;
    const G4int ilow = (G4int)val;
    val    -= ilow;
    lpmGs   =  (gLPMFuncs.fLPMFuncG[ilow+1]-gLPMFuncs.fLPMFuncG[ilow])*val 
             + gLPMFuncs.fLPMFuncG[ilow];
    lpmPhis =  (gLPMFuncs.fLPMFuncPhi[ilow+1]-gLPMFuncs.fLPMFuncPhi[ilow])*val 
             + gLPMFuncs.fLPMFuncPhi[ilow];
  } else {
    G4double ss = sval*sval;
    ss *= ss;
    lpmPhis = 1.0-0.01190476/ss;
    lpmGs   = 1.0-0.0230655/ss;
  }
}

void G4PairProductionRelModel::ComputeLPMfunctions(G4double &funcXiS, 
                 G4double &funcGS, G4double &funcPhiS, const G4double eps, 
                 const G4double egamma, const G4int izet)
{
  //  1. y = E_+/E_{\gamma} with E_+ being the total energy transfered 
  //                        to one of the e-/e+ pair
  //  s' = \sqrt{ \frac{1}{8} \frac{1}{y(1-y)} \frac{E^{KL}_{LPM}}{E_{\gamma}} }
  const G4double varSprime = std::sqrt(0.125*fLPMEnergy/(eps*egamma*(1.0-eps)));
  const G4double condition = gElementData[izet]->fLPMVarS1Cond;
  funcXiS = 2.0;
  if (varSprime > 1.0) {
    funcXiS = 1.0;
  } else if (varSprime > condition) {
    const G4double dum = gElementData[izet]->fLPMILVarS1Cond;
    const G4double funcHSprime = G4Log(varSprime)*dum;
    funcXiS =  1.0 + funcHSprime 
             - 0.08*(1.0-funcHSprime)*funcHSprime*(2.0-funcHSprime)*dum;
  }
  //  2. s=\frac{s'}{\sqrt{\xi(s')}}
  const G4double varShat = varSprime / std::sqrt(funcXiS);
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  if (funcXiS * funcPhiS > 1. || varShat > 0.57) {
    funcXiS = 1. / funcPhiS;
  }
}

