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
// ----------------------------------------------------------------------------
//
// GEANT4 Class implementation file
//
// File name:     G4GoudsmitSaundersonMscModel
//
// Author:        Mihaly Novak / (Omrane Kadri)
//
// Creation date: 20.02.2009
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 12.05.2010 O.Kadri: adding Qn1 and Qn12 as private doubles
// 18.05.2015 M. Novak provide PRELIMINARY verison of the revised class
//            This class has been revised and updated, new methods added.
//            A new version of Kawrakow-Bielajew Goudsmit-Saunderson MSC model
//            based on the screened Rutherford DCS for elastic scattering of
//            electrons/positrons has been introduced[1,2]. The corresponding MSC
//            angular distributions over a 2D parameter grid have been recomputed
//            and the CDFs are now stored in a variable transformed (smooth) form[2,3]
//            together with the corresponding rational interpolation parameters.
//            These angular distributions are handled by the new
//            G4GoudsmitSaundersonTable class that is responsible to sample if
//            it was no, single, few or multiple scattering case and delivers the
//            angular deflection (i.e. cos(theta) and sin(theta)).
//            Two screening options are provided:
//             - if fIsUsePWATotalXsecData=TRUE i.e. SetOptionPWAScreening(TRUE)
//               was called before initialisation: screening parameter value A is
//               determined such that the first transport coefficient G1(A)
//               computed according to the screened Rutherford DCS for elastic
//               scattering will reproduce the one computed from the PWA elastic
//               and first transport mean free paths[4].
//             - if fIsUsePWATotalXsecData=FALSE i.e. default value or
//               SetOptionPWAScreening(FALSE) was called before initialisation:
//               screening parameter value A is computed according to Moliere's
//               formula (by using material dependent parameters \chi_cc2 and b_c
//               precomputed for each material used at initialization in
//               G4GoudsmitSaundersonTable) [3]
//            Elastic and first trasport mean free paths are used consistently.
//            The new version is self-consistent, several times faster, more
//            robust and accurate compared to the earlier version.
//            Spin effects as well as a more accurate energy loss correction and
//            computations of Lewis moments will be implemented later on.
//            [1] A.F.Bielajew, NIMB 111 (1996) 195-208
//            [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//            [3] I.Kawrakow, E.Mainegra-Hing, D.W.O.Rogers, F.Tessier,B.R.B.Walters,
//                NRCC Report PIRS-701 (2013)
//            [4] F.Salvat, A.Jablonski, C.J. Powell, CPC 165(2005) 157-190
// 02.09.2015 M. Novak: first version of new step limit is provided.
//            fUseSafetyPlus corresponds to Urban fUseSafety (default)
//            fUseDistanceToBoundary corresponds to Urban fUseDistanceToBoundary
//            fUseSafety  corresponds to EGSnrc error-free stepping algorithm
//            Range factor can be significantly higher at each case than in Urban.
// 23.08.2017 M. Novak: added corrections to account spin effects (Mott-correction).
//            It can be activated by setting the fIsMottCorrection flag to be true
//            before initialization using the SetOptionMottCorrection() public method.
//            The fMottCorrection member is responsible to handle pre-computed Mott
//            correction (rejection) functions obtained by numerically computing
//            Goudsmit-Saunderson agnular distributions based on a DCS accounting spin
//            effects and screening corrections. The DCS used to compute the accurate
//            GS angular distributions is: DCS_{cor} = DCS_{SR}x[ DCS_{R}/DCS_{Mott}] where :
//               # DCS_{SR} is the relativistic Screened-Rutherford DCS (first Born approximate
//                 solution of the Klein-Gordon i.e. relativistic Schrodinger equation =>
//                 scattering of spinless e- on exponentially screened Coulomb potential)
//                 note: the default (without using Mott-correction) GS angular distributions
//                 are based on this DCS_{SR} with Moliere's screening parameter!
//               # DCS_{R} is the Rutherford DCS which is the same as above but without
//                 screening
//               # DCS_{Mott} is the Mott DCS i.e. solution of the Dirac equation with a bare
//                 Coulomb potential i.e. scattering of particles with spin (e- or e+) on a
//                 point-like unscreened Coulomb potential
//               # moreover, the screening parameter of the DCS_{cor} was determined such that
//                 the DCS_{cor} with this corrected screening parameter reproduce the first
//                 transport cross sections obtained from the corresponding most accurate DCS
//                 (i.e. from elsepa [4])
//            Unlike the default GS, the Mott-corrected angular distributions are particle type
//            (different for e- and e+ <= the DCS_{Mott} and the screening correction) and target
//            (Z and material) dependent.
// 27.10.2017 M. Novak:
//            - Mott-correction flag is set now through the G4EmParameters
//            - new form of PWA correction to integrated quantities and screening (default)
//            - changed step limit flag conventions:
//               # fUseSafety corresponds to Urban's fUseSafety
//               # fUseDistanceToBoundary corresponds to Urban's fUseDistanceToBoundary
//               # fUseSafetyPlus corresponds to the error-free stepping algorithm
// 02.02.2018 M. Novak: implemented CrossSectionPerVolume interface method (used only for testing)
// 26.10.2025 M. Novak: the model has only its accurate stepping and boundary crossing algorithms
//            left as the only option that ensures the expected precision, especially when activating
//            its Mott correction option (that also activates the screeing and scattering power
//            corrections). The model has been used for describing e-/e+ MSC (below 100 MeV kinetic)
//            energy in the option4, Penelope and Livermore EM physics constructors since Geant4 10.6.
//
// References:
//   M. Novak: https://arxiv.org/abs/2410.13361
// -----------------------------------------------------------------------------


#include "G4GoudsmitSaundersonMscModel.hh"

#include "G4GoudsmitSaundersonTable.hh"
#include "G4GSPWACorrections.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleChangeForMSC.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4Track.hh"
#include "G4PhysicsTable.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include <fstream>


G4GoudsmitSaundersonMscModel::G4GoudsmitSaundersonMscModel(const G4String& nam)
  : G4VMscModel(nam) {
  currentMaterialIndex   = -1;
  presafety              = 0.*mm;
  //
  particle               = nullptr;
  currentKinEnergy       = 0.0;
  currentRange           = 0.0;
  currentCouple          = nullptr;
  fParticleChange        = nullptr;
  //
  // Moliere screeing parameter will be used and (by default) corrections are
  // appalied to the integrated quantities (screeing parameter, elastic mfp, first
  // and second moments) derived from the corresponding PWA quantities
  // this PWA correction is ignored if Mott-correction is set to true because
  // Mott-correction contains all these corrections as well
  fIsUsePWACorrection    = true;
  fIsUseMottCorrection   = false;
  // Use the "range-rejection" like optimisation?
  fIsUseOptimisation     = true;
  //
  fLambda0               = 0.0; // elastic mean free path
  fLambda1               = 0.0; // first transport mean free path
  fScrA                  = 0.0; // screening parameter
  fG1                    = 0.0; // first transport coef.
  //
  fMCtoScrA              = 1.0;
  fMCtoQ1                = 1.0;
  fMCtoG2PerG1           = 1.0;
  //
  fTheTrueStepLenght     = 0.;
  fTheZPathLenght        = 0.;

  fTheDisplacementVector.set(0.,0.,0.);
  fTheNewDirection.set(0.,0.,1.);

  fIsEndedUpOnBoundary   = false;
  fIsMultipleScattering  = false;
  fIsSingleScattering    = false;
  fIsNoScatteringInMSC   = false;
  fIsSimplified          = false;

  fGSTable               = nullptr;
  fPWACorrection         = nullptr;
}


G4GoudsmitSaundersonMscModel::~G4GoudsmitSaundersonMscModel() {
  if (IsMaster()) {
    if (fGSTable) {
      delete fGSTable;
      fGSTable = nullptr;
    }
    if (fPWACorrection) {
      delete fPWACorrection;
      fPWACorrection = nullptr;
    }
  }
}


void G4GoudsmitSaundersonMscModel::Initialise(const G4ParticleDefinition* p, const G4DataVector&) {
  SetParticle(p);
  InitialiseParameters(p);
  // -create GoudsmitSaundersonTable and init its Mott-correction member if
  //  Mott-correction was required
  if (IsMaster()) {
    // get the Mott-correction flag from EmParameters
    if (G4EmParameters::Instance()->UseMottCorrection()) {
      fIsUseMottCorrection = true;
    }
    // Mott-correction includes other way of PWA x-section corrections so deactivate it even if it was true
    // when Mott-correction is activated by the user
    if (fIsUseMottCorrection) {
      fIsUsePWACorrection = false;
    }
    // clear GS-table
    if (fGSTable) {
      delete fGSTable;
      fGSTable = nullptr;
    }
    // clear PWA corrections table if any
    if (fPWACorrection) {
      delete fPWACorrection;
      fPWACorrection = nullptr;
    }
    // create GS-table
    G4bool isElectron = true;
    if (p->GetPDGCharge()>0.) {
      isElectron = false;
    }
    fGSTable = new G4GoudsmitSaundersonTable(isElectron);
    // G4GSTable will be initialised:
    // - Screened-Rutherford DCS based GS angular distributions will be loaded only if they are not there yet
    // - Mott-correction will be initialised if Mott-correction was requested to be used
    fGSTable->SetOptionMottCorrection(fIsUseMottCorrection);
    // - set PWA correction (correction to integrated quantites from Dirac-PWA)
    fGSTable->SetOptionPWACorrection(fIsUsePWACorrection);
    // init
    fGSTable->Initialise(LowEnergyLimit(),HighEnergyLimit());
    // create PWA corrections table if it was requested (and not deactivated because active Mott-correction)
    if (fIsUsePWACorrection) {
      fPWACorrection = new G4GSPWACorrections(isElectron);
      fPWACorrection->Initialise();
    }
  }
  fParticleChange = GetParticleChangeForMSC(p);
}


void G4GoudsmitSaundersonMscModel::InitialiseLocal(const G4ParticleDefinition*, G4VEmModel* masterModel) {
   fGSTable               = static_cast<G4GoudsmitSaundersonMscModel*>(masterModel)->GetGSTable();
   fIsUseMottCorrection   = static_cast<G4GoudsmitSaundersonMscModel*>(masterModel)->GetOptionMottCorrection();
   fIsUsePWACorrection    = static_cast<G4GoudsmitSaundersonMscModel*>(masterModel)->GetOptionPWACorrection();
   fPWACorrection         = static_cast<G4GoudsmitSaundersonMscModel*>(masterModel)->GetPWACorrection();
}


// computes macroscopic first transport cross section: used only in testing not during mc transport
G4double G4GoudsmitSaundersonMscModel::CrossSectionPerVolume(const G4Material* mat,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double,
                                         G4double) {
  // use Moliere's screening (with Mott-corretion if it was requested)
  kineticEnergy = std::max(kineticEnergy, 10.*CLHEP::eV);
  // // total mometum square, beta2
  const G4double pt2   = kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2);
  const G4double beta2 = pt2/(pt2 + electron_mass_c2*electron_mass_c2);
  const G4int matindx  = static_cast<G4int>(mat->GetIndex());
  // get the Mott-correcton factors if Mott-correcton was requested by the user
  fMCtoScrA       = 1.0;
  fMCtoQ1         = 1.0;
  fMCtoG2PerG1    = 1.0;
  G4double scpCor = 1.0; // keep this for consistency (this method is not used in normal runs)
  if (fIsUseMottCorrection) {
    fGSTable->GetMottCorrectionFactors(G4Log(kineticEnergy), beta2, matindx, fMCtoScrA, fMCtoQ1, fMCtoG2PerG1);
    // ! no scattering power correction since the current couple is not set before this interface method is called
    // scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  } else if (fIsUsePWACorrection) {
    fPWACorrection->GetPWACorrectionFactors(G4Log(kineticEnergy), beta2, matindx, fMCtoScrA, fMCtoQ1, fMCtoG2PerG1);
    // scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  }
  // screening parameter:
  // - if Mott-corretioncorrection: the Screened-Rutherford times Mott-corretion DCS with this
  //   screening parameter gives back the (elsepa) PWA first transport cross section
  // - if PWA correction: he Screened-Rutherford DCS with this screening parameter
  //   gives back the (elsepa) PWA first transport cross section
  const G4double bc = fGSTable->GetMoliereBc(matindx);
  fScrA    = fGSTable->GetMoliereXc2(matindx)/(4.0*pt2*bc)*fMCtoScrA;
  // elastic mean free path in Geant4 internal lenght units: the neglected (1+screening parameter) term is corrected
  // (if Mott-corretion: the corrected screening parameter is used for this (1+A) correction + Moliere b_c is also
  // corrected with the screening parameter correction)
  fLambda0 = beta2*(1.+fScrA)*fMCtoScrA/(bc*scpCor);
  // first transport coefficient (if Mott-corretion: the corrected screening parameter is used (it will be fully
  // consistent with the one used during the pre-computation of the Mott-correted GS angular distributions))
  fG1      = 2.0*fScrA*((1.0+fScrA)*G4Log(1.0/fScrA+1.0)-1.0);
  // first transport mean free path
  fLambda1 = fLambda0/fG1;
  // return with the macroscopic first transport cross section
  return fLambda1 > 0.0 ? 1.0/fLambda1 : 1.0E+20;
}


// gives back the first transport mean free path in internal G4 units
G4double
G4GoudsmitSaundersonMscModel::GetTransportMeanFreePath(const G4ParticleDefinition* /*partdef*/,
                                                       G4double kineticEnergy) {
  // use Moliere's screening (with Mott-corretion if it was requested)
  kineticEnergy = std::max(kineticEnergy, 10.0*CLHEP::eV);
  // total mometum square, beta2
  const G4double pt2     = kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2);
  const G4double beta2   = pt2/(pt2 + electron_mass_c2*electron_mass_c2);
  const G4int    matindx = static_cast<G4int>(currentCouple->GetMaterial()->GetIndex());
  // get the Mott and scattering power correcton factors
  // (if Mott-correcton was activated)
  fMCtoScrA       = 1.0;
  fMCtoQ1         = 1.0;
  fMCtoG2PerG1    = 1.0;
  G4double scpCor = 1.0;
  if (fIsUseMottCorrection) {
    fGSTable->GetMottCorrectionFactors(G4Log(kineticEnergy), beta2, matindx, fMCtoScrA, fMCtoQ1, fMCtoG2PerG1);
    scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  } else if (fIsUsePWACorrection) {
    fPWACorrection->GetPWACorrectionFactors(G4Log(kineticEnergy), beta2, matindx, fMCtoScrA, fMCtoQ1, fMCtoG2PerG1);
    // scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  }
  // screening parameter:
  // - if Mott-corretioncorrection: the Screened-Rutherford times Mott-corretion DCS with this
  //   screening parameter gives back the (elsepa) PWA first transport cross section
  // - if PWA correction: he Screened-Rutherford DCS with this screening parameter
  //   gives back the (elsepa) PWA first transport cross section
  const G4double bc = fGSTable->GetMoliereBc(matindx);
  fScrA    = fGSTable->GetMoliereXc2(matindx)/(4.0*pt2*bc)*fMCtoScrA;
  // elastic mean free path in Geant4 internal lenght units: the neglected (1+screening parameter) term is corrected
  // (if Mott-corretion: the corrected screening parameter is used for this (1+A) correction + Moliere b_c is also
  // corrected with the screening parameter correction)
  fLambda0 = beta2*(1.+fScrA)*fMCtoScrA/(bc*scpCor);
  // first transport coefficient (if Mott-corretion: the corrected screening parameter is used (it will be fully
  // consistent with the one used during the pre-computation of the Mott-correted GS angular distributions))
  fG1      = 2.0*fScrA*((1.0+fScrA)*G4Log(1.0/fScrA+1.0)-1.0);
  // first transport mean free path
  fLambda1 = fLambda0/fG1;

  return fLambda1;
}


G4double
G4GoudsmitSaundersonMscModel::GetTransportMeanFreePathOnly(const G4ParticleDefinition* /*partdef*/,
                                                           G4double kineticEnergy) {
  // use Moliere's screening (with Mott-corretion if it was requested)
  kineticEnergy = std::max(kineticEnergy, 10.0*CLHEP::eV);
  // total mometum square, beta2
  const G4double pt2     = kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2);
  const G4double beta2   = pt2/(pt2 + electron_mass_c2*electron_mass_c2);
  const G4int    matindx = static_cast<G4int>(currentCouple->GetMaterial()->GetIndex());
  // get the Mott and scattering power correcton factors
  // (if Mott-correcton was activated)
  G4double mctoScrA    = 1.0;
  G4double mctoQ1      = 1.0;
  G4double mctoG2PerG1 = 1.0;
  G4double scpCor      = 1.0;
  if (fIsUseMottCorrection) {
    fGSTable->GetMottCorrectionFactors(G4Log(kineticEnergy), beta2, matindx, mctoScrA, mctoQ1, mctoG2PerG1);
    scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  } else if (fIsUsePWACorrection) {
    fPWACorrection->GetPWACorrectionFactors(G4Log(kineticEnergy), beta2, matindx, mctoScrA, mctoQ1, mctoG2PerG1);
    // scpCor = fGSTable->ComputeScatteringPowerCorrection(currentCouple, kineticEnergy);
  }
  // screening parameter (with corrections)
  const G4double bc      = fGSTable->GetMoliereBc(matindx);
  const G4double scrA    = fGSTable->GetMoliereXc2(matindx)/(4.0*pt2*bc)*mctoScrA;
  // elastic mean free path (under the SR DCS with Mott and scattering power corrections)
  const G4double lambda0 = beta2*(1.+scrA)*mctoScrA/(bc*scpCor);
  // first transport coefficient (under the SR DCS)
  const G4double g1      = 2.0*scrA*((1.0+scrA)*G4Log(1.0/scrA+1.0)-1.0);
  // first transport mean free path
  const G4double lambda1 = lambda0/g1;

  return lambda1;
}


void G4GoudsmitSaundersonMscModel::StartTracking(G4Track* track) {
  SetParticle(track->GetDynamicParticle()->GetDefinition());
}


G4double G4GoudsmitSaundersonMscModel::ComputeTruePathLengthLimit(const G4Track& track,
			                                                            G4double& currentMinimalStep) {
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  G4StepPoint* sp             = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus     = sp->GetStepStatus();
  currentCouple               = track.GetMaterialCutsCouple();
  SetCurrentCouple(currentCouple);
  currentMaterialIndex        = (G4int)currentCouple->GetMaterial()->GetIndex();
  currentKinEnergy            = dp->GetKineticEnergy();
  currentRange                = GetRange(particle,currentKinEnergy,currentCouple,
                                         dp->GetLogKineticEnergy());
  // elastic and first transport mfp, screening parameter and G1 are also set
  // (Mott-correction will be used if it was requested by the user)
  fLambda1 = GetTransportMeanFreePath(particle,currentKinEnergy);
  // Set initial values:
  //  : lengths are initialised to currentMinimalStep  which is the true, minimum
  //    step length from all other physics
  fTheTrueStepLenght    = currentMinimalStep;
  fTheZPathLenght       = currentMinimalStep;  // will need to be converted
  fTheDisplacementVector.set(0.0, 0.0, 0.0);
  fTheNewDirection.set(0.0, 0.0, 1.0);

  // Multiple scattering needs to be sampled ?
  fIsMultipleScattering = false;
  // Single scattering needs to be sampled ?
  fIsSingleScattering   = false;
  // Was zero deflection in multiple scattering sampling ?
  fIsNoScatteringInMSC  = false;

  // === The error-free stepping and boundary crossing
  presafety =  ComputeSafety(sp->GetPosition(), fTheTrueStepLenght);
  // "range-rejection" like optimisation (i.e. deep inside the volume)
  fIsSimplified = false;
  if (fIsUseOptimisation && currentRange < presafety ) {
    // it's assumed that the e-/e+ never leaves the volume
    // --> simplified tracking, no MSC step limit, no displacement only deflection
    fIsSimplified = true;
    SampleMSC();
    return ConvertTrueToGeom(fTheTrueStepLenght, currentMinimalStep);
  } else {
    // the step limit that corresponds to the accurate stepping and boundary
    // crossing algorithm
    // Set skin depth = skin x elastic_mean_free_path
    const G4double skindepth = skin*fLambda0;
    // Check if we can try Single Scattering because we are within skindepth
    // distance from/to a boundary OR the current minimum true-step-length is
    // shorter than skindepth. NOTICE: the latest has only efficieny reasons
    // because the MSC angular sampling is fine for any short steps but much
    // faster to try single scattering in case of short steps.
    if ((stepStatus == fGeomBoundary) || (presafety < skindepth) || (fTheTrueStepLenght < skindepth)) {
      //Try single scattering:
      // - sample distance to next single scattering interaction (sslimit)
      // - compare to current minimum length
      //      == if sslimit is the shorter:
      //          - set the step length to sslimit
      //          - indicate that single scattering needs to be done
      //      == else : nothing to do
      //- in both cases, the step length was very short so geometrical and
      //  true path length are the same
      const G4double sslimit = -1.*fLambda0*G4Log(G4UniformRand());
      // compare to current minimum step length
      if (sslimit < fTheTrueStepLenght) {
        fTheTrueStepLenght  = sslimit;
        fIsSingleScattering = true;
      }
      // short step -> true step length equal to geometrical path length
      fTheZPathLenght = fTheTrueStepLenght;
      // Set that everything was done in step-limit phase (so no MSC call later)
      // We will check if we need to perform the single-scattering angular
      // sampling i.e. if single elastic scattering was the winer!
    } else {
      // After checking we know that we cannot try single scattering so we will
      // need to make an MSC step
      // Indicate that we need to make and MSC step.
      fIsMultipleScattering = true;
      // limit from range factor (keep this to allow stricter control)
      fTheTrueStepLenght    = std::min(fTheTrueStepLenght, facrange*currentRange);
      // never let the particle go further than the safety if we are out of the skin
      // if we are here we are out of the skin, presafety > 0.
      if (fTheTrueStepLenght > presafety) {
        fTheTrueStepLenght = std::min(fTheTrueStepLenght, presafety);
      }
      // make sure that we are still within the aplicability of condensed histry model
      // i.e. true step length is not longer than 1/2 first transport mean free path.
      // We schould take into account energy loss along the 1/2*lambda_transport1
      // step length as well but we don't. So let it just 0.5*lambda_transport1
      fTheTrueStepLenght = std::min(fTheTrueStepLenght, fLambda1*0.5);
    }
  }
  // === end of step limit
  // performe single scattering/multiple scattering
  if (fIsMultipleScattering) {
    // sample multiple scattering (including zero, one, few and multiple)
    // the final position will also be calculated using the accurate LLCA
    // leading to the transport vector, its projection to the current direction
    // (`fTheZPathLenght`) and the perpendicular plane (`fTheDisplacementVector`)
    // The MSC scattering able is also used to set the  `fTheNewDirection`
    SampleMSC();
  } else if (fIsSingleScattering) {
    // Sample single scattering angular deflection and claculate the new
    // direction that is set in `fTheNewDirection`
    // Note: `fTheZPathLenght` set to `true-step-length` while `fTheDisplacementVector`
    //       stays `(0, 0, 0)` vector in this case.
    const G4double lekin  = G4Log(currentKinEnergy);
    const G4double pt2    = currentKinEnergy*(currentKinEnergy+2.0*CLHEP::electron_mass_c2);
    const G4double beta2  = pt2/(pt2+CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);
    G4double cost   = fGSTable->SingleScattering(1.0, fScrA, lekin, beta2, currentMaterialIndex);
    cost = std::max(-1.0, std::min(cost, 1.0));
    // compute sint
    const G4double dum    = 1.0 - cost;
    const G4double sint   = std::sqrt(dum*(2.0 - dum));
    const G4double phi    = CLHEP::twopi*G4UniformRand();
    const G4double sinPhi = std::sin(phi);
    const G4double cosPhi = std::cos(phi);
    fTheNewDirection.set(sint*cosPhi, sint*sinPhi, cost);
  } // otherwise: single scattering case because within skin but not elastic
    //            scattering proposed the shortest step (but most likely geometry)
 return ConvertTrueToGeom(fTheTrueStepLenght, currentMinimalStep);
}

G4double G4GoudsmitSaundersonMscModel::ComputeGeomPathLength(G4double) {
  // convert true -> geom (before transportation still at the pre-step point)
  // (called from the step limitation ComputeTruePathLengthLimit)
  // "range-rejection" like optimisation (just give whatever geometrical length)
  if (fIsSimplified) {
    fTheZPathLenght = std::min(fTheTrueStepLenght, fLambda1);
  }
  // Note: as we computed the post-step point, transport vector and its projection
  //       along the pre-step point direction at the step limit phase in
  //      `ComputeTruePathLengthLimit` above, we just return that projection here
  return fTheZPathLenght;
}

G4double G4GoudsmitSaundersonMscModel::ComputeTrueStepLength(G4double geomStepLength) {
  // convert geom -> true (after transportation already at the post-step point)
  // Note: we computed the post-step point, transport vector and its projection
  //       along the pre-step point direction at the step limit phase in
  //      `ComputeTruePathLengthLimit` above. That computation was bsed on the
  //       actual true step length that was stored. Morover, we ensured that
  //       that post-step point will be either within the safety or used single
  //       scattering within the skin. Therefore, the post-step point is either
  //       exactly as we expected (`geomStepLength == fTheZPathLenght` as
  //       transportation could move the track with `fTheZPathLenght`) or the
  //       boundary was reached (`geomStepLength < fTheZPathLenght`) but then
  //       it was a single scattering step. The true step length is the one
  //       stored in the first case while equal to the geometrical one in the
  //       second case (as it was a single scattering step anyway).
  fIsEndedUpOnBoundary = false;
  if (geomStepLength < fTheZPathLenght) {
    // reached boundary (or very last step) otherwise
    fIsEndedUpOnBoundary = true;
    fTheZPathLenght      = geomStepLength;
    fTheTrueStepLenght   = geomStepLength;
  }
  return fTheTrueStepLenght;
}

G4ThreeVector&
G4GoudsmitSaundersonMscModel::SampleScattering(const G4ThreeVector& oldDirection, G4double) {
  // do nothing on the boundary (reached that in single scattering mode)
  if (!fIsEndedUpOnBoundary) {
    // "range-rejection" like optimisation
    if (fIsSimplified && !fIsNoScatteringInMSC) {
       fTheNewDirection.rotateUz(oldDirection);
       fParticleChange->ProposeMomentumDirection(fTheNewDirection);
       return fTheDisplacementVector;
    }
    // if single scattering mode (i.e. within skin) and it's actually happened
    if (fIsSingleScattering) {
      fTheNewDirection.rotateUz(oldDirection);
      fParticleChange->ProposeMomentumDirection(fTheNewDirection);
      return fTheDisplacementVector;
    }
    // if multiple scattering mode (i.e. out of skin) and it's actually happened
    if (fIsMultipleScattering && !fIsNoScatteringInMSC) {
       fTheNewDirection.rotateUz(oldDirection);
       fTheDisplacementVector.rotateUz(oldDirection);
       fParticleChange->ProposeMomentumDirection(fTheNewDirection);
    }
  }
  // on the boundary: reached in single scattering (or optimisation) -> displacement is (0,0,0))
  return fTheDisplacementVector;
}

void G4GoudsmitSaundersonMscModel::SampleMSC() {
  fIsNoScatteringInMSC = false;
  // Energy loss correction: calculate effective energy and step length
  // (Eqs.(77) in my notes)
  //
  // `currentKinEnergy` is the pre-step point kinetic energy `E_0`
  // (mean) energy loss (estimate only assuming a step length = `fTheTrueStepLenght`)
  const G4double eloss   = currentKinEnergy - GetEnergy(particle, currentRange - fTheTrueStepLenght, currentCouple);
  const G4double midEkin = currentKinEnergy - 0.5*eloss; // mid-step energy `\tilde{E} := E_0 - 0.5*\Delta E`
  const G4double tau     = midEkin/electron_mass_c2;     // `\tau := \tilde{E}/(mc^2)`
  const G4double tau2    = tau*tau;
  const G4double eps0    = eloss/currentKinEnergy; // energy loss fraction to pre-step energy: `\epsilon := \Delta E/E_0`
  const G4double epsm    = eloss/midEkin; // energy loss fraction to the mid-step energy: `\epsilon := \Delta E/\tilde{E}`
  const G4double epsm2   = epsm*epsm;
  const G4double effEner = midEkin*(1.0 - epsm2*(6.0 + 10.0*tau + 5.0*tau2)/(24.0*tau2 + 72.0*tau + 48.0));
  const G4double dum     = 1.0/((tau + 1.0)*(tau + 2.0));
  const G4double effStep = fTheTrueStepLenght*(1.0 - 0.166666*epsm2*dum*dum*(4.0 + tau*(6.0 + tau*(7.0 + tau*(4.0 + tau)))));
  // compute elastic mfp, first transport mfp, screening parameter, and G1 at this `E_eff`
  // (with Mott-correction if it was requested by the user)
  fLambda1 = GetTransportMeanFreePath(particle, effEner);
  // s/lambda_el i.e. mean number elastic scattering along the step
  // (using the effective step length and energy)
  const G4double lambdan = fLambda0 > 0.0 ? effStep/fLambda0 : 0.0;
  // do nothing in case of very small mean elastic scattering
  if (lambdan <= 1.0E-12) {
    fTheZPathLenght      = fTheTrueStepLenght;
    fIsNoScatteringInMSC = true;
    return;
  }
  // first moment: 2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);
  // (which is equal to Q1 = s*G1/\lambda = s/\lambda_1)
  G4double Qn1 = lambdan *fG1;
  // sample scattering angles: divide the step into two and combine at the end
  // NOTE: this shouldn't happen in normal case as we limit the step such that
  //       `s_max <= 0.5*\lambda_1` leading to `Q1_max = s/|lambda_1 = 0.5`.
  //       But we can see `Q1 > 0.5` in some corner case (e.g. due to the energy
  //       loss along the step as we assumed `\lambda_1` to be constant).
  //       Anyway, our first GS angular distribution table covers the [0.001,0.99]
  //       which is already well above the validity of the MSC theory as well as
  //       the LLCA algorithm that gives accurate post-step positions only up to
  //       Q1 < ~ 0.5. This is why we have that `s_max = 0.5*\lambda_1` step limit.
  //       (while we have a second GS angular distribution table up to Q1=7.99
  //        that is a smooth convergence to isotropic angular distributions)
  G4double cosTheta1 = 1.0, sinTheta1 = 0.0, cosTheta2 = 1.0, sinTheta2 = 0.0;
  if (0.5*Qn1 < 7.0) {
    // (expected to have Qn1 that are usually < ~0.5)
    // sample 2 scattering cost1, sint1, cost2 and sint2 for half of the step
    const G4double lekin  = G4Log(effEner);
    const G4double pt2    = effEner*(effEner + 2.0*CLHEP::electron_mass_c2);
    const G4double beta2  = pt2/(pt2 + CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);
    // backup GS angular dtr pointer (kinetic energy and delta index in case of Mott-correction)
    // if the first was an msc sampling (the same will be used if the second is also an msc step)
    G4GoudsmitSaundersonTable::GSMSCAngularDtr *gsDtr = nullptr;
    G4int mcEkinIdx    = -1;
    G4int mcDeltIdx    = -1;
    G4double transfPar = 0.0;
    G4bool isMsc = fGSTable->Sampling(0.5*lambdan, 0.5*Qn1, fScrA, cosTheta1, sinTheta1, lekin, beta2,
                                     currentMaterialIndex, &gsDtr, mcEkinIdx, mcDeltIdx, transfPar,
                                     true);
    fGSTable->Sampling(0.5*lambdan, 0.5*Qn1, fScrA, cosTheta2, sinTheta2, lekin, beta2,
                      currentMaterialIndex, &gsDtr, mcEkinIdx, mcDeltIdx, transfPar, !isMsc);
    if (cosTheta1 + cosTheta2 == 2.0) { // no scattering happened
        fTheZPathLenght      = fTheTrueStepLenght;
        fIsNoScatteringInMSC = true;
        return;
    }
  } else {
    // corner case e.g. last step or "range-rejection" like opt.--> isotropic scattering
    cosTheta1 = 1.-2.*G4UniformRand();
    sinTheta1 = std::sqrt((1.-cosTheta1)*(1.+cosTheta1));
    cosTheta2 = 1.-2.*G4UniformRand();
    sinTheta2 = std::sqrt((1.-cosTheta2)*(1.+cosTheta2));
  }
  // sample 2 azimuthal angles
  const G4double phi1    = CLHEP::twopi*G4UniformRand();
  const G4double sinPhi1 = std::sin(phi1);
  const G4double cosPhi1 = std::cos(phi1);
  const G4double phi2    = CLHEP::twopi*G4UniformRand();
  const G4double sinPhi2 = std::sin(phi2);
  const G4double cosPhi2 = std::cos(phi2);
  // compute final direction after the two scatterings (from the initial 0,0,1 dir)
  const G4double u2 = sinTheta2*cosPhi2;
  const G4double v2 = sinTheta2*sinPhi2;
  const G4double tm = cosTheta1*u2 + sinTheta1*cosTheta2;
  const G4double uf = tm*cosPhi1 - v2*sinPhi1;
  const G4double vf = tm*sinPhi1 + v2*cosPhi1;
  const G4double wf = cosTheta1*cosTheta2 - sinTheta1*u2;
  // set the new direction (still in the scattering frame that is 0,0,1)
  fTheNewDirection.set(uf, vf, wf);
  //
  // "range-rejection" like optimisation --> no dispalcement
  if (fIsSimplified) {
    return;
  }
  //
  // == Compute final position (using the accurate LLCA)
  //    [NIMB 142(1998) Eq.76-80 with energy loss correction from Med.Phys. 27 (2000))
  // apply correction to Q1
  Qn1 *=  fMCtoQ1;
  // energy loss correction (gamma, delta, eta with alpha1 and alpha2 for eloss cor.)
  // note: quantities are based on the SR DCS, gamma and delta with (Mott, screening, DPWA) corrections
  //       while the energy loss correction (alpha1, alpha2) is pure SR DCS based.
  // compute alpha1, delta + alpha2 and gamma for energy loss correction
  const G4double tp2  = tau + 2.0;
  const G4double tp1  = tau + 1.0;
  const G4double loga = G4Log(1.0 + 1.0/fScrA);
  const G4double lgf1 = loga*(1.0 + fScrA) - 1.0;
  const G4double lgf2 = loga*(1.0 + 2.0*fScrA) - 2.0;
  // energy loss correction factors \alpha_1 and 1.0-alpha1
  const G4double alpha1 = ((2.0 + tau*tp2)/tp1 - tp1/lgf1)*epsm/tp2;
  // \alpha2 (plus the -0.25.. contrib from higer order)
  const G4double alpha2 = 0.4082483*(eps0*tp1/(tp2*lgf1*lgf2) - 0.25*alpha1*alpha1);
  // then gamma and delta (with energy loss correction: \delta -> \delta + \alpha_2)
  const G4double gamma  = (6.0*fScrA*lgf2*(1.0 + fScrA)/fG1)*fMCtoG2PerG1;
  // calculate \delta -->  \delta + alpha_2
  const G4double delta  = 0.9082483 - (0.1020621 - 0.0263747*gamma)*Qn1 + alpha2;
  //
  // ready to calculate the final position:
  // sample \eta from p(eta)=2*eta i.e. P(eta) = eta_square ;-> P(eta) = rand --> eta = sqrt(rand)
  const G4double eta  = std::sqrt(G4UniformRand());
  const G4double eta0 = 0.5*(1 - eta)*(1.0 + alpha1); // \eta_0 --> \eta_0(1 + \alpha_1)
  const G4double eta1 = eta*delta;       // \eta_1
  const G4double eta2 = eta*(1.0 - delta); // \eta_2
  // calculate the post-step point: xf, yf, zf are the final position coordinates
  // divided by the (true) step length `s` and as given in NIMB 142 (1998)
  // with the energy loss correction derived in Med.Phys. 27 (2000)
  const G4double a1   = 0.5*(1 - eta)*(1.0 - alpha1);
  const G4double w1v2 = cosTheta1*v2;
  const G4double dum1 = eta1*sinTheta1;
  const G4double xf   = dum1*cosPhi1 + eta2*(cosPhi1*u2 - sinPhi1*w1v2) + a1*uf;
  const G4double yf   = dum1*sinPhi1 + eta2*(sinPhi1*u2 + cosPhi1*w1v2) + a1*vf;
  const G4double zf   = eta0 + eta1*cosTheta1 + eta2*cosTheta2          + a1*wf;
  // final position relative to the pre-step point in the scattering frame
  // (rx, ry, rz) = (xf, yf, zf)*s
  const G4double rx = xf*fTheTrueStepLenght;
  const G4double ry = yf*fTheTrueStepLenght;
  const G4double rz = zf*fTheTrueStepLenght;
  // calculate the transport distance (with protection: tr_distance <= true_step_length)
  const G4double transportDistance = std::min(std::sqrt(rx*rx + ry*ry + rz*rz), fTheTrueStepLenght);
  // set the z-path length i.e. the geometrical length to be the transport
  // distance the displacement such that the final post-step point will
  // eventually be at (rx, ry, rz) relative to the pre-step point (with rotations)
  fTheZPathLenght = transportDistance;
  fTheDisplacementVector.set(rx, ry, rz - fTheZPathLenght);
}
