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
// GEANT4 Class header file
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
// 18.05.2015 M. Novak provide PLERIMINARYY version of updated class.
//            All algorithms of the class were revised and updated, new methods added.
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
//             - if fgIsUsePWATotalXsecData=TRUE i.e. SetOptionPWAScreening(TRUE)
//               was called before initialisation: screening parameter value A is
//               determined such that the first transport coefficient G1(A)
//               computed according to the screened Rutherford DCS for elastic
//               scattering will reproduce the one computed from the PWA elastic
//               and first transport mean free paths[4].
//             - if fgIsUsePWATotalXsecData=FALSE i.e. default value or
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
// 02.02.2018 M. Novak: implemented CrossSectionPerVolume interface method (used only for testing)
//
// Class description:
//   Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened Rutherford DCS
//   for elastic scattering of e-/e+. Option, to include (Mott) correction (see above), is
//   also available now (SetOptionMottCorrection(true)). An EGSnrc like error-free stepping
//   algorithm (UseSafety) is available beyond the usual Geant4 step limitation algorithms
//   and true to geomerty and geometry to true step length computations that were adopted
//   from the Urban model[5]. The most accurate setting: error-free stepping (UseSafety)
//   with Mott-correction (SetOptionMottCorrection(true)).
//
// References:
//   [1] A.F.Bielajew, NIMB 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//   [3] I.Kawrakow, E.Mainegra-Hing, D.W.O.Rogers, F.Tessier,B.R.B.Walters, NRCC
//       Report PIRS-701 (2013)
//   [4] F.Salvat, A.Jablonski, C.J. Powell, CPC 165(2005) 157-190
//   [5] L.Urban, Preprint CERN-OPEN-2006-077 (2006)
//
// -----------------------------------------------------------------------------

#ifndef G4GoudsmitSaundersonMscModel_h
#define G4GoudsmitSaundersonMscModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "globals.hh"


class G4DataVector;
class G4ParticleChangeForMSC;
class G4LossTableManager;
class G4GoudsmitSaundersonTable;
class G4GSPWACorrections;

class G4GoudsmitSaundersonMscModel : public G4VMscModel
{
public:

  G4GoudsmitSaundersonMscModel(const G4String& nam = "GoudsmitSaunderson");

  ~G4GoudsmitSaundersonMscModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition* p, G4VEmModel* masterModel) override;

  G4ThreeVector& SampleScattering(const G4ThreeVector&, G4double safety) override;

  G4double ComputeTruePathLengthLimit(const G4Track& track, G4double& currentMinimalStep) override;

  G4double ComputeGeomPathLength(G4double truePathLength) override;

  G4double ComputeTrueStepLength(G4double geomStepLength) override;

  // method to compute first transport cross section per Volume (i.e. macroscropic first transport cross section; this 
  // method is used only for testing and not during a normal simulation) 
  G4double CrossSectionPerVolume(const G4Material*, const G4ParticleDefinition*, G4double kineticEnergy, G4double cutEnergy = 0.0, G4double maxEnergy = DBL_MAX) override;

  void     StartTracking(G4Track*) override;

  void     SampleMSC();

  G4double GetTransportMeanFreePath(const G4ParticleDefinition*, G4double);

  void SetOptionPWACorrection(G4bool opt)    { fIsUsePWACorrection = opt; }

  G4bool GetOptionPWACorrection() const      { return fIsUsePWACorrection; }

  void   SetOptionMottCorrection(G4bool opt) { fIsUseMottCorrection = opt; }

  G4bool GetOptionMottCorrection() const     { return fIsUseMottCorrection; }

  G4GoudsmitSaundersonTable* GetGSTable()          { return fGSTable; }

  G4GSPWACorrections*        GetPWACorrection()    { return fPWACorrection; }

  //  hide assignment operator
  G4GoudsmitSaundersonMscModel & operator=(const  G4GoudsmitSaundersonMscModel &right) = delete;
  G4GoudsmitSaundersonMscModel(const  G4GoudsmitSaundersonMscModel&) = delete;

private:
  inline void     SetParticle(const G4ParticleDefinition* p);

  inline G4double GetLambda(G4double);

  G4double GetTransportMeanFreePathOnly(const G4ParticleDefinition*,G4double);

  inline G4double Randomizetlimit();

private:
  CLHEP::HepRandomEngine* rndmEngineMod;
  //
  G4double currentKinEnergy;
  G4double currentRange;
  //
  G4double fr;
  G4double rangeinit;
  G4double geombig;
  G4double geomlimit;
  G4double tlimit;
  G4double tgeom;
  //
  G4double par1;
  G4double par2;
  G4double par3;
  G4double tlimitminfix2;
  G4double tausmall;
  G4double mass;
  G4double taulim;
  //
  //
  G4double presafety;
  G4double fZeff;
  //
  G4int    charge;
  G4int    currentMaterialIndex;
  //
  G4bool   firstStep;
  //
  G4LossTableManager*         theManager;
  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;
  const G4MaterialCutsCouple* currentCouple;

  G4GoudsmitSaundersonTable*  fGSTable;
  G4GSPWACorrections*         fPWACorrection;

  G4bool   fIsUsePWACorrection;
  G4bool   fIsUseMottCorrection;
  //
  G4double fLambda0; // elastic mean free path
  G4double fLambda1; // first transport mean free path
  G4double fScrA;    // screening parameter
  G4double fG1;      // first transport coef.
  // in case of Mott-correction
  G4double fMCtoScrA;
  G4double fMCtoQ1;
  G4double fMCtoG2PerG1;
  //
  G4double fTheTrueStepLenght;
  G4double fTheTransportDistance;
  G4double fTheZPathLenght;
  //
  G4ThreeVector fTheDisplacementVector;
  G4ThreeVector fTheNewDirection;
  //
  G4bool fIsEndedUpOnBoundary;  // step ended up on boundary i.e. transportation is the winer
  G4bool fIsMultipleSacettring;
  G4bool fIsSingleScattering;
  G4bool fIsEverythingWasDone;
  G4bool fIsNoScatteringInMSC;
  G4bool fIsNoDisplace;
  G4bool fIsInsideSkin;
  G4bool fIsWasOnBoundary;
  G4bool fIsFirstRealStep;
  //
  static G4bool gIsUseAccurate;
  static G4bool gIsOptimizationOn;
};

////////////////////////////////////////////////////////////////////////////////
inline
void G4GoudsmitSaundersonMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    charge = (G4int)(p->GetPDGCharge()/CLHEP::eplus);
    mass = p->GetPDGMass();
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline
G4double G4GoudsmitSaundersonMscModel::Randomizetlimit()
{
  G4double temptlimit;
    do {
         temptlimit = G4RandGauss::shoot(rndmEngineMod,tlimit,0.1*tlimit);
       } while ( (temptlimit<0.) || (temptlimit>2.*tlimit));

  return temptlimit;
}



#endif
