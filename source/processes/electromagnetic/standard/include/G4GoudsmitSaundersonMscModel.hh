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
// $Id: G4GoudsmitSaundersonMscModel.hh 94933 2015-12-18 09:22:52Z gcosmo $
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
//
// Class description:
//   Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened
//   Rutherford DCS for elastic scattering of electrons/positrons. Step limitation
//   algorithm as well as true to geomerty and geometry to true step length
//   computations are adopted from Urban model[5].
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
class G4PWATotalXsecTable;


class G4GoudsmitSaundersonMscModel : public G4VMscModel
{
public:

  G4GoudsmitSaundersonMscModel(const G4String& nam = "GoudsmitSaunderson");

  virtual ~G4GoudsmitSaundersonMscModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  void StartTracking(G4Track*);

  G4double GetTransportMeanFreePath(const G4ParticleDefinition*, G4double);
  void SingleScattering(G4double &cost, G4double &sint);
  void SampleMSC();

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&,
					  G4double safety);

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4double& currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

  void SetOptionPWAScreening(G4bool opt){fIsUsePWATotalXsecData=opt;}

private:
  inline void SetParticle(const G4ParticleDefinition* p);

  inline G4double GetLambda(G4double);

  //  hide assignment operator
  G4GoudsmitSaundersonMscModel & operator=(const  G4GoudsmitSaundersonMscModel &right);
  G4GoudsmitSaundersonMscModel(const  G4GoudsmitSaundersonMscModel&);
  G4double GetTransportMeanFreePathOnly(const G4ParticleDefinition*,G4double);

  inline G4double Randomizetlimit();

private:
  CLHEP::HepRandomEngine*     rndmEngineMod;


  G4double lowKEnergy;
  G4double highKEnergy;
  G4double currentKinEnergy;
  G4double currentRange;

  G4double fr,rangeinit,geombig,geomlimit;
  G4double lambdalimit,tlimit,tgeom;
  G4int    charge,currentMaterialIndex;

  G4bool   firstStep;

  G4double par1,par2,par3,tlimitminfix2,tausmall,mass,taulim;

  G4LossTableManager*         theManager;
  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;
  const G4MaterialCutsCouple* currentCouple;

  static G4GoudsmitSaundersonTable* fgGSTable;
  static G4PWATotalXsecTable*       fgPWAXsecTable;

  G4bool fIsUsePWATotalXsecData;

  G4double presafety;
  G4double fZeff;

  //
  G4double fLambda0; // elastic mean free path
  G4double fLambda1; // first transport mean free path
  G4double fScrA;    // screening parameter
  G4double fG1;      // first transport coef.
  //
  G4double      fTheTrueStepLenght;
  G4double      fTheTransportDistance;
  G4double      fTheZPathLenght;
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
  static G4bool fgIsUseAccurate;
  static G4bool fgIsOptimizationOn;
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
  G4double temptlimit = tlimit;
    do {
         temptlimit = G4RandGauss::shoot(rndmEngineMod,tlimit,0.3*tlimit);
       } while ( (temptlimit<0.) || (temptlimit > 2.*tlimit));

  return temptlimit;
}



#endif
