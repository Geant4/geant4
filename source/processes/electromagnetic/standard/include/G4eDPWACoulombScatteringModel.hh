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
// GEANT4 Class header file
//
//
// File name:     G4eDPWACoulombScatteringModel
//
// Author:        Mihaly Novak
//
// Creation date: 02.07.2020
//
// Modifications:
//
// Class Description:
//
// e-/e+ Coulomb scattering model based on numerical Differential Cross Sections
// (DCS) obtained by Dirac Partial Wave Analysis (DPWA) and supplied by the
// G4eDPWAElasticDCS class.
// The model contains the possibility to incorporate the effects of angular
// deflections of sub-threshold ionisation intercations when it's described by
// the condensed history model. Note, this must be inactivated (by setting the
// `isscpcor` input argument of the CTR to false) when ionisation is described
// with a classical, event by event based simulation model instead of usign the
// condensed history approach (otherwise, the corresponding angular defelctions
// will be "double counted").
//
// -------------------------------------------------------------------



#ifndef G4eDPWACoulombScatteringModel_h
#define G4eDPWACoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "globals.hh"

class G4eDPWAElasticDCS;
class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4DataVector;

class G4eDPWACoulombScatteringModel : public G4VEmModel {

public:

  /**
   * Constructor.
   *
   * @param[in] ismixed  Indicates if the model is for mixed or for pure single
   *                     Coulomb scattering. Different type of tables are pre-
   *                     pared for sampling polar angle of Coulomb scattering
   *                     for mixed and for pure single scattering models: cosine
   *                     of the polar scattering angle can be sampled in a
   *                     restriced inteval (see mumin input parameter below).
   * @param[in] isscpcor Indicates if scattering power correction should be used.
   *                     Note, scattering power correction accounts the effects
   *                     angular deflections due to sub-threshold ionisations
   *                     when ionisation is described by using condensed history
   *                     model (should be active only in this case).
   * @param[in] mumin    When the model is used for mixed simulation, Coulomb
   *                     scatterings, resulting in a minimum t_c polar angular
   *                     deflection, modelled explicitly. Therefore, cross
   *                     sections are computed, and angular deflections are
   *                     sampled ina resricted [\theta_c,\pi] interval. The
   *                     minimum of this interval is determined by the mumin
   *                     parameter as:
   *                     \mu_{min} = \mu(\theta_c)=0.5[1-\cos(\theta_c)]
   */
  G4eDPWACoulombScatteringModel(G4bool ismixed=false, G4bool isscpcor=true,
                                G4double mumin=0.0);

  ~G4eDPWACoulombScatteringModel() override;

  //
  // Interface methods:

  void     Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void     InitialiseLocal(const G4ParticleDefinition*, G4VEmModel*) override;

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*, G4double ekin,
                                      G4double Z, G4double A, G4double prodcut,
                                      G4double emax) override;

  void     SampleSecondaries(std::vector<G4DynamicParticle*>*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                             G4double tmin,
                             G4double maxEnergy) override;

  G4double MinPrimaryEnergy(const G4Material*, const G4ParticleDefinition*, 
                            G4double) override { return 10.0*CLHEP::eV; }

  void     SetTheDCS(G4eDPWAElasticDCS* theDCS) { fTheDCS = theDCS; }

  G4eDPWAElasticDCS* GetTheDCS() { return fTheDCS; }


private:

  // Indicates if the model is mixed: MSC for soft (theta<theta_c), Singe
  // Scattering(SS) for hard scatterings(theta>theta_c). SS otherwise.
  // Note, that while the model provides restricted (elastic and transport)
  // cross sections, it's responsible to handle, i.e. provide final state,
  // only for the Singe Scattering part in case of a mixed model.
  G4bool                     fIsMixedModel;
  // indicates if scattering power correction should be applied: correction due
  // to deflection in case of sub-threshold, inelastic interactions -> only in
  // case of condensed history simulation of inonisation!
  G4bool                     fIsScpCorrection;
  // mu(theta)=0.5[1-cos(theta)]: the model porvides final states \in [fMuMin,1]
  G4double                   fMuMin;
  // the object that provides cross sections and polar angle of scattering
  G4eDPWAElasticDCS*         fTheDCS;
  // particle change
  G4ParticleChangeForGamma*  fParticleChange;

};

#endif
