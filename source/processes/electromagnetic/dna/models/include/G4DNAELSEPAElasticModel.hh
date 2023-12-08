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
// Created on 2016/01/18
//
// Authors: D. Sakata, S. Incerti
//
// Based on a recent release of the ELSEPA code 
// developed and provided kindly by F. Salvat et al. 
// See
// Computer Physics Communications, 165(2), 157-190. (2005)
// http://dx.doi.org/10.1016/j.cpc.2004.09.006
//

#ifndef G4DNAELSEPAElasticModel_h
#define G4DNAELSEPAElasticModel_h 1

#include <map>
#include "G4DNACrossSectionDataSet.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4DNAELSEPAElasticModel : public G4VEmModel
{

public:

  G4DNAELSEPAElasticModel(const G4ParticleDefinition* particle = nullptr,
                          const G4String& nam = "DNAELSEPAElasticModel");

  ~G4DNAELSEPAElasticModel() override;

  G4DNAELSEPAElasticModel & operator=(const G4DNAELSEPAElasticModel &right) = delete;
  G4DNAELSEPAElasticModel(const G4DNAELSEPAElasticModel&) = delete;

  void Initialise(
                   const G4ParticleDefinition* particle, const G4DataVector&) override;

  G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* particle,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy) override;

  void SetMaximumEnergy (G4double input)
                   {fhighEnergyLimit = input; SetHighEnergyLimit(input);};
  
  void     SetKillBelowThreshold (G4double threshold);
  
  G4double GetKillBelowThreshold() {return fkillBelowEnergy_Au;}

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
 
  const std::vector<G4double>* fpMolDensity=nullptr;
  std::vector <G4double> kIntersectionEnergySR;

  G4double fkillBelowEnergy_Au=0.;
  G4double flowEnergyLimit=0.;
  G4double fhighEnergyLimit=0.;
  
  G4bool isInitialised=false;
  G4int verboseLevel=0;

  G4DNACrossSectionDataSet* fpData_Au=nullptr;
  G4DNACrossSectionDataSet* fpData_H2O=nullptr;

  G4double Theta(G4int Z, G4ParticleDefinition * aParticleDefinition,
                 G4double k,
                 G4double integrDiff);

  G4double LinLinInterpolate(G4double e1,
                             G4double e2,
                             G4double e,
                             G4double xs1,
                             G4double xs2);

  G4double LinLogInterpolate(G4double e1,
                             G4double e2,
                             G4double e,
                             G4double xs1,
                             G4double xs2);

  G4double LogLinInterpolate(G4double e1,
                             G4double e2,
                             G4double e,
                             G4double xs1,
                             G4double xs2);

  G4double LogLogInterpolate(G4double e1,
                             G4double e2,
                             G4double e,
                             G4double xs1,
                             G4double xs2);

  G4double QuadInterpolator(G4double e11,
                            G4double e12,
                            G4double e21,
                            G4double e22,
                            G4double x11,
                            G4double x12,
                            G4double x21,
                            G4double x22,
                            G4double t1,
                            G4double t2,
                            G4double t,
                            G4double e);

  G4double RandomizeCosTheta(G4int Z, G4double k);

  using TriDimensionMapZ = std::map<G4int, std::map<G4double, std::map<G4double, G4double>>>;
  TriDimensionMapZ fAngleDataZ;

  std::map <G4int, std::vector<G4double> > eEdummyVecZ;

  using VecMap = std::map<G4double, std::vector<G4double>>;
  VecMap eCum_Au;
  VecMap eCum_H2O;

  using TriDimensionMap = std::map<G4double, std::map<G4double, G4double>>;
  TriDimensionMap fAngleData_Au;
  TriDimensionMap fAngleData_H2O;

  std::vector<G4double> eEdummyVec_Au;
  std::vector<G4double> eEdummyVec_H2O;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
