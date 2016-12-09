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
// $Id: $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4GenericBiasingPhysics_h
#define G4GenericBiasingPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4GenericBiasingPhysics : public G4VPhysicsConstructor
{
public:
  
  G4GenericBiasingPhysics(const G4String& name = "BiasingP");
  virtual ~G4GenericBiasingPhysics();

public:
  // ------------------------------
  // -- Biasing activation methods:
  // ------------------------------
  // -- Used to select particles and processes to be under biasing:
  // ---- Put under biasing all physics processes of given particleName:
  void PhysicsBias(const G4String& particleName);
  // ---- Put under biasing processes in processToBiasNames of given particleName:
  void PhysicsBias(const G4String& particleName, const std::vector< G4String >& processToBiasNames);
  // ---- Allow for non physics biasing for particle:
  void NonPhysicsBias(const G4String& particleName);
  // ---- Put under biasing all physics processes and allow for non physics biasing:
  void Bias(const G4String& particleName);
  // ---- Put under biasing processes in processToBiasNames of given particleName:
  void Bias(const G4String& particleName, const std::vector< G4String >& processToBiasNames);

  // -- Bias groups of particles:
  // --  - particles which have been setup by names with above methods are not affected
  // --  - particles can be specified by PDG range
  // --  - particles can be specified by the charged ou neutral nature
  // --     - particles specified by name and PDG range are unaffected
  // -- Add a PDG range for particle to bias, anti-particles are included by default:
  void    PhysicsBiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle = true );
  void NonPhysicsBiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle = true );
  void           BiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle = true );
  // -- Will bias all charged particles:
  void    PhysicsBiasAllCharged( G4bool includeShortLived = false );
  void NonPhysicsBiasAllCharged( G4bool includeShortLived = false );
  void           BiasAllCharged( G4bool includeShortLived = false );
  // -- Will bias all neutral particles:
  void    PhysicsBiasAllNeutral( G4bool includeShortLived = false );
  void NonPhysicsBiasAllNeutral( G4bool includeShortLived = false );
  void           BiasAllNeutral( G4bool includeShortLived = false );

  
  // -------------------------------------------------------------
  // -- Activation of parallel geometries used by generic biasing:
  // -------------------------------------------------------------
  // -- Each method can be called several times:
  // --    - on a same particle type :
  // --        myBiasingPhysics->AddParallelGeometry("neutron", "geometry1");
  // --        myBiasingPhysics->AddParallelGeometry("neutron", "geometry2");
  // --    - on a range of PDG particle:
  // --        myBiasingPhysics->AddParallelGeometry(PDG1, PDG2, "geometryXX");
  // --        myBiasingPhysics->AddParallelGeometry(PDG3, PDG4, vectorOfGeometries);
  // -- etc.
  void AddParallelGeometry( const G4String& particleName, const G4String&                parallelGeometryName  );
  void AddParallelGeometry( const G4String& particleName, const std::vector< G4String >& parallelGeometryNames );
  void AddParallelGeometry( G4int PDGlow, G4int PDGhigh,  const G4String&                parallelGeometryName , G4bool includeAntiParticle = true );
  void AddParallelGeometry( G4int PDGlow, G4int PDGhigh,  const std::vector< G4String >& parallelGeometryNames, G4bool includeAntiParticle = true );
  void AddParallelGeometryAllCharged(                     const G4String&                parallelGeometryName , G4bool includeShortLived = false );
  void AddParallelGeometryAllCharged(                     const std::vector< G4String >& parallelGeometryNames, G4bool includeShortLived = false );
  void AddParallelGeometryAllNeutral(                     const G4String&                parallelGeometryName , G4bool includeShortLived = false );
  void AddParallelGeometryAllNeutral(                     const std::vector< G4String >& parallelGeometryNames, G4bool includeShortLived = false );
  


  // -- Information about biased particles:
  void BeVerbose() { fVerbose = true; }
  
public:
  
  // This method is dummy for physics
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
  
private:
  
  // hide assignment operator
  G4GenericBiasingPhysics & operator=(const G4GenericBiasingPhysics &right);
  G4GenericBiasingPhysics(const G4GenericBiasingPhysics&);

  // -- Particles under biasing:
  std::vector< G4String >  fBiasedParticles;
  std::vector< G4bool >   fBiasAllProcesses;
  // -- Related biased processes:
  std::vector< std::vector< G4String > > fBiasedProcesses;
  // -- non physics biased particles:
  std::vector< G4String > fNonPhysBiasedParticles;

  // -- Group of particles under biasing:
  std::vector< G4int >    fPhysBiasByPDGRangeLow,    fPhysBiasByPDGRangeHigh;
  std::vector< G4int > fNonPhysBiasByPDGRangeLow, fNonPhysBiasByPDGRangeHigh;
  G4bool fPhysBiasAllCharged, fNonPhysBiasAllCharged;
  G4bool fPhysBiasAllChargedISL, fNonPhysBiasAllChargedISL;
  G4bool fPhysBiasAllNeutral,    fNonPhysBiasAllNeutral;
  G4bool fPhysBiasAllNeutralISL, fNonPhysBiasAllNeutralISL;

  
  // -- Particles associated with parallel geometries:
  std::vector< G4String >                       fParticlesWithParallelGeometries;
  std::map< G4String, std::vector< G4String > > fParallelGeometriesForParticle;
  std::vector< G4int >                          fPDGlowParallelGeometries, fPDGhighParallelGeometries;
  std::map< G4int,    std::vector< G4String > > fPDGrangeParallelGeometries;
  std::vector< G4String >                       fParallelGeometriesForCharged,    fParallelGeometriesForNeutral;
  std::vector< G4bool >                         fAllChargedParallelGeometriesISL, fAllNeutralParallelGeometriesISL;

  
  void AssociateParallelGeometries();

  
  // -- Report:
  G4bool fVerbose;


  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
