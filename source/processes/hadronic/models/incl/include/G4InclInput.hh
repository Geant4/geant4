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
#ifndef G4INCLINPUT_HH
#define G4INCLINPUT_HH 1

#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ParticleTable.hh"

#define FSIZE 15
/**
 * Initial values of a hadronic cascade problem.
 */
class G4InclInput {
public:
  G4InclInput() {
    isExtended = false;
    breakupThreshold = 10;
    fTargetA = 0;
    fTargetZ = 0;
    fBulletType = 0;
    fBulletE = 0.0;
    fTimeScale = 1.0;
    fNuclearPotential = 45.0; // Nuclear potential
    icoup = 0;
    
    theExtendedProjectileA = 0;
    theExtendedProjectileZ = 0;
    isExtended = false;

    fMinProtonE = 0.0;
    fNuclearPotential = 45.0;
    fTimeScale = 1.0;
    fMinNeutronEnergy = 0.0;

    usingInverseKinematics = false;
  };

  G4InclInput(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus, G4bool inverseKinematics);

  ~G4InclInput();
  
  void printInfo();

  static void printProjectileTargetInfo(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus);

  static G4bool canUseInverseKinematics(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus);

  G4double bulletE() {
    return fBulletE;
  }

  G4int getClusterOption() { return 0; }; // No clusters (and in 4.2 there never will be!)

  G4int bulletType() {
    return fBulletType;
  };

  void setExtendedProjectileInfo(const G4ParticleDefinition *pd);

  G4int getBulletType(const G4ParticleDefinition *pd);
  static G4ParticleDefinition* getParticleDefinition(G4int inclParticleCode);

  G4bool isInverseKinematics() {
    return usingInverseKinematics;
  };

  G4int targetA() { return fTargetA; };
  G4int targetZ() { return fTargetZ; };

  G4int extendedProjectileA() { return theExtendedProjectileA; };
  G4int extendedProjectileZ() { return theExtendedProjectileZ; };
  G4bool isExtendedProjectile() { return isExtended; };
  void isExtendedProjectile(G4bool ext) { isExtended = ext; };

  G4double getPotential() { return fNuclearPotential; };

  G4int getBreakupThreshold() { return breakupThreshold; };
  G4double getTimeScale() { return fTimeScale; };

private:
  G4int theExtendedProjectileA;
  G4int theExtendedProjectileZ;
  G4bool isExtended;

  G4int breakupThreshold;
  /**
   * Here f is an array containing the following initial values:
   * - f[0] : target mass number
   * - f[1] : target charge number
   */
  G4int fTargetA, fTargetZ;

  /*
   * - f[2] : bullet energy
   */
  G4double fBulletE;

  /*
   * - f[3] : minimum proton energy to leave the target (default: 0.0)
   */
  G4double fMinProtonE;

  /*
   * - f[4] : nuclear potential (default: 45.0 MeV)
   */
  G4double fNuclearPotential;

  /*
   * - f[5] : time scale (default: 1.0)
   */
  G4double fTimeScale;

  /*
   * - f[6] : bullet type (1: proton, 2: neutron, 3: pi+, 4: pi0 5: pi-, 6:H2, 7: H3, 8: He3, 9: He4
   */
  G4int fBulletType;

  /*
   * - f[7] : minimum neutron energy to leave the target (default: 0.0)
   */
  G4double fMinNeutronEnergy;

  /*
   * - f[8] : target material identifier (G4Mat)
   * - f[9] : not used
   * - f[10] : not used
   * - f[11] : not used
   * - f[12] : not used
   * - f[13] : not used
   * - f[14] : not used
   */
  //  G4double f[FSIZE];

  /**
   * Number of events to be processed.
   */
  G4int icoup;

  G4bool usingInverseKinematics;
};

#endif
