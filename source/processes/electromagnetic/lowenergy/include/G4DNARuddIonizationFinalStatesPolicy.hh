//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DNARuddIonizationTotalCrossSectionPolicy, 2005/09/13 10:17:45 francis
// GEANT4 tag $Name: emlowen-V07-01-06

#ifndef  G4DNARuddIonizationFinalStatesPolicy_HH
#define  G4DNARuddIonizationFinalStatesPolicy_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNARuddIonizationFinalStatesPolicy : public EnergyLimitsPolicy
 {
  protected:
                                        G4DNARuddIonizationFinalStatesPolicy() {}
                                       ~G4DNARuddIonizationFinalStatesPolicy() {}


  G4bool                                KillIncomingParticle(G4double energy) const;
  void                                  BuildFinalStatesData(void) const;
  G4double                              RandomizeEjectedElectronEnergy(G4double incomingParticleEnergy, G4int shell) const;
  void                                  RandomizeEjectedElectronDirection(G4double incomingParticleEnergy, G4double
                                        outgoingParticleEnergy, G4double & cosTheta, G4double & phi ) const;
  G4double                              EnergyConstant(G4int ionizationLevel) const;

  private:
  G4double                             DifferentialCrossSection(G4double k, G4double energyTransfer, G4int shell) const;
  G4double                             ShellBindingEnergy(G4int ionizationLevelIndex) const;

  // Hides default constructor and assignment operator as private
                                        G4DNARuddIonizationFinalStatesPolicy(const G4DNARuddIonizationFinalStatesPolicy & copy);
   G4DNARuddIonizationFinalStatesPolicy & operator=(const G4DNARuddIonizationFinalStatesPolicy & right);
 };

 #include "G4DNARuddIonizationFinalStatesPolicy.icc"
#endif /* G4DNARuddIonizationFinalStatesPolicy_HH */


