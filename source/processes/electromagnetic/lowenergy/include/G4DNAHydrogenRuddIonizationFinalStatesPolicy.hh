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
// $Id: G4DNAHydrogenRuddIonizationTotalCrossSectionPolicy, 2005/09/16 9:25:45 Ziad FRANCIS
// GEANT4 tag $Name: emlowen-V07-01-09

#ifndef  G4DNAHydrogenRuddIonizationFinalStatesPolicy_HH
#define  G4DNAHydrogenRuddIonizationFinalStatesPolicy_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNAHydrogenRuddIonizationFinalStatesPolicy : public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAHydrogenRuddIonizationFinalStatesPolicy() {}
                                       ~G4DNAHydrogenRuddIonizationFinalStatesPolicy() {}


  G4bool                                KillIncomingParticle(G4double energy) const;
  void                                  BuildFinalStatesData(void) const;
  G4double                              RandomizeEjectedElectronEnergy(G4double incomingParticleEnergy, G4int shell) const;
  void                                  RandomizeEjectedElectronDirection(G4double incomingParticleEnergy, G4double
                                        outgoingParticleEnergy, G4double & cosTheta, G4double & phi ) const;
  G4double                              EnergyConstant(G4int ionizationLevel) const;

  private:
  G4double                             DifferentialCrossSection(G4double k, G4double energyTransfer, G4int shell) const;
 
  // Hides default constructor and assignment operator as private
                                        G4DNAHydrogenRuddIonizationFinalStatesPolicy(const G4DNAHydrogenRuddIonizationFinalStatesPolicy & copy);
   G4DNAHydrogenRuddIonizationFinalStatesPolicy & operator=(const G4DNAHydrogenRuddIonizationFinalStatesPolicy & right);
 };

 #include "G4DNAHydrogenRuddIonizationFinalStatesPolicy.icc"
#endif /* G4DNAHydrogenRuddIonizationFinalStatesPolicy_HH */


