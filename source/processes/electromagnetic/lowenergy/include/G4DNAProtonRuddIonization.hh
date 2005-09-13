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
// $Id: G4DNAProtonRuddIonization.hh,v 1.1 2005-09-13 08:59:06 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAPROTONRUDDIONIZATION_HH
 #define  G4DNAPROTONRUDDIONIZATION_HH 1

 #include "G4DNAIonizationInWater.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNARuddIonizationFinalStatesPolicy.hh"
 #include "G4LogLogInterpolation.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 #include "G4DNARuddIonizationTotalCrossSectionPolicy.hh"

 class G4DNAProtonRuddIonizationEnergyLimitsPolicy
 {
  protected:
                      G4DNAProtonRuddIonizationEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAProtonRuddIonizationIncomingParticlePolicy
 {
  protected:
                                        G4DNAProtonRuddIonizationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };
 
 class G4DNAProtonRuddDataFilePolicy
 {
  public :
                                        G4DNAProtonRuddDataFilePolicy();
   const G4double                       lowEnergyLimit;
   const G4bool                         zeroBelowLowEnergyLimit;
   const G4double                       highEnergyLimit;
   const G4bool                         zeroAboveHighEnergyLimit;
   const G4double                       dataFileEnergyUnit;
   const G4double                       dataFileCrossSectionUnit;
   const char * const                   dataFileName;
 };
 
 class G4DNAProtonRuddIonization : public G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonRuddIonizationIncomingParticlePolicy, G4DNAProtonRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAProtonRuddIonizationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAProtonRuddIonization(const G4String & name = "G4DNAProtonRuddIonization") : G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonRuddIonizationIncomingParticlePolicy, G4DNAProtonRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAProtonRuddIonizationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAProtonRuddIonization() {}
 };
#endif /* G4DNAPROTONRUDDIONIZATION_HH */
