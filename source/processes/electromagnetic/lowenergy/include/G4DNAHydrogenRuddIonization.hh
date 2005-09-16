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
// $Id: G4DNAHydrogenRuddIonization.hh,v 1.1 2005-09-16 08:41:52 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAHYDROGENRUDDIONIZATION_HH
 #define  G4DNAHYDROGENRUDDIONIZATION_HH 1

 #include "G4DNAIonizationInWater.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNAHydrogenRuddIonizationFinalStatesPolicy.hh"
 #include "G4LogLogInterpolation.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

 class G4DNAHydrogenRuddIonizationEnergyLimitsPolicy
 {
  protected:
                      G4DNAHydrogenRuddIonizationEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAHydrogenRuddIonizationIncomingParticlePolicy
 {
  protected:
                                        G4DNAHydrogenRuddIonizationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };

 class G4DNAHydrogenRuddDataFilePolicy
 {
  public :
                                        G4DNAHydrogenRuddDataFilePolicy();
   const G4double                       lowEnergyLimit;
   const G4bool                         zeroBelowLowEnergyLimit;
   const G4double                       highEnergyLimit;
   const G4bool                         zeroAboveHighEnergyLimit;
   const G4double                       dataFileEnergyUnit;
   const G4double                       dataFileCrossSectionUnit;
   const char * const                   dataFileName;
 };
 
 class G4DNAHydrogenRuddIonization : public G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAHydrogenRuddIonizationIncomingParticlePolicy, G4DNAHydrogenRuddDataFilePolicy, G4LogLogInterpolation>, G4DNAHydrogenRuddIonizationFinalStatesPolicy<G4DNAHydrogenRuddIonizationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAHydrogenRuddIonization(const G4String & name = "G4DNAHydrogenRuddIonization") : G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAHydrogenRuddIonizationIncomingParticlePolicy, G4DNAHydrogenRuddDataFilePolicy, G4LogLogInterpolation>, G4DNAHydrogenRuddIonizationFinalStatesPolicy<G4DNAHydrogenRuddIonizationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAHydrogenRuddIonization() {}
 };
#endif /* G4DNAHYDROGENRUDDIONIZATION_HH */
