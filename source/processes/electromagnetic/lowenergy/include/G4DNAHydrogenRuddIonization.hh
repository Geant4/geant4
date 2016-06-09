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
// $Id: G4DNAHydrogenRuddIonization.hh,v 1.3 2006/06/29 19:34:36 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#ifndef   G4DNAHYDROGENRUDDIONIZATION_HH
 #define  G4DNAHYDROGENRUDDIONIZATION_HH 1

 #include "G4DNAIonizationInWater.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNARuddIonizationFinalStatesPolicy.hh"
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
 
 class G4DNAHydrogenRuddIonization : public G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAHydrogenRuddIonizationIncomingParticlePolicy, G4DNAHydrogenRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAHydrogenRuddIonizationEnergyLimitsPolicy, G4DNAHydrogenRuddIonizationIncomingParticlePolicy> >
 {
  public:
                                         G4DNAHydrogenRuddIonization(const G4String & name = "G4DNAHydrogenRuddIonization") : G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAHydrogenRuddIonizationIncomingParticlePolicy, G4DNAHydrogenRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAHydrogenRuddIonizationEnergyLimitsPolicy, G4DNAHydrogenRuddIonizationIncomingParticlePolicy> > (name) {}
   virtual                              ~G4DNAHydrogenRuddIonization() {}
 };
#endif /* G4DNAHYDROGENRUDDIONIZATION_HH */
