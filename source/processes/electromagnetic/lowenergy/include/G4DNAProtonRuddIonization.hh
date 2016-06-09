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
// $Id: G4DNAProtonRuddIonization.hh,v 1.4 2006/06/29 19:35:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#ifndef   G4DNAPROTONRUDDIONIZATION_HH
 #define  G4DNAPROTONRUDDIONIZATION_HH 1

 #include "G4DNAIonizationInWater.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNARuddIonizationFinalStatesPolicy.hh"
 #include "G4LogLogInterpolation.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

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
 
 class G4DNAProtonRuddIonization : public G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonRuddIonizationIncomingParticlePolicy, G4DNAProtonRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAProtonRuddIonizationEnergyLimitsPolicy, G4DNAProtonRuddIonizationIncomingParticlePolicy> >
 {
  public:
                                         G4DNAProtonRuddIonization(const G4String & name = "G4DNAProtonRuddIonization") : G4DNAIonizationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonRuddIonizationIncomingParticlePolicy, G4DNAProtonRuddDataFilePolicy, G4LogLogInterpolation>, G4DNARuddIonizationFinalStatesPolicy<G4DNAProtonRuddIonizationEnergyLimitsPolicy, G4DNAProtonRuddIonizationIncomingParticlePolicy> > (name) {}
   virtual                              ~G4DNAProtonRuddIonization() {}
 };
#endif /* G4DNAPROTONRUDDIONIZATION_HH */
