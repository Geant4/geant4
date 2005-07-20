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
// $Id: G4DNAExcitationInWater.hh,v 1.1 2005-07-20 10:01:54 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef   G4DNAEXCITATIONINWATER_HH
 #define  G4DNAEXCITATIONINWATER_HH 1
 
 #include "G4VDNAProcessInWater.hh"
 
 // TotalCrossSectionPolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void)
 //  - [protected] G4double TotalCrossSection(G4double k, G4int z)
 //  - [protected] void BuildTotalCrossSection(void)
 //  - [protected] G4int RandomizePartialCrossSection(G4double k, G4int z)
 //  - [protected] G4double EnergyConstant(G4int excitationLevelIndex)
 
 // FinalStatesPolicy must provide:
 //  - [protected] G4bool KillIncomingParticle(G4double k)
 //  - [protected] void BuildFinalStatesData(void)
 
 template<typename TotalCrossSectionPolicy, typename FinalStatesPolicy>
 class G4DNAExcitationInWater : public G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>
 {
  public:
                                         G4DNAExcitationInWater(const G4String & name) : G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>(name) {}
   virtual                              ~G4DNAExcitationInWater() {}
 
   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);

  private:
   // Hides default constructor and assignment operator as private 
                                         G4DNAExcitationInWater(const G4DNAExcitationInWater & copy);
   G4DNAExcitationInWater &              operator=(const G4DNAExcitationInWater & right);
 };
 
 #include "G4DNAExcitationInWater.icc"
#endif /* G4DNAEXCITATIONINWATER_HH */

