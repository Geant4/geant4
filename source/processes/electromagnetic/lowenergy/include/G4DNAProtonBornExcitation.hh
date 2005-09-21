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
// $Id: G4DNAProtonBornExcitation.hh,v 1.1 2005-09-21 09:18:35 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAProtonBornExcitation_HH
 #define  G4DNAProtonBornExcitation_HH 1

 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNABornExcitationFinalStatesPolicy.hh"
 #include "G4LogLogInterpolation.hh"

 class G4DNAProtonBornExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAProtonBornExcitationEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAProtonBornExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAProtonBornExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };

 class G4DNAProtonBornExcitationDataFilePolicy
 {
  public :
                                        G4DNAProtonBornExcitationDataFilePolicy();
   const G4double                       lowEnergyLimit;
   const G4bool                         zeroBelowLowEnergyLimit;
   const G4double                       highEnergyLimit;
   const G4bool                         zeroAboveHighEnergyLimit;
   const G4double                       dataFileEnergyUnit;
   const G4double                       dataFileCrossSectionUnit;
   const char * const                   dataFileName;
 };


 class G4DNAProtonBornExcitation : public G4DNAExcitationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonBornExcitationIncomingParticlePolicy, G4DNAProtonBornExcitationDataFilePolicy, G4LogLogInterpolation>, G4DNABornExcitationFinalStatesPolicy<G4DNAProtonBornExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAProtonBornExcitation(const G4String & name = "G4DNAProtonBornExcitation") : G4DNAExcitationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAProtonBornExcitationIncomingParticlePolicy, G4DNAProtonBornExcitationDataFilePolicy, G4LogLogInterpolation>, G4DNABornExcitationFinalStatesPolicy<G4DNAProtonBornExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAProtonBornExcitation() {}
 };
#endif /* G4DNAProtonBornExcitation_HH */
