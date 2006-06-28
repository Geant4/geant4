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
// $Id: G4DNAElectronBornExcitation.hh,v 1.2 2006-06-28 13:58:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAElectronBornExcitation_HH
 #define  G4DNAElectronBornExcitation_HH 1

 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 #include "G4DNATotalCrossSectionFromFilePolicy.hh"
 #include "G4DNABornExcitationFinalStatesPolicy.hh"
 #include "G4LogLogInterpolation.hh"

 class G4DNAElectronBornExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAElectronBornExcitationEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAElectronBornExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAElectronBornExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };
 
 class G4DNAElectronBornExcitationDataFilePolicy
 {
 public : G4DNAElectronBornExcitationDataFilePolicy();
 
 const G4double lowEnergyLimit;
 const G4bool   zeroBelowLowEnergyLimit;
 const G4double highEnergyLimit;
 const G4bool   zeroAboveHighEnergyLimit;
 const G4double dataFileEnergyUnit;
 const G4double dataFileCrossSectionUnit;
 const char * const dataFileName;
 };

 class G4DNAElectronBornExcitation : public G4DNAExcitationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAElectronBornExcitationIncomingParticlePolicy, G4DNAElectronBornExcitationDataFilePolicy, G4LogLogInterpolation>, G4DNABornExcitationFinalStatesPolicy<G4DNAElectronBornExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAElectronBornExcitation(const G4String & name = "G4DNAElectronBornExcitation") : G4DNAExcitationInWater<G4DNATotalCrossSectionFromFilePolicy<G4DNAElectronBornExcitationIncomingParticlePolicy, G4DNAElectronBornExcitationDataFilePolicy, G4LogLogInterpolation>, G4DNABornExcitationFinalStatesPolicy<G4DNAElectronBornExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAElectronBornExcitation() {}
 };
#endif /* G4DNAElectronBornExcitation_HH */
