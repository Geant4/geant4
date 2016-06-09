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
// $Id: G4DNAElectronBornExcitation.hh,v 1.3 2006/06/29 19:33:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

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
