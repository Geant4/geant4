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
// $Id: G4DNAHeliumExcitation.hh,v 1.2 2006/06/29 19:34:24 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#ifndef   G4DNAHELIUMEXCITATION_HH
 #define  G4DNAHELIUMEXCITATION_HH 1
 
 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 
 class G4DNAHeliumExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAHeliumExcitationEnergyLimitsPolicy();
  
   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };
 
 class G4DNAHeliumExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAHeliumExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;

   G4double                             slaterEffectiveCharge[3];
   G4double                             sCoefficient[3];
   const G4double                       kineticEnergyCorrection;
 };
 
 class G4DNAHeliumExcitation : public G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAHeliumExcitationIncomingParticlePolicy, G4DNAHeliumExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAHeliumExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAHeliumExcitation(const G4String & name = "G4DNAHeliumExcitation") : G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAHeliumExcitationIncomingParticlePolicy, G4DNAHeliumExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAHeliumExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAHeliumExcitation() {}
 };
#endif /* G4DNAHELIUMEXCITATION_HH */
