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
// $Id: G4DNAElectronElasticScatteringInWater.hh,v 1.1 2005-05-31 09:58:40 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Nucl. Instr. Meth. 155 (1978) 145-156
// J. Phys. D: Appl. Phys. 33 (2000) 932-944
// Phys. Med. Biol. 45 (2000) 3171-3194
// Phys. Med. Biol. 29 N.4 (1983) 443-447

#ifndef G4DNAElectronElasticScatteringInWater_hh
 #define G4DNAElectronElasticScatteringInWater_hh 1
 
 #include "G4VDNAProcessInWater.hh"
 
 class G4DNAElectronElasticScatteringInWater : public G4VDNAProcessInWater
 {
  public:
                                         G4DNAElectronElasticScatteringInWater(const G4String & name="DNAElectronElasticScatteringInWater");
   virtual                              ~G4DNAElectronElasticScatteringInWater() {}
 
   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);
   virtual G4bool                        IsApplicable(const G4ParticleDefinition & aParticleDefinition);

  protected:
   virtual G4double                      GetMeanFreePath(const G4Track & aTrack, G4double previousStepSize, G4ForceCondition * condition);

  private:
   G4double                              RutherfordTotalCrossSection(G4double k, G4int z);
   G4double                              ScreeningFactor(G4double k, G4int z);
   G4double                              GenerateCosThetaElasticEmfietzoglou(G4double k, G4int z);
   G4double                              GenerateCosThetaElasticBrenner(G4double k);
   G4double                              CalulatePolynomial(G4double k, const G4double *vector, G4int size);

   // Hides default constructor and assignment operator as private 
   G4DNAElectronElasticScatteringInWater & operator=(const G4DNAElectronElasticScatteringInWater & right);



   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   const G4double intrinsicLowEnergyLimit;
   const G4double intrinsicHighEnergyLimit;
 };

#endif /* G4DNAElectronElasticScatteringInWater_hh */

