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
// $Id: G4DNAElectronElasticBrenner.hh,v 1.1 2005-06-02 15:02:54 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Phys. Med. Biol. 29 N.4 (1983) 443-447

#ifndef G4DNAElectronElasticBrenner_hh
 #define G4DNAElectronElasticBrenner_hh 1
 
 #include "G4VDNAElectronElasticScatteringInWater.hh"
 
 class G4DNAElectronElasticBrenner : public G4VDNAElectronElasticScatteringInWater
 {
  public:
                                         G4DNAElectronElasticBrenner(const G4String & name="DNAElectronElasticBrenner");
   virtual                              ~G4DNAElectronElasticBrenner() {}

   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);
 
  protected:
   virtual G4double                      RandomizeCosTheta(G4double k, G4int z);
   virtual G4double                      TotalCrossSection (G4double k, G4int z);

  private:
   G4double                              CalculatePolynomial(G4double k, const G4double *vector, G4int size) const;

   // Hides default constructor and assignment operator as private 
   G4DNAElectronElasticBrenner &         operator=(const G4DNAElectronElasticBrenner & right);

   const G4double lowEnergyLimit;
   const G4double highEnergyLimit;
 };

#endif /* G4DNAElectronElasticBrenner_hh */

