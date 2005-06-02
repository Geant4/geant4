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
// $Id: G4VDNAElectronElasticScatteringInWater.hh,v 1.1 2005-06-02 15:02:54 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Nucl. Instr. Meth. 155 (1978) 145-156
// J. Phys. D: Appl. Phys. 33 (2000) 932-944
// Phys. Med. Biol. 45 (2000) 3171-3194

#ifndef G4VDNAElectronElasticScatteringInWater_hh
 #define G4VDNAElectronElasticScatteringInWater_hh 1
 
 #include "G4VDNAProcessInWater.hh"
 
 class G4VDNAElectronElasticScatteringInWater : public G4VDNAProcessInWater
 {
  public:
                                         G4VDNAElectronElasticScatteringInWater(const G4String & name);
   virtual                              ~G4VDNAElectronElasticScatteringInWater() {}
 
   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);
   virtual G4bool                        IsApplicable(const G4ParticleDefinition & aParticleDefinition);

  protected:
   virtual G4double                      RandomizeCosTheta(G4double k, G4int z) = 0 ;

   G4double                              RutherfordTotalCrossSection(G4double k, G4int z) const;
   G4double                              ScreeningFactor(G4double k, G4int z) const;

  private:
   // Hides default constructor and assignment operator as private 
                                         G4VDNAElectronElasticScatteringInWater();
   G4VDNAElectronElasticScatteringInWater & operator=(const G4VDNAElectronElasticScatteringInWater & right);
 };

#endif /* G4VDNAElectronElasticScatteringInWater_hh */

