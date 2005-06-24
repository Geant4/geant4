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
// $Id: G4DNAElectronElasticScatteringInWater.hh,v 1.3 2005-06-24 10:07:13 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Nucl. Instr. Meth. 155 (1978) 145-156
// J. Phys. D: Appl. Phys. 33 (2000) 932-944
// Phys. Med. Biol. 45 (2000) 3171-3194

#ifndef   G4DNAELECTRONELASTICSCATTERINGINWATER_HH
 #define  G4DNAELECTRONELASTICSCATTERINGINWATER_HH 1
 
 #include "G4VDNAProcessInWater.hh"
 
 template<typename TotalCrossSectionPolicy, typename FinalStatesPolicy>
 class G4DNAElectronElasticScatteringInWater : public G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>
 {
  public:
                                         G4DNAElectronElasticScatteringInWater(const G4String & name) : G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>(name) {}
   virtual                              ~G4DNAElectronElasticScatteringInWater() {}
 
   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);

  private:
   // Hides default constructor and assignment operator as private 
                                         G4DNAElectronElasticScatteringInWater(const G4DNAElectronElasticScatteringInWater & copy);
   G4DNAElectronElasticScatteringInWater & operator=(const G4DNAElectronElasticScatteringInWater & right);
 };
 
 #include "G4DNAElectronElasticScatteringInWater.icc"
#endif /* G4DNAELECTRONELASTICSCATTERINGINWATER_HH */

