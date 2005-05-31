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
// $Id: G4VDNAProcessInWater.hh,v 1.1 2005-05-31 09:58:40 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VDNAProcessInWater_hh
 #define G4VDNAProcessInWater_hh 1
 
 #include "G4VLowEnergyTestableDiscreteProcess.hh"
 
 class G4VDNAProcessInWater : public G4VLowEnergyTestableDiscreteProcess
 {
  public:
                                         G4VDNAProcessInWater(const G4String & name) : G4VLowEnergyTestableDiscreteProcess(name) {}
   virtual                              ~G4VDNAProcessInWater() {}
 
   virtual G4VParticleChange *           PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);
   void                                  ValidateInWater(const G4Track & aTrack) const;

  private:
   // Hides default constructor and assignment operator as private 
                                         G4VDNAProcessInWater();
   G4VDNAProcessInWater &                operator=(const G4VDNAProcessInWater & right);
 };

#endif /* G4VDNAProcessInWater_hh */

