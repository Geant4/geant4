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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PrimaryGeneratorAction header
// --------------------------------------------------------------

#ifndef DMXPrimaryGeneratorMessenger_h
#define DMXPrimaryGeneratorMessenger_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#ifdef DMXENV_GPS_USE
class G4GeneralParticleSource;
#else
class DMXParticleSource;
#endif
class DMXAnalysisManager;

class G4Event;

class DMXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

   public:
      DMXPrimaryGeneratorAction();
      ~DMXPrimaryGeneratorAction();
  
   public:
      void GeneratePrimaries(G4Event* anEvent);

   private:
#ifdef DMXENV_GPS_USE
  G4GeneralParticleSource* particleGun;
#else
  DMXParticleSource* particleGun;
#endif

   private:
  //      const long* seeds;
  long seeds[2];
  G4double energy_pri;

   public:
  const long* GetEventSeeds()        {return seeds;};
  G4double GetEnergyPrimary()   {return energy_pri;};


};

#endif

