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
// G4AdjointSimMessenger
//
// Class description:
//
// This class represents the Messenger that defines the G4UI macro commands
// allowing the user controlling an adjoint/reverse MC simulation. It calls
// methods of G4AdjointSimManager.
//
// 1) Start an adjoint simulation
// ---------------------------------------------------
// Command:
//  -/adjoint/start_run nb: Start an adjoint simulation with a number of
//                          events given by nb.
// 2) Definition of the external source
// ---------------------------------------------------
// The external source represents the real external source of particles till
// which adjoint particles are tracked in the reverse tracking mode of the
// simulation (see G4AdjointSimManager.hh and G4Application Developer guide for
// more infos). The user can define the source as the external surface of a
// sphere or of G4 volume of the geometry. He can also set the maximum energy
// of the source. If an adjoint particle get an energy higher than this maximum
// energy before reaching the external surface source it is killed without
// being registered.
// Commands:
//  -/adjoint/DefineSphericalExtSource R X Y Z unit_length
//   The external source is set on a sphere with radius R and centered on
//   position (X,Y,Z)
//  -/adjoint/DefineSphericalExtSourceCenteredOnAVolume pvol_name R unit_length
//   The external source is set on a sphere with radius R and with its center
//   position located at the center of the physical volume specified by the
//   name pvol_name.
//  -/adjoint/DefineExtSourceOnExtSurfaceOfAVolume pvol_name
//   The external surface is set as the external boundary of a the physical
//   volume with name pvol_name.
//  -/adjoint/SetExtSourceEmax Emax energy_unit
//   Set the maximum energy of the external source.
//
// 3) Definition of the adjoint source
// ---------------------------------------------------
// The adjoint source represents the source from which adjoint primary
// particles are generated (see G4AdjointSimManager.hh and G4Application
// Developer guide for more infos).
// The user can define the source as the external surface of a sphere or of
// G4 volume of the geometry. He sets the minimum maximum energy of the
// source and defines which type of adjoint primary particles should be
// considered.
// Commands:
//  -/adjoint/DefineSphericalAdjSource R X Y Z unit_length
//  The adjoint source is set on a sphere with radius R and centered on
//  position (X,Y,Z)
// -/adjoint/DefineSphericalAdjSourceCenteredOnAVolume pvol_name R unit_length
//  The external source is set on a sphere with radius R and with its center
//  position located at the center of the physical volume specified by the
//  name pvol_name.
// -/adjoint/DefineAdjSourceOnExtSurfaceOfAVolume pvol_name
//  The external surface is set as the external boundary of a the
//  physical volume with name pvol_name
// -/adjoint/SetAdjSourceEmin Emin energy_unit
//  Set the minimum energy of the external source
// -/adjoint/SetAdjSourceEmax Emax energy_unit
//   Set the maximum energy of the external source
// -/adjoint/ConsiderAsPrimary particle_name
//  The type of particle specified by  "particle_name" will be added in
//  the list of primary adjoint particles. The list of candidates depends on the
//  reverse physics processes considered in the simulation. At the most the 					potential
//  candidates are (e-, gamma, proton, ion).
// -/adjoint/NeglectAsPrimary particle_name
//  The type  of particle specified by  "particle_name" will be removed
// from the list of primary adjoint particles. The list of candidates depends
// on the reverse physics processes considered in the simulation. At the most
// the potential candidates are (e-, gamma, proton, ion).

// --------------------------------------------------------------------
//   Class Name:   G4AdjointSimMessenger
//   Author:       L. Desorgher, 2007-2009
//   Organisation: SpaceIT GmbH
//   Contract:     ESA contract 21435/08/NL/AT
//   Customer:     ESA/ESTEC
// --------------------------------------------------------------------
#ifndef G4AdjointSimMessenger_hh
#define G4AdjointSimMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4AdjointSimManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;

// --------------------------------------------------------------------

class G4AdjointSimMessenger : public G4UImessenger
{
  public:
    G4AdjointSimMessenger(G4AdjointSimManager*);
    ~G4AdjointSimMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    G4AdjointSimManager* theAdjointRunManager;

    G4UIdirectory* AdjointSimDir = nullptr;
    G4UIcommand* beamOnCmd = nullptr;

    G4UIcommand* DefineSpherExtSourceCmd = nullptr;
    G4UIcommand* DefineSpherExtSourceCenteredOnAVolumeCmd = nullptr;
    G4UIcmdWithAString* DefineExtSourceOnAVolumeExtSurfaceCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* setExtSourceEMaxCmd = nullptr;

    G4UIcommand* DefineSpherAdjSourceCmd = nullptr;
    G4UIcommand* DefineSpherAdjSourceCenteredOnAVolumeCmd = nullptr;
    G4UIcmdWithAString* DefineAdjSourceOnAVolumeExtSurfaceCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* setAdjSourceEminCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* setAdjSourceEmaxCmd = nullptr;

    G4UIcmdWithAString* ConsiderParticleAsPrimaryCmd = nullptr;
    G4UIcmdWithAString* NeglectParticleAsPrimaryCmd = nullptr;

    G4UIcmdWithAnInteger* setNbOfPrimaryFwdGammasPerEventCmd = nullptr;
    G4UIcmdWithAnInteger* setNbOfPrimaryAdjGammasPerEventCmd = nullptr;
    G4UIcmdWithAnInteger* setNbOfPrimaryAdjElectronsPerEventCmd = nullptr;
};

#endif
