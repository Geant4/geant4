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
// $Id: G4AdjointSimMessenger.hh 98735 2016-08-09 10:54:06Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class Name:		G4AdjointSimMessenger.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	-1st January 2007 creation by L. Desorgher
//		-November-December 2009 Some cleaning and adaptation for the first Release in the Geant4 toolkit, L. Desorgher   
//				
//
//-------------------------------------------------------------
//	Documentation:
//		This class represents the Messenger that defined the G4UI macro comands allowing the
//		user contreol an adjoint/reverse MC simulation. It calls methods of  G4AdjointSimManager
//	List of commands
//	-----------------
//		1)Start an adjoint simulation
//		--------------------------------------------
//			Command:
//			-/adjoint/start_run nb:	Start an adjoint simulation with a number of events given by nb.
//		2)Definition of the external source
//		---------------------------------------------------
//		The external source represents the real external source of particles till which adjoint particles are tracked in the reverse tracking mode
//		of the simulation (see  G4AdjointSimManager.hh and G4Application Developer guide for more infos).
//		The user can define the source as the external surface of a sphere or of G4 volume of the geometry. He can also set the maximum energy of the
//		source. If an adjoint particle get an energy higher than this maximum energy before reaching the external surface source it is killed without being registered.
//			Commands:
//			-/adjoint/DefineSphericalExtSource R X Y Z unit_length:
//					The external source is set on a sphere with radius R and centered on position (X,Y,Z) 
//				 
//			-/adjoint/DefineSphericalExtSourceCenteredOnAVolume phys_vol_name R unit_length
//					The external source is set on a sphere with radius R and with its center position located at the center of the 
//					the physical volume specified by the name phys_vol_name.
//			-/adjoint/DefineExtSourceOnExtSurfaceOfAVolume phys_vol_name 
//					The external surface is set as the external boundary of a the physical volume with name phys_vol_name
//			-/adjoint/SetExtSourceEmax  Emax energy_unit 
//					Set the maximum  energy of the external source
//
//
//		3)Definition of the adjoint source
//		---------------------------------------------------
//		The adjoint source represents the source from which adjoint primary particles are generated.(see  G4AdjointSimManager.hh and G4Application Developer guide for more infos)
//		The user can define the source as the external surface of a sphere or of G4 volume of the geometry. He set the minimum maximum energy of the
//		source and define which type of adjoint primary particles should be considered. 
//			Commands:
//			-/adjoint/DefineSphericalAdjSource R X Y Z unit_length:
//					The adjoint source is set on a sphere with radius R and centered on position (X,Y,Z) 
//				 
//			-/adjoint/DefineSphericalAdjSourceCenteredOnAVolume phys_vol_name R unit_length
//					The external source is set on a sphere with radius R and with its center position located at the center of the 
//					the physical volume specified by the name phys_vol_name.
//			-/adjoint/DefineAdjSourceOnExtSurfaceOfAVolume phys_vol_name 
//					The external surface is set as the external boundary of a the physical volume with name phys_vol_name
//		 	
//			-/adjoint/SetAdjSourceEmin  Emin energy_unit 
//					Set the minimum  energy of the external source
//		
//			-/adjoint/SetAdjSourceEmax  Emax energy_unit 
//					Set the maximum  energy of the external source
//			
//			-/adjoint/ConsiderAsPrimary  particle_name 
//					The type  of particle specified by  "particle_name" will be added in the list of primary adjoint particles. 
//					The list of candidates depends on the reverse physics processes considered in the simulation. At the most the 
//					potential candidates are (e-, gamma, proton , ion) 						 
//			
//			-/adjoint/NeglectAsPrimary  particle_name 
//					The type  of particle specified by  "particle_name" will be removed from the list of primary adjoint particles. 
//					The list of candidates depends on the reverse physics processes considered in the simulation. At the most the 
//					potential candidates are (e-, gamma, proton , ion) 
//			
//

#ifndef G4AdjointSimMessenger_h
#define G4AdjointSimMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4AdjointSimManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
/*
#ifdef G4MULTITHREADED
class G4MTAdjointSimManager;
#endif
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4AdjointSimMessenger: public G4UImessenger
{
  public:
    G4AdjointSimMessenger(G4AdjointSimManager* );
/*
#ifdef G4MULTITHREADED
    G4AdjointSimMessenger(G4MTAdjointSimManager* );
#endif
*/

   ~G4AdjointSimMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    G4AdjointSimManager* theAdjointRunManager;
/*
#ifdef G4MULTITHREADED
    G4MTAdjointSimManager* theMTAdjointRunManager;
#endif
*/
    
    G4UIdirectory*             AdjointSimDir;
    G4UIcommand *               beamOnCmd;
    
    G4UIcommand *  DefineSpherExtSourceCmd;
    G4UIcommand *  DefineSpherExtSourceCenteredOnAVolumeCmd;
    G4UIcmdWithAString *  DefineExtSourceOnAVolumeExtSurfaceCmd;
    G4UIcmdWithADoubleAndUnit*  setExtSourceEMaxCmd;
    
    G4UIcommand *  DefineSpherAdjSourceCmd;
    G4UIcommand *  DefineSpherAdjSourceCenteredOnAVolumeCmd;
    G4UIcmdWithAString *  DefineAdjSourceOnAVolumeExtSurfaceCmd;
    
    G4UIcmdWithADoubleAndUnit*  setAdjSourceEminCmd;
    G4UIcmdWithADoubleAndUnit*  setAdjSourceEmaxCmd;

     
    G4UIcmdWithAString*  ConsiderParticleAsPrimaryCmd;
    G4UIcmdWithAString*  NeglectParticleAsPrimaryCmd;

    G4UIcmdWithAnInteger*  setNbOfPrimaryFwdGammasPerEventCmd;
    G4UIcmdWithAnInteger*  setNbOfPrimaryAdjGammasPerEventCmd;
    G4UIcmdWithAnInteger*  setNbOfPrimaryAdjElectronsPerEventCmd;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

