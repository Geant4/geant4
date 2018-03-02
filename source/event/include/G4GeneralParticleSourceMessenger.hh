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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSourceMessenger.hh
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    After changes to version 1.1 as in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - old commonds have been retained for backward compatibility, will be
//       removed in the future. 
//
// Version 3.0, Aug-Oct 2014, Andrea Dotti
//    Transformations for thread safety and use in MT application
//    Messenger is now a singleton w/ explicit Destroy() method for deletion
//    Note the following: the class should be instantiated only once
//    by a worker thread. It relies on a new feature of basic messenger class
//    that allows for UI commands to be created by worker threads but being
//    executed by master thread. For this reason the messenger itself should
//    be created once, form here the singleton pattern
///////////////////////////////////////////////////////////////////////////////
//
//
// Class Description:
//
// The function of the G4GeneralParticleSourceMessenger is to allow the user to
// enter commands either in interactive command line mode or through macros to
// control the G4GeneralParticleSource. 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4GeneralParticleSourceMessenger(G4GeneralParticleSource *fPtclGun)
//     Constructor:  Sets up commands.
//
// ~G4GeneralParticleSourceMessenger()
//     Destructor:  Deletes commands.
//
// void SetParticleGun(G4SingleParticleSource *fpg) { fParticleGun = fpg; } ;
//     To selecte the particle gun to be defined/modified. 
// void SetNewValue(G4UIcommand *command, G4String newValues)
//     Uses the appropriate methods in the G4GeneralParticleSource to carry out
//     the user commands.
// G4String GetCurrentValue(G4UIcommand *command)
//     Allows the user to retrieve the current values of parameters.
//     Not implemented yet.
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef G4GeneralParticleSourceMessenger_h
#define G4GeneralParticleSourceMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class G4SingleParticleSource;
class G4GeneralParticleSource;

/** Andrea Dotti Feb 2015
 * GPS messenger design requires some explanation for what distributions
 * parameters are concerned : Each thread has its own GPS
 * since primary generation is a user action.
 * However to save memory the underlying structures that provide the
 * GPS functionalities ( the G4SPS*Distribution classes and the
 * G4SPSRandomGenerator class)
 * are shared among threads. This implies that modifying parameters of sources
 * requires some attention:
 * 1- Only one thread should change source parameters.
 * 2- Changing of parameters can happen only between runs, when is guaranteed
 *    that no thread is accessing them
 * 2- UI commands require that even if messenger is instantiated in a thread
 *    the commands are executed in the master (this is possible since V10.1)
 * The simplest solution is to use UI commands to change GPS parameters and
 * avoid C++ APIs. If this is inevitable a simple solution is to instantiate
 * an instance of G4GeneralParticleSource explicitly in the master thread
 * (for example in G4VUserActionInitialization::BuildForMaster() and set the
 * defaults parameter there).
 */

class G4GeneralParticleSourceMessenger: public G4UImessenger
{
  public:

    void SetParticleGun(G4SingleParticleSource *fpg) { fParticleGun = fpg; } ;
    // Select the particle gun to be defined/modified
   
    void SetNewValue(G4UIcommand *command, G4String newValues);
    // Identifies the command which has been invoked by the user, extracts the
    // parameters associated with that command (held in newValues), and uses
    // these values with the appropriate member function of G4GeneralParticleSource.

    G4String GetCurrentValue(G4UIcommand *command);

    static G4GeneralParticleSourceMessenger* GetInstance(G4GeneralParticleSource*);
    static void Destroy();
 private:
    G4GeneralParticleSourceMessenger(G4GeneralParticleSource*);
    ~G4GeneralParticleSourceMessenger();

    void IonCommand(G4String newValues);
    void IonLvlCommand(G4String newValues);

  private:
    G4GeneralParticleSource* fGPS;
    G4SingleParticleSource* fParticleGun;
    G4ParticleTable* particleTable;
    G4String histtype;
    
  private: //commands
    G4UIdirectory* gpsDirectory;

    // multiple source control commands
    G4UIdirectory              *sourceDirectory;
    G4UIcmdWithADouble         *addsourceCmd;
    G4UIcmdWithoutParameter    *listsourceCmd;
    G4UIcmdWithoutParameter    *clearsourceCmd;
    G4UIcmdWithoutParameter    *getsourceCmd;
    G4UIcmdWithAnInteger       *setsourceCmd;  
    G4UIcmdWithADouble         *setintensityCmd;
    G4UIcmdWithAnInteger       *deletesourceCmd;
    G4UIcmdWithABool           *multiplevertexCmd;
    G4UIcmdWithABool           *flatsamplingCmd;

    // positional commands
    G4UIdirectory              *positionDirectory;
    G4UIcmdWithAString         *typeCmd1;
    G4UIcmdWithAString         *shapeCmd1;
    G4UIcmdWith3VectorAndUnit  *centreCmd1;
    G4UIcmdWith3Vector         *posrot1Cmd1;
    G4UIcmdWith3Vector         *posrot2Cmd1;
    G4UIcmdWithADoubleAndUnit  *halfxCmd1;
    G4UIcmdWithADoubleAndUnit  *halfyCmd1;
    G4UIcmdWithADoubleAndUnit  *halfzCmd1;
    G4UIcmdWithADoubleAndUnit  *radiusCmd1;
    G4UIcmdWithADoubleAndUnit  *radius0Cmd1;
    G4UIcmdWithADoubleAndUnit  *possigmarCmd1;
    G4UIcmdWithADoubleAndUnit  *possigmaxCmd1;
    G4UIcmdWithADoubleAndUnit  *possigmayCmd1;
    G4UIcmdWithADoubleAndUnit  *paralpCmd1;
    G4UIcmdWithADoubleAndUnit  *partheCmd1;
    G4UIcmdWithADoubleAndUnit  *parphiCmd1;  
    G4UIcmdWithAString         *confineCmd1;
         
//  //old ones, will be reomved soon
//  G4UIcmdWithAString         *typeCmd;
//  G4UIcmdWithAString         *shapeCmd;
//  G4UIcmdWith3VectorAndUnit  *centreCmd;
//  G4UIcmdWith3Vector         *posrot1Cmd;
//  G4UIcmdWith3Vector         *posrot2Cmd;
//  G4UIcmdWithADoubleAndUnit  *halfxCmd;
//  G4UIcmdWithADoubleAndUnit  *halfyCmd;
//  G4UIcmdWithADoubleAndUnit  *halfzCmd;
//  G4UIcmdWithADoubleAndUnit  *radiusCmd;
//  G4UIcmdWithADoubleAndUnit  *radius0Cmd;
//  G4UIcmdWithADoubleAndUnit  *possigmarCmd;
//  G4UIcmdWithADoubleAndUnit  *possigmaxCmd;
//  G4UIcmdWithADoubleAndUnit  *possigmayCmd;
//  G4UIcmdWithADoubleAndUnit  *paralpCmd;
//  G4UIcmdWithADoubleAndUnit  *partheCmd;
//  G4UIcmdWithADoubleAndUnit  *parphiCmd;  
//  G4UIcmdWithAString         *confineCmd; 
    
    // angular commands
    G4UIdirectory* angularDirectory;
    G4UIcmdWithAString         *angtypeCmd1;
    G4UIcmdWith3Vector         *angrot1Cmd1;
    G4UIcmdWith3Vector         *angrot2Cmd1;
    G4UIcmdWithADoubleAndUnit  *minthetaCmd1;
    G4UIcmdWithADoubleAndUnit  *maxthetaCmd1;
    G4UIcmdWithADoubleAndUnit  *minphiCmd1;
    G4UIcmdWithADoubleAndUnit  *maxphiCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmarCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmaxCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmayCmd1;
    G4UIcmdWith3VectorAndUnit  *angfocusCmd;
    G4UIcmdWithABool           *useuserangaxisCmd1;
    G4UIcmdWithABool           *surfnormCmd1;

  // old ones, will be removed soon
//  G4UIcmdWithAString         *angtypeCmd;
//  G4UIcmdWith3Vector         *angrot1Cmd;
//  G4UIcmdWith3Vector         *angrot2Cmd;
//  G4UIcmdWithADoubleAndUnit  *minthetaCmd;
//  G4UIcmdWithADoubleAndUnit  *maxthetaCmd;
//  G4UIcmdWithADoubleAndUnit  *minphiCmd;
//  G4UIcmdWithADoubleAndUnit  *maxphiCmd;
//  G4UIcmdWithADoubleAndUnit  *angsigmarCmd;
//  G4UIcmdWithADoubleAndUnit  *angsigmaxCmd;
//  G4UIcmdWithADoubleAndUnit  *angsigmayCmd;
//  G4UIcmdWithABool           *useuserangaxisCmd;
//  G4UIcmdWithABool           *surfnormCmd;

    // energy commands
    G4UIdirectory* energyDirectory;
    G4UIcmdWithAString         *energytypeCmd1;
    G4UIcmdWithADoubleAndUnit  *eminCmd1;
    G4UIcmdWithADoubleAndUnit  *emaxCmd1;
    G4UIcmdWithADoubleAndUnit  *monoenergyCmd1;
    G4UIcmdWithADoubleAndUnit  *engsigmaCmd1;
    G4UIcmdWithADouble         *alphaCmd1;
    G4UIcmdWithADouble         *tempCmd1;
    G4UIcmdWithADouble         *ezeroCmd1;
    G4UIcmdWithADouble         *gradientCmd1;
    G4UIcmdWithADouble         *interceptCmd1;
    G4UIcmdWithADouble         *arbeintCmd1;
    G4UIcmdWithoutParameter    *calculateCmd1;
    G4UIcmdWithABool           *energyspecCmd1;
    G4UIcmdWithABool           *diffspecCmd1;

    // old ones, will be removed soon
//    G4UIcmdWithAString         *energytypeCmd;
//    G4UIcmdWithADoubleAndUnit  *eminCmd;
//    G4UIcmdWithADoubleAndUnit  *emaxCmd;
//    G4UIcmdWithADoubleAndUnit  *monoenergyCmd;
//    G4UIcmdWithADoubleAndUnit  *engsigmaCmd;
//    G4UIcmdWithADouble         *alphaCmd;
//    G4UIcmdWithADouble         *tempCmd;
//    G4UIcmdWithADouble         *ezeroCmd;
//    G4UIcmdWithADouble         *gradientCmd;
//    G4UIcmdWithADouble         *interceptCmd;
//    G4UIcmdWithoutParameter    *calculateCmd;
//    G4UIcmdWithABool           *energyspecCmd;
//    G4UIcmdWithABool           *diffspecCmd;

    // histogram commands
    G4UIdirectory              *histDirectory;
    // old ones, will be removed soon
    //G4UIcmdWith3Vector         *histpointCmd;
    //G4UIcmdWithAString         *histnameCmd;
    //G4UIcmdWithAString         *arbintCmd;
    //G4UIcmdWithAString         *resethistCmd;
    G4UIcmdWith3Vector         *histpointCmd1;
    G4UIcmdWithAString         *histfileCmd1;
    G4UIcmdWithAString         *histnameCmd1;
    G4UIcmdWithAString         *arbintCmd1;
    G4UIcmdWithAString         *resethistCmd1;

    G4UIcmdWithAnInteger* verbosityCmd;

    // Commands from G4ParticleGun
    G4UIcommand* ionCmd;
    G4UIcommand* ionLvlCmd;
    G4UIcmdWithAString* particleCmd;
    G4UIcmdWithADoubleAndUnit* timeCmd;
    G4UIcmdWith3Vector* polCmd;
    G4UIcmdWithAnInteger* numberCmd;
    G4UIcmdWith3VectorAndUnit* positionCmd;
    G4UIcmdWith3Vector* directionCmd;
    G4UIcmdWithADoubleAndUnit* energyCmd;
    G4UIcmdWithoutParameter* listCmd;

  private: // for ion shooting
    G4bool   fShootIon; 

    G4int    fAtomicNumber;
    G4int    fAtomicMass;
    G4int    fIonCharge;
    G4double fIonExciteEnergy;

    G4int    fAtomicNumberL;
    G4int    fAtomicMassL;
    G4int    fIonChargeL;
    G4int    fIonEnergyLevel;
};

#endif

