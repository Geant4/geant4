#ifndef G4GeneralParticleSourceMessenger_h
#define G4GeneralParticleSourceMessenger_h 1
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSourceMessenger.hh
//
// Version:      1.1
// Date:         19/10/00
// Author:       C Ferguson, F Lei and P Truscott
// Organisation: University of Southampton / DERA
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// DESCRIPTION
// -----------
//
// The function of the G4GeneralParticleSourceMessenger is to allow the user to
// enter commands either in interactive command line mode or through macros to
// control the G4GeneralParticleSource. The G4GeneralParticleSourceMessenger
// class is based on G4ParticleGunMessenger.
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 28 February 2000, C Ferguson, Created.
//
// Version 1.1, 19 October 2000, Modified to inherit from G4VPrimaryGenerator.
// New name at the request of M. Asai.
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
// void SetNewValue(G4UIcommand *command, G4String newValues)
//     Uses the appropriate methods in the G4GeneralParticleSource to carry out
//     the user commands.
//
// G4String GetCurrentValue(G4UIcommand *command)
//     Allows the user to retrieve the current values of parameters.
//     Not implemented yet.
//
///////////////////////////////////////////////////////////////////////////////
#include "G4UImessenger.hh"
#include "globals.hh"
//#include "UIcmdWithNucleusAndUnit.hh"

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

class G4GeneralParticleSource;

class G4GeneralParticleSourceMessenger: public G4UImessenger
{
public:
  G4GeneralParticleSourceMessenger(G4GeneralParticleSource *fPtclGun);
  ~G4GeneralParticleSourceMessenger();
    
  void SetNewValue(G4UIcommand *command, G4String newValues);
  //    Identifies the command which has been invoked by the user, extracts the
  //    parameters associated with that command (held in newValues), and uses
  //    these values with the appropriate member function of G4GeneralParticleSource.
  G4String GetCurrentValue(G4UIcommand *command);

private:
  G4GeneralParticleSource *fParticleGun;
  G4ParticleTable *particleTable;
  G4String histtype;
    
private: //commands
  G4UIdirectory              *gpsDirectory;

  G4UIcmdWithAString         *typeCmd;
  G4UIcmdWithAString         *shapeCmd;
  G4UIcmdWith3VectorAndUnit  *centreCmd;
  G4UIcmdWith3Vector         *posrot1Cmd;
  G4UIcmdWith3Vector         *posrot2Cmd;
  G4UIcmdWithADoubleAndUnit  *halfxCmd;
  G4UIcmdWithADoubleAndUnit  *halfyCmd;
  G4UIcmdWithADoubleAndUnit  *halfzCmd;
  G4UIcmdWithADoubleAndUnit  *radiusCmd;
  G4UIcmdWithADoubleAndUnit  *radius0Cmd;
  G4UIcmdWithADoubleAndUnit  *paralpCmd;
  G4UIcmdWithADoubleAndUnit  *partheCmd;
  G4UIcmdWithADoubleAndUnit  *parphiCmd;  
  G4UIcmdWithAString         *confineCmd;         

  G4UIcmdWithAString         *angtypeCmd;
  G4UIcmdWith3Vector         *angrot1Cmd;
  G4UIcmdWith3Vector         *angrot2Cmd;
  G4UIcmdWithADoubleAndUnit  *minthetaCmd;
  G4UIcmdWithADoubleAndUnit  *maxthetaCmd;
  G4UIcmdWithADoubleAndUnit  *minphiCmd;
  G4UIcmdWithADoubleAndUnit  *maxphiCmd;
  G4UIcmdWithABool           *surfnormCmd;

  G4UIcmdWithAString         *energytypeCmd;
  G4UIcmdWithADoubleAndUnit  *eminCmd;
  G4UIcmdWithADoubleAndUnit  *emaxCmd;
  G4UIcmdWithADoubleAndUnit  *monoenergyCmd;
  G4UIcmdWithADouble         *alphaCmd;
  G4UIcmdWithADouble         *tempCmd;
  G4UIcmdWithADouble         *ezeroCmd;
  G4UIcmdWithADouble         *gradientCmd;
  G4UIcmdWithADouble         *interceptCmd;
  G4UIcmdWithoutParameter    *calculateCmd;
  G4UIcmdWithABool           *energyspecCmd;
  G4UIcmdWithABool           *diffspecCmd;

  G4UIcmdWith3Vector         *histpointCmd;
  G4UIcmdWithAString         *histnameCmd;
  G4UIcmdWithAString         *arbintCmd;

  G4UIcmdWithAnInteger       *verbosityCmd;

  G4UIcommand                *ionCmd;

  G4UIcmdWithAString         *particleCmd;
  G4UIcmdWithADoubleAndUnit  *timeCmd;
  G4UIcmdWith3Vector         *polCmd;
  G4UIcmdWithAnInteger       *numberCmd;

  private: // for ion shooting
    G4bool   fShootIon; 
    G4int    fAtomicNumber;
    G4int    fAtomicMass;
    G4int    fIonCharge;
    G4double fIonExciteEnergy;

};

#endif
