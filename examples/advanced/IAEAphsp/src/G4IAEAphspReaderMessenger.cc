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
// Author: M.A. Cortes-Giraldo, Universidad de Sevilla
//
// History changelog prior creation of this example:
// - 17/10/2009: version 1.0
// - 20/11/2009: version 1.1 before publishing:
//   - Changed some names by more suitable ones
// - 02/08/2010: version 1.2-dev:
//   - Added possibility of applying axial symmetries
// - 14/09/2023: version 2.0
//   - Following Geant4 coding guidelines
// - 18/10/2025: version 3.0
//   - Creation of IAEASourceIdRegistry for thread-safe source_id assignation
//


#include "G4IAEAphspReaderMessenger.hh"
#include "G4IAEAphspReader.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"


G4IAEAphspReaderMessenger::G4IAEAphspReaderMessenger(G4IAEAphspReader* reader)
  :fIAEAphspReader(reader)
{
  fPhaseSpaceDir = new G4UIdirectory("/IAEAphspReader/");
  fPhaseSpaceDir
    ->SetGuidance("Commands for the IAEA phase-space file management.");


  fVerboseCmd = new G4UIcmdWithAnInteger("/IAEAphspReader/verbose", this);
  fVerboseCmd->SetGuidance("Set verbose level of G4IAEAphspReader class");
  fVerboseCmd->SetParameterName("value", false);
  fVerboseCmd->SetRange("value >= 0");
  fVerboseCmd->AvailableForStates(G4State_Idle);

  fNofParallelRunsCmd =
    new G4UIcmdWithAnInteger("/IAEAphspReader/numberOfParallelRuns", this);
  fNofParallelRunsCmd
    ->SetGuidance("Select the number of fragments N in which the phase-space");
  fNofParallelRunsCmd
    ->SetGuidance(" file is divided into.");
  fNofParallelRunsCmd->SetParameterName("N", false);
  fNofParallelRunsCmd->SetRange("N > 0");
  fNofParallelRunsCmd->AvailableForStates(G4State_Idle);

  fParallelRunCmd =
    new G4UIcmdWithAnInteger("/IAEAphspReader/parallelRun", this);
  fParallelRunCmd->
    SetGuidance("Use the fragment F (of a total of N) from which particles");
  fParallelRunCmd->SetGuidance(" are extracted. (1 <= F <= N).");
  fParallelRunCmd->SetParameterName("frag", false);
  fParallelRunCmd->SetRange("frag > 0");
  fParallelRunCmd->AvailableForStates(G4State_Idle);

  fTimesRecycledCmd =
    new G4UIcmdWithAnInteger("/IAEAphspReader/recycling", this);
  fTimesRecycledCmd
    ->SetGuidance("Select the number of times that each particle is reused.");
  fTimesRecycledCmd
    ->SetGuidance("(Caution: 1 means that each particle is used twice.)");
  fTimesRecycledCmd->SetParameterName("choice", false);
  fTimesRecycledCmd->SetRange("choice >= 0");
  fTimesRecycledCmd->AvailableForStates(G4State_Idle);

  fPhspGlobalTranslationCmd = 
    new G4UIcmdWith3VectorAndUnit("/IAEAphspReader/translate", this);
  fPhspGlobalTranslationCmd->SetGuidance("Set the translation components.");
  fPhspGlobalTranslationCmd->SetParameterName("x0", "y0", "z0", false);
  fPhspGlobalTranslationCmd->SetDefaultUnit("cm");
  fPhspGlobalTranslationCmd->SetUnitCandidates("mm cm m");
  fPhspGlobalTranslationCmd->AvailableForStates(G4State_Idle);

  fPhspRotationOrderCmd = 
    new G4UIcmdWithAnInteger("/IAEAphspReader/rotationOrder", this);
  fPhspRotationOrderCmd
    ->SetGuidance("Select the order in which the rotations are performed.");
  fPhspRotationOrderCmd
    ->SetGuidance("1 means X axis, 2 means Y axis and 3 means Z axis.");
  fPhspRotationOrderCmd
    ->SetGuidance("The argument must be a 3-digit integer without repetition.");
  fPhspRotationOrderCmd->SetParameterName("choice", false);
  fPhspRotationOrderCmd->AvailableForStates(G4State_Idle);

  fRotXCmd = new G4UIcmdWithADoubleAndUnit("/IAEAphspReader/rotateX", this);
  fRotXCmd->SetGuidance("Set the rotation angle around the global X axis.");
  fRotXCmd->SetParameterName("angle", false);
  fRotXCmd->SetDefaultUnit("deg");
  fRotXCmd->SetUnitCandidates("deg rad");
  fRotXCmd->AvailableForStates(G4State_Idle);

  fRotYCmd = new G4UIcmdWithADoubleAndUnit("/IAEAphspReader/rotateY", this);
  fRotYCmd->SetGuidance("Set the rotation angle around the global Y axis.");
  fRotYCmd->SetParameterName("angle", false);
  fRotYCmd->SetDefaultUnit("deg");
  fRotYCmd->SetUnitCandidates("deg rad");
  fRotYCmd->AvailableForStates(G4State_Idle);

  fRotZCmd = new G4UIcmdWithADoubleAndUnit("/IAEAphspReader/rotateZ", this);
  fRotZCmd->SetGuidance("Set the rotation angle around the global Z axis.");
  fRotZCmd->SetParameterName("angle", false);
  fRotZCmd->SetDefaultUnit("deg");
  fRotZCmd->SetUnitCandidates("deg rad");
  fRotZCmd->AvailableForStates(G4State_Idle);

  fIsocenterPosCmd = 
    new G4UIcmdWith3VectorAndUnit("/IAEAphspReader/isocenterPosition", this);
  fIsocenterPosCmd->SetGuidance("Set the isocenter position.");
  fIsocenterPosCmd->SetParameterName("Xic", "Yic", "Zic", false);
  fIsocenterPosCmd->SetDefaultUnit("cm");
  fIsocenterPosCmd->SetUnitCandidates("mm cm m");
  fIsocenterPosCmd->AvailableForStates(G4State_Idle);

  fCollimatorRotAxisCmd = 
    new G4UIcmdWith3Vector("/IAEAphspReader/collimatorRotationAxis", this);
  fCollimatorRotAxisCmd
    ->SetGuidance("Set the rotation axis of the collimator.");
  fCollimatorRotAxisCmd->SetGuidance("It has to be a unit vector.");
  fCollimatorRotAxisCmd->SetParameterName("Ucol", "Vcol", "Wcol", false);
  fCollimatorRotAxisCmd->SetRange("Ucol != 0 || Vcol != 0 || Wcol != 0");
  fCollimatorRotAxisCmd->AvailableForStates(G4State_Idle);

  fCollimatorAngleCmd = 
    new G4UIcmdWithADoubleAndUnit("/IAEAphspReader/collimatorAngle", this);
  fCollimatorAngleCmd
    ->SetGuidance("Set the rotation angle of the phase space plane around ");
  fCollimatorAngleCmd
    ->SetGuidance("the rotation axis of the collimator.");
  fCollimatorAngleCmd->SetParameterName("angle", false);
  fCollimatorAngleCmd->SetDefaultUnit("deg");
  fCollimatorAngleCmd->SetUnitCandidates("deg rad");
  fCollimatorAngleCmd->AvailableForStates(G4State_Idle);

  fGantryRotAxisCmd = 
    new G4UIcmdWith3Vector("/IAEAphspReader/gantryRotationAxis", this);
  fGantryRotAxisCmd->SetGuidance("Set the rotation axis of the gantry.");
  fGantryRotAxisCmd->SetGuidance("It has to be a unit vector.");
  fGantryRotAxisCmd->SetParameterName("Ugan", "Vgan", "Wgan", false);
  fGantryRotAxisCmd->SetRange("Ugan != 0 || Vgan != 0 || Wgan != 0");
  fGantryRotAxisCmd->AvailableForStates(G4State_Idle);

  fGantryAngleCmd = 
    new G4UIcmdWithADoubleAndUnit("/IAEAphspReader/gantryAngle", this);
  fGantryAngleCmd
    ->SetGuidance("Set the rotation angle of the phase space plane around ");
  fGantryAngleCmd
    ->SetGuidance("the gantry rotation axis.");
  fGantryAngleCmd->SetParameterName("angle", false);
  fGantryAngleCmd->SetDefaultUnit("deg");
  fGantryAngleCmd->SetUnitCandidates("deg rad");
  fGantryAngleCmd->AvailableForStates(G4State_Idle);

  fAxialSymmetryXCmd =
    new G4UIcmdWithABool("/IAEAphspReader/axialSymmetryX", this);
  fAxialSymmetryXCmd
    ->SetGuidance("Command to take into account rotational symmetry around X");
  fAxialSymmetryXCmd->SetParameterName("choice", true);
  fAxialSymmetryXCmd->SetDefaultValue(true);
  fAxialSymmetryXCmd->AvailableForStates(G4State_Idle);

  fAxialSymmetryYCmd =
    new G4UIcmdWithABool("/IAEAphspReader/axialSymmetryY", this);
  fAxialSymmetryYCmd
    ->SetGuidance("Command to take into account rotational symmetry around Y");
  fAxialSymmetryYCmd->SetParameterName("choice", true);
  fAxialSymmetryYCmd->SetDefaultValue(true);
  fAxialSymmetryYCmd->AvailableForStates(G4State_Idle);

  fAxialSymmetryZCmd =
    new G4UIcmdWithABool("/IAEAphspReader/axialSymmetryZ", this);
  fAxialSymmetryZCmd
    ->SetGuidance("Command to take into account rotational symmetry around Z");
  fAxialSymmetryZCmd->SetParameterName("choice", true);
  fAxialSymmetryZCmd->SetDefaultValue(true);
  fAxialSymmetryZCmd->AvailableForStates(G4State_Idle);
}


G4IAEAphspReaderMessenger::~G4IAEAphspReaderMessenger()
{ 
  delete fPhaseSpaceDir;
  delete fVerboseCmd;
  delete fNofParallelRunsCmd;
  delete fParallelRunCmd;
  delete fTimesRecycledCmd;
  delete fPhspGlobalTranslationCmd;
  delete fPhspRotationOrderCmd;
  delete fRotXCmd;
  delete fRotYCmd;
  delete fRotZCmd;
  delete fIsocenterPosCmd;
  delete fCollimatorRotAxisCmd;
  delete fCollimatorAngleCmd;
  delete fGantryRotAxisCmd;
  delete fGantryAngleCmd;
  delete fAxialSymmetryXCmd;
  delete fAxialSymmetryYCmd;
  delete fAxialSymmetryZCmd;
}


void G4IAEAphspReaderMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if( command == fVerboseCmd )
    fIAEAphspReader->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));

  else if( command == fNofParallelRunsCmd )
    fIAEAphspReader
      ->SetTotalParallelRuns(fNofParallelRunsCmd->GetNewIntValue(newValue));

  else if( command == fParallelRunCmd )
    fIAEAphspReader->SetParallelRun(fParallelRunCmd->GetNewIntValue(newValue));

  else if( command == fTimesRecycledCmd )
    fIAEAphspReader
      ->SetTimesRecycled(fTimesRecycledCmd->GetNewIntValue(newValue) );

  else if( command == fPhspGlobalTranslationCmd )
    fIAEAphspReader
      ->SetGlobalPhspTranslation(fPhspGlobalTranslationCmd
                                 ->GetNew3VectorValue(newValue) );

  else if( command == fPhspRotationOrderCmd )
    fIAEAphspReader
      ->SetTimesRecycled(fPhspRotationOrderCmd->GetNewIntValue(newValue) );

  else if( command == fRotXCmd )
    fIAEAphspReader->SetRotationX( fRotXCmd->GetNewDoubleValue(newValue) );

  else if( command == fRotYCmd )
    fIAEAphspReader->SetRotationY( fRotYCmd->GetNewDoubleValue(newValue) );

  else if( command == fRotZCmd )
    fIAEAphspReader->SetRotationZ( fRotZCmd->GetNewDoubleValue(newValue) );

  else if( command == fIsocenterPosCmd )
    fIAEAphspReader
      ->SetIsocenterPosition(fIsocenterPosCmd->GetNew3VectorValue(newValue));

  else if( command == fCollimatorRotAxisCmd )
    fIAEAphspReader
      ->SetCollimatorRotationAxis(fCollimatorRotAxisCmd
                                  ->GetNew3VectorValue(newValue) );

  else if( command == fCollimatorAngleCmd )
    fIAEAphspReader
      ->SetCollimatorAngle(fCollimatorAngleCmd->GetNewDoubleValue(newValue));

  else if( command == fGantryRotAxisCmd )
    fIAEAphspReader
      ->SetGantryRotationAxis(fGantryRotAxisCmd->GetNew3VectorValue(newValue));

  else if( command == fGantryAngleCmd )
    fIAEAphspReader
      ->SetGantryAngle( fGantryAngleCmd->GetNewDoubleValue(newValue) );

  else if( command == fAxialSymmetryXCmd )
    fIAEAphspReader
      ->SetAxialSymmetryX( fAxialSymmetryXCmd->GetNewBoolValue(newValue));

  else if( command == fAxialSymmetryYCmd )
    fIAEAphspReader
      ->SetAxialSymmetryY( fAxialSymmetryYCmd->GetNewBoolValue(newValue) );

  else if( command == fAxialSymmetryZCmd )
    fIAEAphspReader
      ->SetAxialSymmetryZ( fAxialSymmetryZCmd->GetNewBoolValue(newValue) );

}
