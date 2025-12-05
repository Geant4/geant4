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
// G4QSSMessenger implementation
//
// Author: Leandro Gomez Vidal (Univ. Buenos Aires), October 2021
// --------------------------------------------------------------------

#include "G4QSSMessenger.hh"
#include "G4QSSParameters.hh"

G4QSSMessenger::G4QSSMessenger()
{
  commandsShouldBeInMaster = true;

  qssCmdDir = new G4UIdirectory("/QSS/", false);
  qssCmdDir->SetGuidance("G4QSStepper configuration.");

  dQMinCmd = new G4UIcmdWithADoubleAndUnit("/QSS/dQMin",this);
  dQMinCmd->SetDefaultUnit("mm");
  dQMinCmd->SetParameterName("dQMinCmd",false);
  dQMinCmd->SetUnitCategory("Length");

  dQRelCmd = new G4UIcmdWithADouble("/QSS/dQRel",this);
  dQRelCmd->SetGuidance("Default is 1e-5");
  dQRelCmd->SetParameterName("dQRelCmd",false);

  stepperSelectorCmd = new G4UIcmdWithAString("/QSS/selectStepper",this);
  stepperSelectorCmd->SetGuidance("Select stepper.");
  stepperSelectorCmd->SetParameterName("choice", false);
  stepperSelectorCmd->SetCandidates("TemplatedDoPri OldRK45 G4QSS2");

  maxSubstepsCmd = new G4UIcmdWithAnInteger("/QSS/maxSubsteps",this);
  maxSubstepsCmd->SetGuidance("Default is 5000");
  maxSubstepsCmd->SetDefaultValue(5000);
  maxSubstepsCmd->SetParameterName("maxSubstepsCmd", false);

}

G4QSSMessenger::~G4QSSMessenger()
{
  delete qssCmdDir;
  delete dQMinCmd;
  delete dQRelCmd;
  delete stepperSelectorCmd;
  delete maxSubstepsCmd;
}

G4QSSMessenger* G4QSSMessenger::instance()
{
  static G4QSSMessenger theSingleMessengerInstance;
  return &theSingleMessengerInstance;
}

G4bool  G4QSSMessenger::Set_dQMin( G4double dvalue )
{ 
  G4bool good_value= true;
  good_value= G4QSSParameters::Instance()->Set_dQMin( dvalue );
  return good_value;
}

G4bool G4QSSMessenger::Set_dQRel( G4double dvalue )
{
  G4bool good_value= true;
  good_value= G4QSSParameters::Instance()->Set_dQRel( dvalue );
  return good_value;
}

G4bool G4QSSMessenger::SetMaxSubsteps( G4int value )
{
  G4bool good_value= true;
  
  good_value = G4QSSParameters::Instance()->SetMaxSubsteps( value );  

  return good_value;  
}

void G4QSSMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  G4bool good_value= true;

  if ( command == dQMinCmd )
  {
    good_value= Set_dQMin( dQMinCmd->GetNewDoubleValue(newValue) );
  }

  if ( command == dQRelCmd )
  {
    good_value= Set_dQRel( dQRelCmd->GetNewDoubleValue(newValue) );
  }

  if (command == maxSubstepsCmd)
  {
    G4int oldMaxSubstp = GetMaxSubsteps();
    good_value= SetMaxSubsteps ( maxSubstepsCmd->GetNewIntValue(newValue) );
    if( good_value )
    {
      G4cout << "G4QSSMessenger: changed maxSubsteps = " << GetMaxSubsteps()
             << "  ( old value = " << oldMaxSubstp << " )." << G4endl;
    }
  }

  if( !good_value ) 
  {
    G4cerr << " G4QSSMessenger: Command FAILED - parameter was not accepted." << G4endl;
  }

  if(command == stepperSelectorCmd)
  {
    this->selectStepper(newValue);
  }
}

void G4QSSMessenger::selectStepper(const std::string &newValue)
{
  const std::map<std::string, StepperSelection> stepperMapping =
    {
      {"G4QSS2", G4QSS2}, {"G4QSS3", G4QSS3}
    };
  _selectedStepper = stepperMapping.at(newValue);

  if( (StepperSelection::None <= _selectedStepper ) 
      && (_selectedStepper < StepperSelection::NumMethods ) )
  {
    G4cout << "G4QSSMessenger: Selecting stepper " << newValue 
        <<   " ( which is method # " << _selectedStepper << " ) " << G4endl;

    G4int new_order= 0;
    if( _selectedStepper == G4QSS2 )
    {
      new_order= 2;
    }
    else
    {
      if( _selectedStepper == G4QSS3 )
      {
        new_order= 3;
      }
    }
    SetQssOrder( new_order ); // Will check requested value
  }
  else
  {
    G4cerr << "G4QSSMessenger::selectStepper : Invalid type of stepper requested" 
           << " Cannot find the requested stepper type '" << newValue << G4endl;
  }
}

G4bool G4QSSMessenger::SetQssOrder(G4int value)
{
  G4bool good_value= false;
  
  if( 2 <= value && value <= 3 )
  {
    if( value == 2)
    { 
       _selectedStepper = G4QSS2;
    }
    else
    {
       _selectedStepper = G4QSS3;
    }
    good_value= G4QSSParameters::Instance()->SetQssOrder( value );
  }
  else
  {
    G4cerr << "G4QSSMessenger::SetQssOrder : requested invalid *order* of QSS-stepper"
           << " Asked for " << value << " . Value must be either 2 or 3." << G4endl;
  }
  return good_value;  
}

G4QSSMessenger::StepperSelection G4QSSMessenger::selectedStepper()
{
  return _selectedStepper;
}
