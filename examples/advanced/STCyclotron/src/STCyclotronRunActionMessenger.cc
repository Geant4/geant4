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
//  Author: F. Poignant, floriane.poignant@gmail.com
//


#include "STCyclotronRunActionMessenger.hh"
#include "STCyclotronRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4SystemOfUnits.hh"

/////////////////////////////////////////////////////////////////////////////
STCyclotronRunActionMessenger::STCyclotronRunActionMessenger(STCyclotronRunAction* run)
  :fG4Run(run)
{
    // change the time of irradiation
    fIrradiationTime = new G4UIdirectory("/setTimeOfIrradiation/");
    fIrradiationTime -> SetGuidance("Change time of irradiation, in hour(s).");

    fChangeIrradiationTimeCmd = new G4UIcmdWithADouble("/setTimeOfIrradiation/time", this);
    fChangeIrradiationTimeCmd -> SetGuidance("Change the value of the time of irradiation (in hours)"
	                                   "\nDefault value is 6 hours");
    fChangeIrradiationTimeCmd -> SetParameterName("TimeOfIrradiation", true);
    fChangeIrradiationTimeCmd -> SetRange("TimeOfIrradiation > 0.");
    fChangeIrradiationTimeCmd -> SetDefaultValue(6.);
    
   }

/////////////////////////////////////////////////////////////////////////////
STCyclotronRunActionMessenger::~STCyclotronRunActionMessenger()
{
    delete fIrradiationTime; 
    delete fChangeIrradiationTimeCmd; 
    
}

/////////////////////////////////////////////////////////////////////////////
void STCyclotronRunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  
  if( command == fChangeIrradiationTimeCmd)
    {
      G4double updatedValue  = fChangeIrradiationTimeCmd -> GetNewDoubleValue(newValue);
      fG4Run -> SetIrradiationTime(updatedValue);
    }
  
}
