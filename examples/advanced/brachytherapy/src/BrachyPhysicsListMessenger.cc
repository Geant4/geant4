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

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/qgsm/include/G4ASCCrossSection.hh
class G4ASCCrossSection : public G4VAnnihilationCrossSection
{
public:
	G4ASCCrossSection(G4int, G4int, G4double,  G4double, G4double, G4double);
	G4bool InCharge(G4int aCode, G4int bCode);
	G4double GetXsec(G4double S);
private:

	G4int theCode1;
	G4int theCode2;
	G4double theX;
	G4double theY;
	G4double theEta;
	G4double theEps;
};


inline G4bool G4ASCCrossSection::
InCharge(G4int aCode, G4int bCode)
{
	G4bool result;
	result = (aCode==theCode1&&bCode==theCode2)||(aCode==theCode2&&bCode==theCode1);
	return result;
}

inline G4double G4ASCCrossSection::
GetXsec(G4double S)
{
	G4double result = theX*G4Pow::GetInstance()->powA(S, theEps) + 
                          theY*G4Pow::GetInstance()->powA(S, -theEta);
	return result;
}

#endif
=======
#include "BrachyPhysicsListMessenger.hh"
#include "BrachyPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

BrachyPhysicsListMessenger::BrachyPhysicsListMessenger(BrachyPhysicsList* pPhys)
:G4UImessenger(),
 fPhysicsList(pPhys),fPhysDir(0),fListCmd(0)
{ 
  fPhysDir = new G4UIdirectory("/testem/phys/");
  fPhysDir->SetGuidance("physics list commands");


  fListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);  
}

BrachyPhysicsListMessenger::~BrachyPhysicsListMessenger()
{
  delete fListCmd;
  delete fPhysDir;
}

void BrachyPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == fListCmd )
   { fPhysicsList->AddPhysicsList(newValue);}
}
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/advanced/brachytherapy/src/BrachyPhysicsListMessenger.cc
