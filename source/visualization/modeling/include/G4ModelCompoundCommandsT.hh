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
// $Id: G4ModelCompoundCommandsT.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay September 2006
//
// Compound model commands.
//
#ifndef G4MODELCOMPOUNDCOMMANDST_HH
#define G4MODELCOMPOUNDCOMMANDST_HH

#include "G4ModelApplyCommandsT.hh"
#include "G4ModelCommandUtils.hh"
#include "G4ModelCommandsT.hh"
#include "G4String.hh"

////////////////////////////////////////////////////////////////////////
// Set interval context
template <typename M>
class G4ModelCmdAddIntervalContext: public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdAddIntervalContext(M* model, const G4String& placement, 
			    const G4String& cmdName="addInterval")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Add interval.");
  }
  
  virtual ~G4ModelCmdAddIntervalContext() {
    std::vector<G4UImessenger*>::iterator iter = fMessengers.begin();
    
    while (iter != fMessengers.end()) {
      delete *iter;
      iter++;
    }
  }
  
protected:
  
  virtual void Apply(const G4String& param) {
    G4String myString(param);

    G4String name;
    std::istringstream is(param);
    
    is >> name;

    myString.erase(0, name.size());

    G4String dir = G4VModelCommand<M>::Placement()+"/"+G4VModelCommand<M>::Model()->Name();
    
    G4VisTrajContext* context = new G4VisTrajContext(name);
    
    G4ModelCommandUtils::AddContextMsgrs(context, fMessengers, dir);
    G4VModelCommand<M>::Model()->AddIntervalContext(myString, context);
  }

private:
  
  std::vector<G4UImessenger*> fMessengers;

};

////////////////////////////////////////////////////////////////////////
// Set value context
template <typename M>
class G4ModelCmdAddValueContext: public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdAddValueContext(M* model, const G4String& placement, 
			    const G4String& cmdName="addValue")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Add value.");
  }
  
  virtual ~G4ModelCmdAddValueContext() {
    std::vector<G4UImessenger*>::iterator iter = fMessengers.begin();
    
    while (iter != fMessengers.end()) {
      delete *iter;
      iter++;
    }
  }

protected:

  virtual void Apply(const G4String& param) {
    G4String myString(param);

    G4String name;
    std::istringstream is(param);
    
    is >> name;

    myString.erase(0, name.size());

    G4String dir = G4VModelCommand<M>::Placement()+"/"+G4VModelCommand<M>::Model()->Name();
    
    G4VisTrajContext* context = new G4VisTrajContext(name);
    
    G4ModelCommandUtils::AddContextMsgrs(context, fMessengers, dir);
    G4VModelCommand<M>::Model()->AddValueContext(myString, context);
  }

private:
  
  std::vector<G4UImessenger*> fMessengers;

};

#endif
