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
// $Id: G4ModelCommandsT.hh,v 1.7 2006/06/29 21:30:24 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Generic model messenges. 
//
// Jane Tinslay March 2006
//
#ifndef G4MODELCOMMANDST_HH
#define G4MODELCOMMANDST_HH

#include "G4ModelApplyCommandsT.hh"
#include "G4Polymarker.hh"
#include "G4UIdirectory.hh"
#include <sstream>

////////////////////////////////////////////////////////////////////////
// Set parameter colour
template <typename M>
class G4ModelCmdSetStringColour : public G4ModelCmdApplyStringColour<M> {

public: // With description

  G4ModelCmdSetStringColour(M* model, const G4String& placement)
    :G4ModelCmdApplyStringColour<M>(model, placement, "set") {}
  
  virtual ~G4ModelCmdSetStringColour() {}

protected:

  virtual void Apply(const G4String& param, const G4Colour& colour) {
    G4VModelCommand<M>::Model()->Set(param, colour);
  }

};

////////////////////////////////////////////////////////////////////////
// Set default colour
template <typename M>
class G4ModelCmdSetDefaultColour : public G4ModelCmdApplyColour<M> {

public: // With description

  G4ModelCmdSetDefaultColour(M* model, const G4String& placement)
    :G4ModelCmdApplyColour<M>(model, placement, "setDefault") {}
  
  virtual ~G4ModelCmdSetDefaultColour() {}

protected:

  virtual void Apply(const G4Colour& colour) {
    G4VModelCommand<M>::Model()->SetDefault(colour);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Add string command
template <typename M>
class G4ModelCmdAddString : public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdAddString(M* model, const G4String& placement)
    :G4ModelCmdApplyString<M>(model, placement, "add") 
  {
    G4ModelCmdApplyString<M>::Command()->SetGuidance("Add command");
  }

  virtual ~G4ModelCmdAddString() {}

protected:
  
  virtual void Apply(const G4String& newValue) {
    G4VModelCommand<M>::Model()->Add(newValue);
  }
  
};

////////////////////////////////////////////////////////////////////////
//Add integer command
template <typename M>
class G4ModelCmdAddInt : public G4ModelCmdApplyInteger<M> {

public: // With description

  G4ModelCmdAddInt(M* model, const G4String& placement)
  :G4ModelCmdApplyInteger<M>(model, placement, "add") 
  {
    G4ModelCmdApplyInteger<M>::Command()->SetGuidance("Add command");    
  }
  
  virtual ~G4ModelCmdAddInt() {}

protected:

  virtual void Apply(const G4int& newValue) {
    G4VModelCommand<M>::Model()->Add(newValue);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Invert command
template <typename M>
class G4ModelCmdInvert : public G4ModelCmdApplyBool<M> {
  
public: // With description
  
  G4ModelCmdInvert(M* model, const G4String& placement) 
    :G4ModelCmdApplyBool<M>(model, placement, "invert") 
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Invert command");
  }

  virtual ~G4ModelCmdInvert() {}

protected:

  virtual void Apply(const G4bool& newValue) {
    G4VModelCommand<M>::Model()->SetInvert(newValue);
  }

};

////////////////////////////////////////////////////////////////////////
// Active command
template <typename M>
class G4ModelCmdActive : public G4ModelCmdApplyBool<M> {

public: // With description

  G4ModelCmdActive(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "active") 
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Active command");
  }
  
  virtual ~G4ModelCmdActive() {}

protected:

  virtual void Apply(const G4bool& newValue) {
    G4VModelCommand<M>::Model()->SetActive(newValue);
  }

};

////////////////////////////////////////////////////////////////////////
// Verbose command
template <typename M>
class G4ModelCmdVerbose : public G4ModelCmdApplyBool<M> {

public: // With description

  G4ModelCmdVerbose(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "verbose")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Verbose command");
  }

  virtual ~G4ModelCmdVerbose() {}

protected:

  virtual void Apply(const G4bool& newValue) {
    G4VModelCommand<M>::Model()->SetVerbose(newValue);
  }

};

////////////////////////////////////////////////////////////////////////
// Reset command
template <typename M>
class G4ModelCmdReset : public G4ModelCmdApplyNull<M> {

public: // With description

  G4ModelCmdReset(M* model, const G4String& placement)
    :G4ModelCmdApplyNull<M>(model, placement, "reset") 
  {
    G4ModelCmdApplyNull<M>::Command()->SetGuidance("Reset command");    
  }
  
  virtual ~G4ModelCmdReset() {}
  
protected:

  virtual void Apply() {
    G4VModelCommand<M>::Model()->Reset();
  }
 
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points colour command
template <typename M>
class G4ModelCmdSetAuxPtsColour : public G4ModelCmdApplyColour<M> {
  
public:

  G4ModelCmdSetAuxPtsColour(M* model, const G4String& placement)
    :G4ModelCmdApplyColour<M>(model, placement, "setAuxPtsColour") {}
  
protected:

  void Apply(const G4Colour& colour) {
    G4VModelCommand<M>::Model()->SetAuxPtsColour(colour);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set set points colour command
template <typename M>
class G4ModelCmdSetStepPtsColour : public G4ModelCmdApplyColour<M> {
  
public:

  G4ModelCmdSetStepPtsColour(M* model, const G4String& placement)
    :G4ModelCmdApplyColour<M>(model, placement, "setStepPtsColour") {}
  
protected:
  
  void Apply(const G4Colour& colour) {
    G4VModelCommand<M>::Model()->SetStepPtsColour(colour);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set draw line command
template <typename M>
class G4ModelCmdSetDrawLine : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetDrawLine(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setDrawLine")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set draw line command");
  }
  
protected:

  void Apply(const G4bool& myBool) {
    G4VModelCommand<M>::Model()->SetDrawLine(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set line visibility command
template <typename M>
class G4ModelCmdSetLineVisible : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetLineVisible(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setLineVisible") 
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set line visibility command");
  }
  
protected:

   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetLineVisible(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set draw auxiliary points command
template <typename M>
class G4ModelCmdSetDrawAuxPts : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetDrawAuxPts(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setDrawAuxPts")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set draw auxiliary points command");
  }
  
protected:

   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetDrawAuxPts(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points visibility
template <typename M>
class G4ModelCmdSetAuxPtsVisible : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetAuxPtsVisible(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setAuxPtsVisible")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set auxiliary points visibility command");
  }
  
protected:

   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetAuxPtsVisible(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set draw step points command
template <typename M>
class G4ModelCmdSetDrawStepPts : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetDrawStepPts(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setDrawStepPts")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set draw step points command");
  }
  
protected:
  
   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetDrawStepPts(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points visible command
template <typename M>
class G4ModelCmdSetStepPtsVisible : public G4ModelCmdApplyBool<M> {
  
public:

  G4ModelCmdSetStepPtsVisible(M* model, const G4String& placement)
    :G4ModelCmdApplyBool<M>(model, placement, "setStepPtsVisible")
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set step points colour command");
  }
  
protected:

   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetStepPtsVisible(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points size command
template <typename M>
class G4ModelCmdSetAuxPtsSize : public G4ModelCmdApplyDouble<M> {
  
public:

  G4ModelCmdSetAuxPtsSize(M* model, const G4String& placement)
    :G4ModelCmdApplyDouble<M>(model, placement, "setAuxPtsSize")
  {
    G4ModelCmdApplyDouble<M>::Command()->SetGuidance("Set auxiliary points size command");
  }
  
protected:

   void Apply(const G4double& myDouble) {
     G4VModelCommand<M>::Model()->SetAuxPtsSize(myDouble);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points size command
template <typename M>
class G4ModelCmdSetStepPtsSize : public G4ModelCmdApplyDouble<M> {
  
public:

  G4ModelCmdSetStepPtsSize(M* model, const G4String& placement)
    :G4ModelCmdApplyDouble<M>(model, placement, "setStepPtsSize")
  {
    G4ModelCmdApplyDouble<M>::Command()->SetGuidance("Set step points colour command");
  }
  
protected:

   void Apply(const G4double& myDouble) {
     G4VModelCommand<M>::Model()->SetStepPtsSize(myDouble);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points type command
template <typename M>
class G4ModelCmdSetStepPtsType : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsType(M* model, const G4String& placement)
    :G4ModelCmdApplyString<M>(model, placement, "setStepPtsType")
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set step points type.");
    cmd->SetCandidates("dots circles squares");
  }
  
protected:

  void Apply(const G4String& type) {
    G4Polymarker::MarkerType myType;
    
    if (type == "dots") {myType = G4Polymarker::dots;}
    else if (type == "circles") {myType = G4Polymarker::circles;}
    else if (type == "squares") {myType = G4Polymarker::squares;}
    else {
      std::ostringstream o;
      o << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetStepPtsType::Apply",
	 "InvalidArgument", JustWarning, o.str().c_str());
      return;
    }
    G4VModelCommand<M>::Model()->SetStepPtsType(myType);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points type command
template <typename M>
class G4ModelCmdSetAuxPtsType : public G4ModelCmdApplyString<M> {
  
public:
  
  G4ModelCmdSetAuxPtsType(M* model, const G4String& placement)
    :G4ModelCmdApplyString<M>(model, placement, "setAuxPtsType")
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();

    cmd->SetGuidance("Set auxiliary points type.");
    cmd->SetCandidates("dots circles squares");
  }
  
protected:
  
  void Apply(const G4String& type) {
    G4Polymarker::MarkerType myType;
    
    if (type == "dots") {myType = G4Polymarker::dots;}
    else if (type == "circles") {myType = G4Polymarker::circles;}
    else if (type == "squares") {myType = G4Polymarker::squares;}
    else {
      std::ostringstream o;
      o << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetAuxPtsType::Apply",
	 "InvalidArgument", JustWarning, o.str().c_str());
      return;
    }
    
    G4VModelCommand<M>::Model()->SetAuxPtsType(myType);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points fill style command
template <typename M>
class G4ModelCmdSetStepPtsFillStyle : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsFillStyle(M* model, const G4String& placement)
    :G4ModelCmdApplyString<M>(model, placement, "setStepPtsFillStyle")
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set step fill style type.");
    cmd->SetCandidates("noFill hashed filled");
  }
  
protected:

  void Apply(const G4String& type) {
    G4VMarker::FillStyle myFillStyle;
    
    if (type == "noFill") {myFillStyle = G4VMarker::noFill;}
    else if (type == "hashed") {myFillStyle = G4VMarker::hashed;}
    else if (type == "filled") {myFillStyle = G4VMarker::filled;}
    else {
      std::ostringstream o;
      o << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetStepPtsFillStyle::Apply",
	 "InvalidArgument", JustWarning, o.str().c_str());
      return;
    }
    G4VModelCommand<M>::Model()->SetStepPtsFillStyle(myFillStyle);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points fill style command
template <typename M>
class G4ModelCmdSetAuxPtsFillStyle : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetAuxPtsFillStyle(M* model, const G4String& placement)
    :G4ModelCmdApplyString<M>(model, placement, "setAuxPtsFillStyle")
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set auxiliary fill style.");
    cmd->SetCandidates("noFill hashed filled");
  }
  
protected:

  void Apply(const G4String& type) {
    G4VMarker::FillStyle myFillStyle;
    
    if (type == "noFill") {myFillStyle = G4VMarker::noFill;}
    else if (type == "hashed") {myFillStyle = G4VMarker::hashed;}
    else if (type == "filled") {myFillStyle = G4VMarker::filled;}
    else {
      std::ostringstream o;
      o << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetAuxPtsFillStyle::Apply",
	 "InvalidArgument", JustWarning, o.str().c_str());
      return;
    }
    G4VModelCommand<M>::Model()->SetAuxPtsFillStyle(myFillStyle);
  }
  
};

////////////////////////////////////////////////////////////////////////
// SetLineColour command
template <typename M>
class G4ModelCmdSetLineColour : public G4ModelCmdApplyColour<M> {
  
public:

  G4ModelCmdSetLineColour(M* model, const G4String& placement)
    :G4ModelCmdApplyColour<M>(model, placement, "setLineColour") {}
  
protected:

  void Apply(const G4Colour& colour) {
    G4VModelCommand<M>::Model()->SetLineColour(colour);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Create context directory command
template <typename M>
class G4ModelCmdCreateContextDir : public G4UImessenger {
  
public:

  G4ModelCmdCreateContextDir(M* model, const G4String& placement) {
    G4String title = placement+"/"+model->Name()+"/";
    cmd = new G4UIdirectory(title);
    
    cmd->SetGuidance("Commands for default configuration");
  }

  virtual ~G4ModelCmdCreateContextDir() {
    delete cmd;
  }

protected:

  G4UIdirectory* cmd;
  
};
#endif      

