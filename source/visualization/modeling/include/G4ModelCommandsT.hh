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
// $Id: G4ModelCommandsT.hh 66373 2012-12-18 09:41:34Z gcosmo $
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

  G4ModelCmdSetStringColour(M* model, const G4String& placement, 
			    const G4String& cmdName="set")
    :G4ModelCmdApplyStringColour<M>(model, placement, cmdName) {}
  
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

  G4ModelCmdSetDefaultColour(M* model, const G4String& placement, 
			     const G4String& cmdName="setDefault")
    :G4ModelCmdApplyColour<M>(model, placement, cmdName) {}
  
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

  G4ModelCmdAddString(M* model, const G4String& placement,
		      const G4String& cmdName="add")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
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

  G4ModelCmdAddInt(M* model, const G4String& placement,
		   const G4String& cmdName="add")
  :G4ModelCmdApplyInteger<M>(model, placement, cmdName)
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
  
  G4ModelCmdInvert(M* model, const G4String& placement,
		   const G4String& cmdName="invert")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdActive(M* model, const G4String& placement,
		   const G4String& cmdName="active")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdVerbose(M* model, const G4String& placement,
		    const G4String& cmdName="verbose")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdReset(M* model, const G4String& placement,
			    const G4String& cmdName="reset") 
    :G4ModelCmdApplyNull<M>(model, placement, cmdName)
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

  G4ModelCmdSetAuxPtsColour(M* model, const G4String& placement,
			    const G4String& cmdName="setAuxPtsColour") 
    :G4ModelCmdApplyColour<M>(model, placement, cmdName) {}

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

  G4ModelCmdSetStepPtsColour(M* model, const G4String& placement,
			     const G4String& cmdName="setStepPtsColour")
    :G4ModelCmdApplyColour<M>(model, placement, cmdName) {}

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

  G4ModelCmdSetDrawLine(M* model, const G4String& placement,
			const G4String& cmdName="setDrawLine")    
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdSetLineVisible(M* model, const G4String& placement,
			   const G4String& cmdName="setLineVisible") 
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdSetDrawAuxPts(M* model, const G4String& placement,
			  const G4String& cmdName="setDrawAuxPts")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdSetAuxPtsVisible(M* model, const G4String& placement,
			     const G4String& cmdName="setAuxPtsVisible")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdSetDrawStepPts(M* model, const G4String& placement,
			   const G4String& cmdName="setDrawStepPts")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
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

  G4ModelCmdSetStepPtsVisible(M* model, const G4String& placement,
			      const G4String& cmdName="setStepPtsVisible")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Set step points visible command");
  }
  
protected:

   void Apply(const G4bool& myBool) {
     G4VModelCommand<M>::Model()->SetStepPtsVisible(myBool);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points size command
template <typename M>
class G4ModelCmdSetAuxPtsSize : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetAuxPtsSize(M* model, const G4String& placement,
			   const G4String& cmdName="setAuxPtsSize")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
  {
    G4ModelCmdApplyString<M>::Command()->SetGuidance("Set auxiliary points size command");
  }
  
protected:

   void Apply(const G4String& sizeString) {
     std::istringstream iss(sizeString);
     G4double size;
     G4String unit;
     iss >> size >> unit;
     if (G4VModelCommand<M>::Model()->GetAuxPtsSizeType() == G4VMarker::world)
       {
	 G4double myDouble = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(sizeString);
	 G4VModelCommand<M>::Model()->SetAuxPtsSize(myDouble);
       }
     else  // none or screen
       {
	 G4VModelCommand<M>::Model()->SetAuxPtsSize(size);
       }
   }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points size command
template <typename M>
class G4ModelCmdSetStepPtsSize : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsSize(M* model, const G4String& placement,
			   const G4String& cmdName="setStepPtsSize")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
  {
    G4ModelCmdApplyString<M>::Command()->SetGuidance("Set step points size command");
  }
  
protected:

   void Apply(const G4String& sizeString) {
     std::istringstream iss(sizeString);
     G4double size;
     G4String unit;
     iss >> size >> unit;
     if (G4VModelCommand<M>::Model()->GetStepPtsSizeType() == G4VMarker::world)
       {
	 G4double myDouble = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(sizeString);
	 G4VModelCommand<M>::Model()->SetStepPtsSize(myDouble);
       }
     else  // none or screen
       {
	 G4VModelCommand<M>::Model()->SetStepPtsSize(size);
       }
   }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points type command
template <typename M>
class G4ModelCmdSetStepPtsType : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsType(M* model, const G4String& placement,
			   const G4String& cmdName="setStepPtsType")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
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
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetStepPtsType::Apply",
	 "modeling0109", JustWarning, ed);
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
  
  G4ModelCmdSetAuxPtsType(M* model, const G4String& placement,
			  const G4String& cmdName="setAuxPtsType")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
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
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetAuxPtsType::Apply",
	 "modeling0110", JustWarning, ed);
      return;
    }
    
    G4VModelCommand<M>::Model()->SetAuxPtsType(myType);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points size type command
template <typename M>
class G4ModelCmdSetStepPtsSizeType : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsSizeType(M* model, const G4String& placement,
				const G4String& cmdName="setStepPtsSizeType")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set step size type.");
    cmd->SetCandidates("none world screen");
  }
  
protected:

  void Apply(const G4String& type) {
    G4VMarker::SizeType mySizeType;
    
    if (type == "none") {mySizeType = G4VMarker::none;}
    else if (type == "world") {mySizeType = G4VMarker::world;}
    else if (type == "screen") {mySizeType = G4VMarker::screen;}
    else {
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetStepPtsSizeType::Apply",
	 "modeling0111", JustWarning, ed);
      return;
    }
    G4VModelCommand<M>::Model()->SetStepPtsSizeType(mySizeType);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set auxiliary points size type command
template <typename M>
class G4ModelCmdSetAuxPtsSizeType : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetAuxPtsSizeType(M* model, const G4String& placement,
			       const G4String& cmdName="setAuxPtsSizeType")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set auxiliary size type.");
    cmd->SetCandidates("none world screen");
  }
  
protected:

  void Apply(const G4String& type) {
    G4VMarker::SizeType mySizeType;
    
    if (type == "none") {mySizeType = G4VMarker::none;}
    else if (type == "world") {mySizeType = G4VMarker::world;}
    else if (type == "screen") {mySizeType = G4VMarker::screen;}
    else {
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetAuxPtsSizeType::Apply",
	 "modeling0112", JustWarning, ed);
      return;
    }
    G4VModelCommand<M>::Model()->SetAuxPtsSizeType(mySizeType);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set step points fill style command
template <typename M>
class G4ModelCmdSetStepPtsFillStyle : public G4ModelCmdApplyString<M> {
  
public:

  G4ModelCmdSetStepPtsFillStyle(M* model, const G4String& placement,
				const G4String& cmdName="setStepPtsFillStyle")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
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
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetStepPtsFillStyle::Apply",
	 "modeling0113", JustWarning, ed);
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

  G4ModelCmdSetAuxPtsFillStyle(M* model, const G4String& placement,
			       const G4String& cmdName="setAuxPtsFillStyle")
    :G4ModelCmdApplyString<M>(model, placement, cmdName)
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
      G4ExceptionDescription ed;
      ed << "Invalid argument. See command guidance for options.";
      G4Exception
	("G4ModelCmdSetAuxPtsFillStyle::Apply",
	 "modeling0114", JustWarning, ed);
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

  G4ModelCmdSetLineColour(M* model, const G4String& placement,
			  const G4String& cmdName="""setLineColour")
    :G4ModelCmdApplyColour<M>(model, placement, cmdName){}
  
protected:

  void Apply(const G4Colour& colour) {
    G4VModelCommand<M>::Model()->SetLineColour(colour);
  }
  
};

////////////////////////////////////////////////////////////////////////
// Set time slice interval command
template <typename M>
class G4ModelCmdSetTimeSliceInterval : public G4ModelCmdApplyDoubleAndUnit<M> {

public:

  G4ModelCmdSetTimeSliceInterval(M* model, const G4String& placement,
				 const G4String& cmdName = "setTimeSliceInterval")
    :G4ModelCmdApplyDoubleAndUnit<M>(model, placement, cmdName)
  {
    G4UIcmdWithADoubleAndUnit* cmd = G4ModelCmdApplyDoubleAndUnit<M>::Command();
    cmd->SetGuidance
      ("Set time slice interval.  Give unit, e.g., \"0.1 ns\"");
    cmd->SetUnitCategory("Time");
  }

protected:

   void Apply(const G4double& myDouble) {
     G4VModelCommand<M>::Model()->SetTimeSliceInterval(myDouble);
  }

};

////////////////////////////////////////////////////////////////////////
// Create context directory command
template <typename M>
class G4ModelCmdCreateContextDir : public G4UImessenger {
  
public:

  G4ModelCmdCreateContextDir(M* model, const G4String& placement) 
  {
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

////////////////////////////////////////////////////////////////////////
// Draw command
template <typename M>
class G4ModelCmdDraw : public G4ModelCmdApplyBool<M> {

public: // With description

  G4ModelCmdDraw(M* model, const G4String& placement,
		 const G4String& cmdName="draw")
    :G4ModelCmdApplyBool<M>(model, placement, cmdName)
  {
    G4ModelCmdApplyBool<M>::Command()->SetGuidance("Draw command");
  }

  virtual ~G4ModelCmdDraw() {}

protected:

  virtual void Apply(const G4bool& newValue) {
    if (newValue) G4VModelCommand<M>::Model()->Draw();
  }

};

////////////////////////////////////////////////////////////////////////
// Set interval
template <typename M>
class G4ModelCmdAddInterval : public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdAddInterval(M* model, const G4String& placement, 
			const G4String& cmdName="addInterval")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set interval.");
  }
  
  virtual ~G4ModelCmdAddInterval() {}

protected:

  virtual void Apply(const G4String& param) {
    G4VModelCommand<M>::Model()->AddInterval(param);

  }
};

////////////////////////////////////////////////////////////////////////
// Set value
template <typename M>
class G4ModelCmdAddValue : public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdAddValue(M* model, const G4String& placement, 
		     const G4String& cmdName="addValue")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
  {
    G4UIcmdWithAString* cmd = G4ModelCmdApplyString<M>::Command();
    cmd->SetGuidance("Set value.");
  }
  
  virtual ~G4ModelCmdAddValue() {}
protected:

  virtual void Apply(const G4String& param) {
    G4VModelCommand<M>::Model()->AddValue(param);

  }
};

////////////////////////////////////////////////////////////////////////
// Set string command
template <typename M>
class G4ModelCmdSetString : public G4ModelCmdApplyString<M> {

public: // With description

  G4ModelCmdSetString(M* model, const G4String& placement, 
		      const G4String& cmdName="set")
    :G4ModelCmdApplyString<M>(model, placement, cmdName) 
  {
    G4ModelCmdApplyString<M>::Command()->SetGuidance("Set command");
  }

  virtual ~G4ModelCmdSetString() {}

protected:
  
  virtual void Apply(const G4String& newValue) {
    G4VModelCommand<M>::Model()->Set(newValue);
  }
  
};

#endif      
