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
//
// $Id: G4VisCommandsViewer.hh 102234 2017-01-13 15:59:41Z gcosmo $

// /vis/viewer commands - John Allison  25th October 1998

#ifndef G4VISCOMMANDSVIEWER_HH
#define G4VISCOMMANDSVIEWER_HH

#include "G4VVisCommand.hh"

class G4VViewer;
class G4ViewParameters;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;

class G4VVisCommandViewer: public G4VVisCommand {
public:
  G4VVisCommandViewer ();
  virtual ~G4VVisCommandViewer ();
protected:
  void SetViewParameters(G4VViewer*, const G4ViewParameters&);
  void RefreshIfRequired(G4VViewer*);
private:
  G4VVisCommandViewer (const G4VVisCommandViewer&);
  G4VVisCommandViewer& operator = (const G4VVisCommandViewer&);
};

class G4VisCommandViewerAddCutawayPlane: public G4VVisCommandViewer {
public:
  G4VisCommandViewerAddCutawayPlane ();
  virtual ~G4VisCommandViewerAddCutawayPlane ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerAddCutawayPlane (const G4VisCommandViewerAddCutawayPlane&);
  G4VisCommandViewerAddCutawayPlane& operator = (const G4VisCommandViewerAddCutawayPlane&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerChangeCutawayPlane: public G4VVisCommandViewer {
public:
  G4VisCommandViewerChangeCutawayPlane ();
  virtual ~G4VisCommandViewerChangeCutawayPlane ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerChangeCutawayPlane (const G4VisCommandViewerChangeCutawayPlane&);
  G4VisCommandViewerChangeCutawayPlane& operator = (const G4VisCommandViewerChangeCutawayPlane&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerClear: public G4VVisCommandViewer {
public:
  G4VisCommandViewerClear ();
  virtual ~G4VisCommandViewerClear ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerClear (const G4VisCommandViewerClear&);
  G4VisCommandViewerClear& operator = (const G4VisCommandViewerClear&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerClearCutawayPlanes: public G4VVisCommandViewer {
public:
  G4VisCommandViewerClearCutawayPlanes ();
  virtual ~G4VisCommandViewerClearCutawayPlanes ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerClearCutawayPlanes (const G4VisCommandViewerClearCutawayPlanes&);
  G4VisCommandViewerClearCutawayPlanes& operator = (const G4VisCommandViewerClearCutawayPlanes&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandViewerClearTransients: public G4VVisCommandViewer {
public:
  G4VisCommandViewerClearTransients ();
  virtual ~G4VisCommandViewerClearTransients ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerClearTransients (const G4VisCommandViewerClearTransients&);
  G4VisCommandViewerClearTransients& operator =
  (const G4VisCommandViewerClearTransients&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerClearVisAttributesModifiers: public G4VVisCommandViewer {
public:
  G4VisCommandViewerClearVisAttributesModifiers ();
  virtual ~G4VisCommandViewerClearVisAttributesModifiers ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerClearVisAttributesModifiers (const G4VisCommandViewerClearVisAttributesModifiers&);
  G4VisCommandViewerClearVisAttributesModifiers& operator = (const G4VisCommandViewerClearVisAttributesModifiers&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandViewerClone: public G4VVisCommandViewer {
public:
  G4VisCommandViewerClone ();
  virtual ~G4VisCommandViewerClone ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerClone (const G4VisCommandViewerClone&);
  G4VisCommandViewerClone& operator =
  (const G4VisCommandViewerClone&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerCopyViewFrom: public G4VVisCommandViewer {
public:
  G4VisCommandViewerCopyViewFrom ();
  virtual ~G4VisCommandViewerCopyViewFrom ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerCopyViewFrom (const G4VisCommandViewerCopyViewFrom&);
  G4VisCommandViewerCopyViewFrom& operator =
  (const G4VisCommandViewerCopyViewFrom&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerCreate: public G4VVisCommandViewer {
public:
  G4VisCommandViewerCreate ();
  virtual ~G4VisCommandViewerCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerCreate (const G4VisCommandViewerCreate&);
  G4VisCommandViewerCreate& operator = (const G4VisCommandViewerCreate&);
  G4String NextName ();
  G4UIcommand* fpCommand;
  G4int fId;
};

class G4VisCommandViewerDolly: public G4VVisCommandViewer {
public:
  G4VisCommandViewerDolly ();
  virtual ~G4VisCommandViewerDolly ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerDolly (const G4VisCommandViewerDolly&);
  G4VisCommandViewerDolly& operator = (const G4VisCommandViewerDolly&);
  G4UIcmdWithADoubleAndUnit* fpCommandDolly;
  G4UIcmdWithADoubleAndUnit* fpCommandDollyTo;
  G4double fDollyIncrement;
  G4double fDollyTo;
};

class G4VisCommandViewerFlush: public G4VVisCommandViewer {
public:
  G4VisCommandViewerFlush ();
  virtual ~G4VisCommandViewerFlush ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerFlush (const G4VisCommandViewerFlush&);
  G4VisCommandViewerFlush& operator = (const G4VisCommandViewerFlush&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerInterpolate: public G4VVisCommandViewer {
public:
  G4VisCommandViewerInterpolate ();
  virtual ~G4VisCommandViewerInterpolate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerInterpolate (const G4VisCommandViewerInterpolate&);
  G4VisCommandViewerInterpolate& operator = (const G4VisCommandViewerInterpolate&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerList: public G4VVisCommandViewer {
public:
  G4VisCommandViewerList ();
  virtual ~G4VisCommandViewerList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerList (const G4VisCommandViewerList&);
  G4VisCommandViewerList& operator = (const G4VisCommandViewerList&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerPan: public G4VVisCommandViewer {
public:
  G4VisCommandViewerPan ();
  virtual ~G4VisCommandViewerPan ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerPan (const G4VisCommandViewerPan&);
  G4VisCommandViewerPan& operator = (const G4VisCommandViewerPan&);
  G4UIcommand* fpCommandPan;
  G4UIcommand* fpCommandPanTo;
  G4double fPanIncrementRight, fPanIncrementUp;
  G4double fPanToRight, fPanToUp;
};

class G4VisCommandViewerReset: public G4VVisCommandViewer {
public:
  G4VisCommandViewerReset ();
  virtual ~G4VisCommandViewerReset ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerReset (const G4VisCommandViewerReset&);
  G4VisCommandViewerReset& operator = (const G4VisCommandViewerReset&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerRefresh: public G4VVisCommandViewer {
public:
  G4VisCommandViewerRefresh ();
  virtual ~G4VisCommandViewerRefresh ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerRefresh (const G4VisCommandViewerRefresh&);
  G4VisCommandViewerRefresh& operator = (const G4VisCommandViewerRefresh&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerRebuild: public G4VVisCommandViewer {
public:
  G4VisCommandViewerRebuild ();
  virtual ~G4VisCommandViewerRebuild ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerRebuild (const G4VisCommandViewerRebuild&);
  G4VisCommandViewerRebuild& operator = (const G4VisCommandViewerRebuild&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerSave: public G4VVisCommandViewer {
public:
  G4VisCommandViewerSave ();
  virtual ~G4VisCommandViewerSave ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerSave (const G4VisCommandViewerSave&);
  G4VisCommandViewerSave& operator = (const G4VisCommandViewerSave&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerScale: public G4VVisCommandViewer {
public:
  G4VisCommandViewerScale ();
  virtual ~G4VisCommandViewerScale ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerScale (const G4VisCommandViewerScale&);
  G4VisCommandViewerScale& operator = (const G4VisCommandViewerScale&);
  G4UIcmdWith3Vector* fpCommandScale;
  G4UIcmdWith3Vector* fpCommandScaleTo;
  G4Vector3D fScaleMultiplier;
  G4Vector3D fScaleTo;
};

class G4VisCommandViewerSelect: public G4VVisCommandViewer {
public:
  G4VisCommandViewerSelect ();
  virtual ~G4VisCommandViewerSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerSelect (const G4VisCommandViewerSelect&);
  G4VisCommandViewerSelect& operator = (const G4VisCommandViewerSelect&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerUpdate: public G4VVisCommandViewer {
public:
  G4VisCommandViewerUpdate ();
  virtual ~G4VisCommandViewerUpdate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerUpdate (const G4VisCommandViewerUpdate&);
  G4VisCommandViewerUpdate& operator = (const G4VisCommandViewerUpdate&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerZoom: public G4VVisCommandViewer {
public:
  G4VisCommandViewerZoom ();
  virtual ~G4VisCommandViewerZoom ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerZoom (const G4VisCommandViewerZoom&);
  G4VisCommandViewerZoom& operator = (const G4VisCommandViewerZoom&);
  G4UIcmdWithADouble* fpCommandZoom;
  G4UIcmdWithADouble* fpCommandZoomTo;
  G4double fZoomMultiplier;
  G4double fZoomTo;
};

#endif
