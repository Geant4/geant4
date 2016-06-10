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
// $Id: G4VisCommandsGeometrySet.hh 66373 2012-12-18 09:41:34Z gcosmo $

// /vis/geometry commands - John Allison  31st January 2006

#ifndef G4VISCOMMANDSGEOMETRYSET_HH
#define G4VISCOMMANDSGEOMETRYSET_HH

#include "G4VisCommandsGeometry.hh"

class G4UIcommand;
class G4VisAttributes;

class G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VVisCommandGeometrySetFunction() {}
  virtual void operator()(G4VisAttributes*) const = 0;
};

class G4VisCommandGeometrySetColourFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetColourFunction() {}
  G4VisCommandGeometrySetColourFunction
  (const G4Colour& colour):
    fColour(colour) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetColour(fColour);}
private:
  const G4Colour& fColour;
};

class G4VisCommandGeometrySetDaughtersInvisibleFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetDaughtersInvisibleFunction() {}
  G4VisCommandGeometrySetDaughtersInvisibleFunction
  (G4bool daughtersInvisible):
    fDaughtersInvisible(daughtersInvisible) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetDaughtersInvisible(fDaughtersInvisible);}
private:
  G4bool fDaughtersInvisible;
};

class G4VisCommandGeometrySetForceAuxEdgeVisibleFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetForceAuxEdgeVisibleFunction() {}
  G4VisCommandGeometrySetForceAuxEdgeVisibleFunction
  (G4bool forceAuxEdgeVisible):
    fForceAuxEdgeVisible(forceAuxEdgeVisible) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetForceAuxEdgeVisible(fForceAuxEdgeVisible);}
private:
  G4bool fForceAuxEdgeVisible;
};

class G4VisCommandGeometrySetForceLineSegmentsPerCircleFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetForceLineSegmentsPerCircleFunction() {}
  G4VisCommandGeometrySetForceLineSegmentsPerCircleFunction
  (G4int lineSegmentsPerCircle):
    fLineSegmentsPerCircle(lineSegmentsPerCircle) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetForceLineSegmentsPerCircle(fLineSegmentsPerCircle);}
private:
  G4int fLineSegmentsPerCircle;
};

class G4VisCommandGeometrySetForceSolidFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetForceSolidFunction() {}
  G4VisCommandGeometrySetForceSolidFunction
  (G4bool forceSolid):
    fForceSolid(forceSolid) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetForceSolid(fForceSolid);}
private:
  G4bool fForceSolid;
};

class G4VisCommandGeometrySetForceWireframeFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetForceWireframeFunction() {}
  G4VisCommandGeometrySetForceWireframeFunction
  (G4bool forceWireframe):
    fForceWireframe(forceWireframe) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetForceWireframe(fForceWireframe);}
private:
  G4bool fForceWireframe;
};

class G4VisCommandGeometrySetLineStyleFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetLineStyleFunction() {}
  G4VisCommandGeometrySetLineStyleFunction
  (G4VisAttributes::LineStyle lineStyle):
    fLineStyle(lineStyle) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetLineStyle(fLineStyle);}
private:
  G4VisAttributes::LineStyle fLineStyle;
};

class G4VisCommandGeometrySetLineWidthFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetLineWidthFunction() {}
  G4VisCommandGeometrySetLineWidthFunction
  (G4double lineWidth):
    fLineWidth(lineWidth) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetLineWidth(fLineWidth);}
private:
  G4double fLineWidth;
};

class G4VisCommandGeometrySetVisibilityFunction:
  public G4VVisCommandGeometrySetFunction {
public:
  virtual ~G4VisCommandGeometrySetVisibilityFunction() {}
  G4VisCommandGeometrySetVisibilityFunction
  (G4bool visibility):
    fVisibility(visibility) {}
  void operator()
    (G4VisAttributes* visAtts) const
  {visAtts->SetVisibility(fVisibility);}
private:
  G4bool fVisibility;
};

class G4VVisCommandGeometrySet: public G4VVisCommandGeometry {
protected:
  void Set(G4String logVolName, const G4VVisCommandGeometrySetFunction&,
	   G4int requestedDepth);
  void SetLVVisAtts(G4LogicalVolume*, const G4VVisCommandGeometrySetFunction&,
		    G4int depth, G4int requestedDepth);
};

class G4VisCommandGeometrySetColour: public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetColour ();
  virtual ~G4VisCommandGeometrySetColour ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetColour (const G4VisCommandGeometrySetColour&);
  G4VisCommandGeometrySetColour& operator = (const G4VisCommandGeometrySetColour&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetDaughtersInvisible:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetDaughtersInvisible ();
  virtual ~G4VisCommandGeometrySetDaughtersInvisible ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetDaughtersInvisible
  (const G4VisCommandGeometrySetDaughtersInvisible&);
  G4VisCommandGeometrySetDaughtersInvisible& operator=
  (const G4VisCommandGeometrySetDaughtersInvisible&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetForceAuxEdgeVisible:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetForceAuxEdgeVisible ();
  virtual ~G4VisCommandGeometrySetForceAuxEdgeVisible ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetForceAuxEdgeVisible
  (const G4VisCommandGeometrySetForceAuxEdgeVisible&);
  G4VisCommandGeometrySetForceAuxEdgeVisible& operator=
  (const G4VisCommandGeometrySetForceAuxEdgeVisible&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetForceSolid:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetForceSolid ();
  virtual ~G4VisCommandGeometrySetForceSolid ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetForceSolid
  (const G4VisCommandGeometrySetForceSolid&);
  G4VisCommandGeometrySetForceSolid& operator=
  (const G4VisCommandGeometrySetForceSolid&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetForceLineSegmentsPerCircle:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetForceLineSegmentsPerCircle ();
  virtual ~G4VisCommandGeometrySetForceLineSegmentsPerCircle ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetForceLineSegmentsPerCircle
  (const G4VisCommandGeometrySetForceLineSegmentsPerCircle&);
  G4VisCommandGeometrySetForceLineSegmentsPerCircle& operator=
  (const G4VisCommandGeometrySetForceLineSegmentsPerCircle&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetForceWireframe:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetForceWireframe ();
  virtual ~G4VisCommandGeometrySetForceWireframe ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetForceWireframe
  (const G4VisCommandGeometrySetForceWireframe&);
  G4VisCommandGeometrySetForceWireframe& operator=
  (const G4VisCommandGeometrySetForceWireframe&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetLineStyle:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetLineStyle ();
  virtual ~G4VisCommandGeometrySetLineStyle ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetLineStyle
  (const G4VisCommandGeometrySetLineStyle&);
  G4VisCommandGeometrySetLineStyle& operator=
  (const G4VisCommandGeometrySetLineStyle&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetLineWidth:
  public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetLineWidth ();
  virtual ~G4VisCommandGeometrySetLineWidth ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometrySetLineWidth
  (const G4VisCommandGeometrySetLineWidth&);
  G4VisCommandGeometrySetLineWidth& operator=
  (const G4VisCommandGeometrySetLineWidth&);
  G4UIcommand* fpCommand;
};

class G4VisCommandGeometrySetVisibility: public G4VVisCommandGeometrySet {
public:
  G4VisCommandGeometrySetVisibility ();
  virtual ~G4VisCommandGeometrySetVisibility ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
  void SetNewValueOnLV (G4LogicalVolume* pLV, G4int,G4bool);
private:
  G4VisCommandGeometrySetVisibility (const G4VisCommandGeometrySetVisibility&);
  G4VisCommandGeometrySetVisibility& operator = (const G4VisCommandGeometrySetVisibility&);
  G4UIcommand* fpCommand;
};

#endif
