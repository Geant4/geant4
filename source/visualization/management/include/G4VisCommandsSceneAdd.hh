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
// $Id: G4VisCommandsSceneAdd.hh 104163 2017-05-15 06:52:42Z gcosmo $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENEADD_HH
#define G4VISCOMMANDSSCENEADD_HH

#include "G4VisCommandsScene.hh"

class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Timer.hh"

class G4VisCommandSceneAddArrow: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddArrow ();
  virtual ~G4VisCommandSceneAddArrow ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddArrow (const G4VisCommandSceneAddArrow&);
  G4VisCommandSceneAddArrow& operator = (const G4VisCommandSceneAddArrow&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddArrow2D: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddArrow2D ();
  virtual ~G4VisCommandSceneAddArrow2D ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddArrow2D (const G4VisCommandSceneAddArrow2D&);
  G4VisCommandSceneAddArrow2D& operator = (const G4VisCommandSceneAddArrow2D&);
  struct Arrow2D {
    Arrow2D(G4double x1, G4double y1,
	    G4double x2, G4double y2,
	    G4double width, const G4Colour& colour);
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4Polyline fShaftPolyline;
    G4Polyline fHeadPolyline;
    G4double fWidth;
    G4Colour fColour;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddAxes: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddAxes ();
  virtual ~G4VisCommandSceneAddAxes ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddAxes (const G4VisCommandSceneAddAxes&);
  G4VisCommandSceneAddAxes& operator = (const G4VisCommandSceneAddAxes&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddDate: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddDate ();
  virtual ~G4VisCommandSceneAddDate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddDate (const G4VisCommandSceneAddDate&);
  G4VisCommandSceneAddDate& operator = (const G4VisCommandSceneAddDate&);
  struct Date {
    Date
    (G4VisManager* vm, G4int size,
     G4double x, G4double y, G4Text::Layout layout,
     const G4String& date):
      fpVisManager(vm), fSize(size),
      fX(x), fY(y), fLayout(layout), fDate(date) {}
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4VisManager* fpVisManager;
    G4Timer fTimer;
    G4int fSize;
    G4double fX, fY;
    G4Text::Layout fLayout;
    G4String fDate;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddDigis: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddDigis ();
  virtual ~G4VisCommandSceneAddDigis ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddDigis (const G4VisCommandSceneAddDigis&);
  G4VisCommandSceneAddDigis& operator = (const G4VisCommandSceneAddDigis&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneAddEventID: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddEventID ();
  virtual ~G4VisCommandSceneAddEventID ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddEventID (const G4VisCommandSceneAddEventID&);
  G4VisCommandSceneAddEventID& operator = (const G4VisCommandSceneAddEventID&);
  struct EventID {
    EventID(G4VisManager* vm, G4int size,
	    G4double x, G4double y, G4Text::Layout layout):
      fpVisManager(vm), fSize(size),
      fX(x), fY(y), fLayout(layout) {}
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4VisManager* fpVisManager;
    G4int fSize;
    G4double fX, fY;
    G4Text::Layout fLayout;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddExtent: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddExtent ();
  virtual ~G4VisCommandSceneAddExtent ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddExtent (const G4VisCommandSceneAddExtent&);
  G4VisCommandSceneAddExtent& operator = (const G4VisCommandSceneAddExtent&);
  struct Extent {
    Extent(G4double xmin, G4double xmax,
           G4double ymin, G4double ymax,
           G4double zmin, G4double zmax);
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4VisExtent fExtent;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddFrame: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddFrame ();
  virtual ~G4VisCommandSceneAddFrame ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddFrame (const G4VisCommandSceneAddFrame&);
  G4VisCommandSceneAddFrame& operator = (const G4VisCommandSceneAddFrame&);
  struct Frame {
    Frame(G4double size, G4double width, const G4Colour& colour):
      fSize(size), fWidth(width), fColour(colour) {}
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4double fSize;
    G4double fWidth;
    G4Colour fColour;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddGPS: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddGPS ();
  virtual ~G4VisCommandSceneAddGPS ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddGPS (const G4VisCommandSceneAddGPS&);
  G4VisCommandSceneAddGPS& operator = (const G4VisCommandSceneAddGPS&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddGhosts: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddGhosts ();
  virtual ~G4VisCommandSceneAddGhosts ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddGhosts (const G4VisCommandSceneAddGhosts&);
  G4VisCommandSceneAddGhosts& operator =
  (const G4VisCommandSceneAddGhosts&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddHits: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddHits ();
  virtual ~G4VisCommandSceneAddHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddHits (const G4VisCommandSceneAddHits&);
  G4VisCommandSceneAddHits& operator = (const G4VisCommandSceneAddHits&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneAddLine: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLine ();
  virtual ~G4VisCommandSceneAddLine ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLine (const G4VisCommandSceneAddLine&);
  G4VisCommandSceneAddLine& operator = (const G4VisCommandSceneAddLine&);
  struct Line {
    Line(G4double x1, G4double y1, G4double z1,
	 G4double x2, G4double y2, G4double z2,
	 G4double width, const G4Colour& colour);
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4Polyline fPolyline;
    G4double fWidth;
    G4Colour fColour;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddLine2D: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLine2D ();
  virtual ~G4VisCommandSceneAddLine2D ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLine2D (const G4VisCommandSceneAddLine2D&);
  G4VisCommandSceneAddLine2D& operator = (const G4VisCommandSceneAddLine2D&);
  struct Line2D {
    Line2D(G4double x1, G4double y1,
	 G4double x2, G4double y2,
	 G4double width, const G4Colour& colour);
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4Polyline fPolyline;
    G4double fWidth;
    G4Colour fColour;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddLogicalVolume: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLogicalVolume ();
  virtual ~G4VisCommandSceneAddLogicalVolume ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLogicalVolume (const G4VisCommandSceneAddLogicalVolume&);
  G4VisCommandSceneAddLogicalVolume& operator =
  (const G4VisCommandSceneAddLogicalVolume&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddLogo: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLogo ();
  virtual ~G4VisCommandSceneAddLogo ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLogo (const G4VisCommandSceneAddLogo&);
  G4VisCommandSceneAddLogo& operator = (const G4VisCommandSceneAddLogo&);
  // Direction of outward-facing normal to front face of logo.
  enum Direction {X, minusX, Y, minusY, Z, minusZ};
  struct G4Logo {
    G4Logo(G4double height, const G4VisAttributes&);
    ~G4Logo();
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
  private:
    G4VisAttributes fVisAtts;
    G4Polyhedron *fpG, *fp4;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddLogo2D: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLogo2D ();
  virtual ~G4VisCommandSceneAddLogo2D ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLogo2D (const G4VisCommandSceneAddLogo2D&);
  G4VisCommandSceneAddLogo2D& operator = (const G4VisCommandSceneAddLogo2D&);
  struct Logo2D {
    Logo2D
    (G4VisManager* vm, G4int size,
     G4double x, G4double y, G4Text::Layout layout):
      fpVisManager(vm), fSize(size),
      fX(x), fY(y), fLayout(layout) {}
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4VisManager* fpVisManager;
    G4int fSize;
    G4double fX, fY;
    G4Text::Layout fLayout;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddMagneticField: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddMagneticField ();
  virtual ~G4VisCommandSceneAddMagneticField ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddMagneticField (const G4VisCommandSceneAddMagneticField&);
  G4VisCommandSceneAddMagneticField& operator = (const G4VisCommandSceneAddMagneticField&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddPSHits: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddPSHits ();
  virtual ~G4VisCommandSceneAddPSHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddPSHits (const G4VisCommandSceneAddPSHits&);
  G4VisCommandSceneAddPSHits& operator = (const G4VisCommandSceneAddPSHits&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddScale: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddScale ();
  virtual ~G4VisCommandSceneAddScale ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddScale (const G4VisCommandSceneAddScale&);
  G4VisCommandSceneAddScale& operator = (const G4VisCommandSceneAddScale&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddText: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddText ();
  virtual ~G4VisCommandSceneAddText ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddText (const G4VisCommandSceneAddText&);
  G4VisCommandSceneAddText& operator = (const G4VisCommandSceneAddText&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddText2D: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddText2D ();
  virtual ~G4VisCommandSceneAddText2D ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddText2D (const G4VisCommandSceneAddText2D&);
  G4VisCommandSceneAddText2D& operator = (const G4VisCommandSceneAddText2D&);
  struct G4Text2D {
    G4Text2D(const G4Text&);
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
  private:
    G4Text fText;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddTrajectories: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddTrajectories ();
  virtual ~G4VisCommandSceneAddTrajectories ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddTrajectories (const G4VisCommandSceneAddTrajectories&);
  G4VisCommandSceneAddTrajectories& operator =
  (const G4VisCommandSceneAddTrajectories&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddUserAction: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddUserAction ();
  virtual ~G4VisCommandSceneAddUserAction ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddUserAction (const G4VisCommandSceneAddUserAction&);
  G4VisCommandSceneAddUserAction& operator = (const G4VisCommandSceneAddUserAction&);
  enum ActionType {runDuration, endOfEvent, endOfRun};
  void AddVisAction(const G4String& name,G4VUserVisAction*,
		    G4Scene*,ActionType,G4VisManager::Verbosity);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddVolume: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddVolume ();
  virtual ~G4VisCommandSceneAddVolume ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddVolume (const G4VisCommandSceneAddVolume&);
  G4VisCommandSceneAddVolume& operator = (const G4VisCommandSceneAddVolume&);
  G4UIcommand* fpCommand;
};

#endif
