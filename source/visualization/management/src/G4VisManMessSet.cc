// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessSet.cc,v 1.8 2001-02-01 17:35:43 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Set sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

void G4VisManMessenger::AddCommandSet () {
  
  G4UIcommand* command;
  G4UIparameter* param;

  ///////////////////////////////////////////  /vis~/set/cull_by_density  ////
  //set \hline
  //set /vis~/set/cull\_by\_density & choice, density &
  //set Cull (i.e., do not Draw) geometry objects with density less
  //set than specified --- on/off and density (g / cm3). \\%
  command = new G4UIcommand ("/vis~/set/cull_by_density", this);
  command -> SetGuidance
    (
     "Cull (i.e., do not Draw) geometry objects with density less"
     "\nthan specified - on/off and density (g / cm3)."
     );
  param   =  new G4UIparameter ("Selector", 'c', true);
  param   -> SetDefaultValue  ("?");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("density", 'd', true);
  param   -> SetDefaultValue  (0.01);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/drawing_style  ////
  //set \hline
  //set /vis~/set/drawing\_style & choice &
  //set Style of drawing (wireframe, etc.). \\%
  command = new G4UIcommand ("/vis~/set/drawing_style", this);
  command -> SetGuidance
    (
     "Style of drawing (wireframe, etc.)."
     );
  param   =  new G4UIparameter ("Style selector", 'i', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/marker_choices  ////
  //set \hline
  //set /vis~/set/marker\_choices & choice &
  //set Choices for marker drawing (hidden by surfaces, etc.). \\%
  command = new G4UIcommand ("/vis~/set/marker_choices", this);
  command -> SetGuidance
    (
     "Chioces for marker drawing (hidden by surfaces, etc.)."
     );
  param   =  new G4UIparameter ("Style selector", 'i', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/rep_style  ////
  //set \hline
  //set /vis~/set/rep\_style & choice &
  //set Style of graphics representation for geometrical volumes
  //set (polyhedron, NURBS,...). \\%
  command = new G4UIcommand ("/vis~/set/rep_style", this);
  command -> SetGuidance
    (
     "Style of graphics representation for geometrical volumes "
     "(polyhedron, NURBS,...)."
     );
  param   =  new G4UIparameter ("Style selector", 'i', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/scene  ////
  //set \hline
  //set /vis~/set/scene & choice &
  //set Several ways of setting the scene.
  //set For example, choosing a physical volume. \\%
  command = new G4UIcommand ("/vis~/set/scene", this);
  command -> SetGuidance
    (
     "Several ways of setting the scene."
     "  For example, choosing a physical volume."
     );
  param   =  new G4UIparameter ("Selector", 'i', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/section_plane   ////
  //set \hline
  //set /vis~/set/section & choice, plane &
  //set Set plane for drawing section (DCUT).  Specify plane by
  //set x y z units nx ny nz, e.g., for a y-z plane at x = 1 cm:
  //set /vis~/set/section\_plane on 1 0 0 cm 1 0 0 \\%
  command = new G4UIcommand ("/vis~/set/section_plane", this);
  command -> SetGuidance
    (
     "Set plane for drawing section (DCUT).  Specify plane by"
     "\nx y z units nx ny nz, e.g., for a y-z plane at x = 1 cm:"
     "\n/vis~/set/section_plane on 1 0 0 cm 1 0 0"
     );
  param   =  new G4UIparameter ("Selector", 'c', true);
  param   -> SetDefaultValue  ("?");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("x", 'd', true);
  param   -> SetDefaultValue  (0);
  param   -> SetGuidance      ("Coordinate of point on the plane.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("y", 'd', true);
  param   -> SetDefaultValue  (0);
  param   -> SetGuidance      ("Coordinate of point on the plane.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("z", 'd', true);
  param   -> SetDefaultValue  (0);
  param   -> SetGuidance      ("Coordinate of point on the plane.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("units", 's', true);
  param   -> SetDefaultValue  ("mm");
  param   -> SetGuidance      ("mm, cm, or m.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("nx", 'd', true);
  param   -> SetDefaultValue  (1);
  param   -> SetGuidance      ("Component of plane normal.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("ny", 'd', true);
  param   -> SetDefaultValue  (0);
  param   -> SetGuidance      ("Component of plane normal.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("nz", 'd', true);
  param   -> SetDefaultValue  (0);
  param   -> SetGuidance      ("Component of plane normal.");
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/sides  ////
  //set \hline
  //set /vis~/set/sides & n &
  //set Number of sides per circle in polygon approximation. \\%
  command = new G4UIcommand ("/vis~/set/sides", this);
  command -> SetGuidance
    (
     "Number of sides per circle in polygon approximation."
     );
  param   =  new G4UIparameter ("No. of sides", 'i', true);
  param   -> SetDefaultValue  (24);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/verbose  //////
  //set \hline
  //set /vis~/set/verbose & integer & 
  //set Controls amount of printing (0 = quiet). \\%
  command = new G4UIcommand ("/vis~/set/verbose", this);
  command -> SetGuidance
    (
     "Controls amount of printing (0 = quiet)."
     );
  param   =  new G4UIparameter ("Verbosity level", 'i', true);
  param   -> SetDefaultValue (-1);
  param   -> SetGuidance ("0: quiet, >0: verbose");
  command -> SetParameter (param);
  fCommandList.append (command);

  ///////////////////////////////////////////  /vis~/set/view  /////////
  //set  \hline
  //set /vis~/set/view & choice &
  //set Make this view current. \\%
  command = new G4UIcommand ("/vis~/set/view", this);
  command -> SetGuidance
    (
     "Make this view current."
     );
  param   =  new G4UIparameter ("View selector", 'i', true);
  param   -> SetDefaultValue (-1);
  command -> SetParameter (param);
  fCommandList.append (command);
}

void G4VisManMessenger::DoCommandSet (const G4String& commandPath,
				      G4String& newValues) {

  ///////////////////////////////////////////  /vis~/set/cull_by_density  ////
  if (commandPath == "/vis~/set/cull_by_density") {
    G4String choice;
    G4double density;  // Units in this section are g / cm3 - WARNING!!!
    const char* t = newValues;
    G4std::istrstream is ((char*)t); is >> choice >> density;
    G4int iSelector = -1;
    if (choice.compareTo ("off",G4String::ignoreCase) == 0) iSelector = 0;
    if (choice.compareTo ("on",G4String::ignoreCase) == 0) iSelector = 1;
    if (iSelector < 0) {
      G4cout << "Choice not recognised (on/off)." << G4endl;
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      G4cout << "Culling by density is currently: ";
      if (getVP.IsDensityCulling ()) G4cout << "on";
      else                             G4cout << "off";
      G4cout << "\nDensity cut is currently: "
	   << getVP.GetVisibleDensity () * cm3 / g << " g/cm3.";
      G4cout << G4endl;
    }
    else {
      G4ViewParameters& setVP = fpVMan -> SetCurrentViewParameters ();
      switch (iSelector) {
      case 0:
	setVP.SetDensityCulling (false);
	break;
      case 1:
	setVP.SetDensityCulling (true);
	setVP.SetVisibleDensity (density * g / cm3);
	break;
      default: G4cout << "Choice not recognised (on/off).\n"; break;
      }
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      G4cout << "Culling by density is now: ";
      if (getVP.IsDensityCulling ()) G4cout << "on";
      else                             G4cout << "off";
      G4cout << "\nDensity cut is now: "
	   << getVP.GetVisibleDensity () * cm3 / g << " g/cm3.";
      G4cout << G4endl;
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (getVP);
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ///////////////////////////////////////////  /vis~/set/drawing_style  ////
  if (commandPath == "/vis~/set/drawing_style") {
    G4int iStyle;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> iStyle;
    if (iStyle < 0 || iStyle > 3) {
      G4cout << "Available drawing styles:";
      G4cout << "\n  0) wireframe";
      G4cout << "\n  1) hidden line removal";
      G4cout << "\n  2) solid (hidden surface removal)";
      G4cout << "\n  3) solids with edges (hidden line, hidden surface removal)";
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      //if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "\nDrawing style is currently"
	     << getVP.GetDrawingStyle ();
	G4cout << "\nCulling is currently ";
	if (getVP.IsCulling ()) G4cout << "on";
	else                    G4cout << "off";
	G4cout << "\nCulling of invisible objects is currently ";
	if (getVP.IsCullingInvisible ()) G4cout << "on";
	else                             G4cout << "off";
	G4cout << "\nDensity culling is currently ";
	if (getVP.IsDensityCulling ()) {
	  G4cout << "on - invisible if density less than "
	       << getVP.GetVisibleDensity () / (1. * g / cm3) << " g cm^-3";
	}
	else G4cout << "off";
	G4cout << "\nCovered daughters are ";
	if (getVP.IsCullingCovered ()) {
	  G4cout << "not";
	}
	G4cout << " drawn.";
	G4cout << "\nTo override culling, /vis~/set/culling off."
	  "\nTo set/reset culling of invisible objects,"
	  "\n/vis~/set/cull_invisible_objects on/off."
	  "\nTo reset density culling, /vis~/set/cull_by_density off."
	  "\nTo set density culling, e.g., /vis~/set/cull_by_density on 0.01."
	     << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      //}
      G4cout << "\nChoose by specifying integer parameter.";
      G4cout << G4endl;
    }
    else {
      G4ViewParameters& setVP = fpVMan -> SetCurrentViewParameters ();
      switch (iStyle) {
      default:
      case 0:
	setVP.SetDrawingStyle (G4ViewParameters::wireframe);
	break;
      case 1:
	setVP.SetDrawingStyle (G4ViewParameters::hlr);
	break;
      case 2:
	setVP.SetDrawingStyle (G4ViewParameters::hsr);
	break;
      case 3:
	setVP.SetDrawingStyle (G4ViewParameters::hlhsr);
	break;
      }
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      //if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Drawing style changed to "
	     << getVP.GetDrawingStyle ();
	G4cout << "\nCulling is now ";
	if (getVP.IsCulling ()) G4cout << "on";
	else                    G4cout << "off";
	G4cout << "\nCulling of invisible objects is now ";
	if (getVP.IsCullingInvisible ()) G4cout << "on";
	else                             G4cout << "off";
	G4cout << "\nDensity culling is now ";
	if (getVP.IsDensityCulling ()) {
	  G4cout << "on - invisible if density less than "
	       << getVP.GetVisibleDensity () / (1. * g / cm3) << " g cm^-3";
	}
	else G4cout << "off";
	G4cout << "\nCovered daughters are ";
	if (getVP.IsCullingCovered ()) {
	  G4cout << "not";
	}
	G4cout << " drawn.";
	G4cout << "\nTo override culling, /vis~/set/culling off."
	  "\nTo set/reset culling of invisible objects,"
	  "\n/vis~/set/cull_invisible_objects on/off."
	  "\nTo reset density culling, /vis~/set/cull_by_density off."
	  "\nTo set density culling, e.g., /vis~/set/cull_by_density on 0.01."
	     << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      //}
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (getVP);
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ////////////////////////////////////////  /vis~/set/marker_choices  ////
  if (commandPath == "/vis~/set/marker_choices") {
    G4int iChoice;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> iChoice;
    if (iChoice < 0 || iChoice > 1) {
      G4cout << "Available marker choices:";
      G4cout << "\n  0) not hidden";
      G4cout << "\n  1) hidden by surfaces";
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      //if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "\nMarker choice is currently: ";
	if (getVP.IsMarkerNotHidden ()) G4cout << "not";
	G4cout << " hidden by surfaces.";
      //}
      G4cout << "\nChoose by specifying integer parameter.";
      G4cout << G4endl;
    }
    else {
      G4ViewParameters& setVP = fpVMan -> SetCurrentViewParameters ();
      switch (iChoice) {
      default:
      case 0:
	setVP.SetMarkerNotHidden ();
	break;
      case 1:
	setVP.SetMarkerHidden ();
	break;
      }
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      //if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Marker choice is now: ";
	if (getVP.IsMarkerNotHidden ()) G4cout << "not";
	G4cout << " hidden by surfaces." << G4endl;
      //}
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (getVP);
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }


  ///////////////////////////////////////////  /vis~/set/section_plane  ////
  if (commandPath == "/vis~/set/section_plane") {
    G4String choice, unit;
    G4double x, y, z, nx, ny, nz;
    const char* t = newValues;
    G4std::istrstream is ((char*)t); is >> choice >> x >> y >> z >> unit
				  >> nx >> ny >> nz;
    G4int iSelector = -1;
    if (choice.compareTo ("off",G4String::ignoreCase) == 0) iSelector = 0;
    if (choice.compareTo ("on",G4String::ignoreCase) == 0) iSelector = 1;
    if (iSelector < 0) {
      G4cout << "Choice not recognised (on/off)." << G4endl;
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      G4cout << "Section drawing is currently: ";
      if (getVP.IsSection ()) G4cout << "on";
      else                    G4cout << "off";
      G4cout << "\nSection plane is currently: "
	   << getVP.GetSectionPlane ();
      G4cout << G4endl;
    }
    else {
      G4ViewParameters& setVP = fpVMan -> SetCurrentViewParameters ();
      G4double F;
      switch (iSelector) {
      case 0:
	setVP.UnsetSectionPlane ();
	break;
      case 1:
	F = G4UnitDefinition::GetValueOf (unit);
	x *= F; y *= F; z *= F;
	setVP.SetSectionPlane (G4Plane3D (G4Normal3D (nx, ny, nz),
					  G4Point3D (x, y, z)));
	setVP.SetViewpointDirection (G4Normal3D (nx, ny, nz));
	break;
      default: G4cout << "Choice not recognised (on/off).\n"; break;
      }
      const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
      G4cout << "Section drawing is now: ";
      if (getVP.IsSection ()) G4cout << "on";
      else                    G4cout << "off";
      G4cout << "\nSection plane is now: "
	   << getVP.GetSectionPlane ();
      G4cout << G4endl;
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (getVP);
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ///////////////////////////////////////////  /vis~/set/sides  ////
  if (commandPath == "/vis~/set/sides") {
    G4int nSides;
    const char* t = newValues;
    G4std::istrstream is ((char*)t); is >> nSides;
    G4cout << "Number of sides per circle in polygon approximation is "
	 << nSides << G4endl;
    fpVMan -> SetCurrentViewParameters ().SetNoOfSides (nSides);
    const G4ViewParameters& getVP = fpVMan -> GetCurrentViewParameters ();
    G4VViewer* pView = fpVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (getVP);
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }

  ///////////////////////////////////////////  /vis~/set/rep_style  ////
  if (commandPath == "/vis~/set/rep_style") {
    G4int iStyle;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> iStyle;
    if (iStyle < 0 || iStyle > 1) {
      G4cout << "Available representation styles:";
      G4cout << "\n  0) polyhedron";
      G4cout << "\n  1) NURBS";
      G4cout << "\nChoose by specifying integer parameter.";
      G4cout << G4endl;
    }
    else {
      G4ViewParameters::RepStyle repStyle;
      switch (iStyle) {
      default:
      case 0: repStyle = G4ViewParameters::polyhedron; break;
      case 1: repStyle = G4ViewParameters::nurbs; break;
      }
      fpVMan -> SetCurrentViewParameters ().SetRepStyle (repStyle);
      if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Representation style changed to " << repStyle << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      }
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ///////////////////////////////////  /vis~/set/scene  ////
  if (commandPath == "/vis~/set/scene") {
    static G4int iState = 0;
    static G4int nOptions = 2;
    G4int iSelector;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> iSelector;
    switch (iState) {
    case 0:
      if (iSelector < 0 || iSelector > nOptions) {
	G4cout << "Available options:";
	G4cout << "\n  0) reset selection algorithm";
	G4cout << "\n  1) select physical volume";
	G4cout << "\nChoose by specifying integer parameter.";
	G4cout << G4endl;
      }
      else {
	switch (iSelector) {
	default:
	case 0:
	  G4cout << "Resetting from state " << iState << G4endl;
	  iState = 0; nOptions = 2;
	  break;
	case 1:
	  G4cout << "Option " << iSelector
	       << " selected from state " << iState << G4endl;
	  iState = 1; nOptions = 1;
	  break;
	}
      }
      break;
    case 1:
      if (iSelector < 0 || iSelector > nOptions) {
	G4cout << "Available options:";
	G4cout << "\n  0) reset selection algorithm";
	G4cout << "\nChoose by specifying integer parameter.";
	G4cout << G4endl;
      }
      else {
	switch (iSelector) {
	default:
	case 0:
	  G4cout << "Resetting from state " << iState << G4endl;
	  iState = 0; nOptions = 2;
	  break;
	}
      }
      break;
    }
    /****************************************************
    G4cin.tie (&G4cout);  // Tie streams together so output is
                      // flushed before any input.
    G4cout << "Temporary algorithm..." << G4endl;

    G4VPhysicalVolume* pVPV =
      fpVMan -> GetCurrentScene ().GetPhysicalVolume ();
    G4VPhysicalVolume* pNewVPV = 0;

    if (pVPV) {
      G4cout << "Current volume is " << pVPV -> GetName () << G4endl;

      G4VPhysicalVolume* pMother = pVPV -> GetMother ();
      if (pMother) {
	char c;
	do {
	  G4cout << "Its mother is " << pMother -> GetName ()
	       << ".  Select? [y/n]: ";
	  c = G4cin.get ();
	  if (c != '\n') while (G4cin.get () != '\n');  // Read on to newline.
	} while (c != 'y' && c != 'Y' && 
		 c != 'n' && c != 'N');
	if (c == 'y' || c == 'Y') {
	  pNewVPV = pMother;
	}
      }

      if (!pNewVPV) {
	G4int nDaughters = pVPV -> GetLogicalVolume () -> GetNoDaughters ();
	G4int iDaughter;
	if (nDaughters) {
	  do {
	    G4cout << "It has " << nDaughters << " daughters.  Select [0-"
		 << nDaughters - 1 << "] (<0 to abandon): ";
	    G4cin >> iDaughter;
	    while (G4cin.get () != '\n');  // Read on to newline.
	  } while (iDaughter >= nDaughters);
	  if (iDaughter >= 0) {
	    G4VPhysicalVolume* pDaughter =
	      pVPV -> GetLogicalVolume () -> GetDaughter (iDaughter);
	    pNewVPV = pDaughter;
	  }
	}
      }

      if (pNewVPV) {
	fpVMan -> AddToCurrentSceneData (pNewVPV);
	fpVMan -> Clear ();  // Comprehensive clear.
	fpVMan -> Draw ();
	fpVMan -> Show ();
      }
      else {
	if (fpVMan -> GetVerboseLevel () > 0) {
	  G4cout << "Nothing selected; nothing changed." << G4endl;
	}
      }
    }
    else {
      G4cerr << "No physical volume available." << G4endl;
    }
    if (fpVMan -> GetVerboseLevel () > 1) {
      fpVMan -> PrintCurrentScene ();
    }
    *********************************************************/
  }

  ///////////////////////////////////////////  /vis~/set/verbose  //////
  if (commandPath == "/vis~/set/verbose") {
    G4int vLevel;
    const char* t = newValues;
    G4std::istrstream is ((char*)t); is >> vLevel;
    if (vLevel < 0 ) {
      G4cout << "Verbosity is " << fpVMan -> GetVerboseLevel () << G4endl;
    }
    else {
      fpVMan -> SetVerboseLevel (vLevel);
      G4cout << "Verbosity set to " << fpVMan -> GetVerboseLevel () << G4endl;
    }
  }

  ///////////////////////////////////////////  /vis~/set/view  /////////
  if (commandPath == "/vis~/set/view") {
    // Make List of available views.
    G4RWTPtrOrderedVector<G4VViewer> vList;
    const G4SceneHandlerList& gml = fpVMan -> GetAvailableSceneHandlers ();
    G4int nViewTotal = 0;
    G4int iGM, nScenes = gml.entries ();
    for (iGM = 0; iGM < nScenes; iGM++) {
      const G4ViewerList& views = gml [iGM] -> GetViewerList ();
      int nViews = views.entries ();
      for (int iView = 0; iView < nViews; iView++) {
	G4VViewer* pView = views [iView];
	vList.append (pView);
	nViewTotal++;
      }
    }
    if (nViewTotal == 0) {
      G4cerr << "No views to select." << G4endl;
    }
    else {
      // Select view and make it current.
      G4int iSelect;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> iSelect;
      if (iSelect < 0 || iSelect >= nViewTotal) {
	G4cout << "Available views:";
	for (int iView = 0; iView < nViewTotal; iView++) {
	  const G4VViewer* pView = vList [iView];
	  const G4VSceneHandler* pScene = pView -> GetSceneHandler ();
	  const G4VGraphicsSystem* pSystem = pScene -> GetGraphicsSystem ();
	  G4cout << "\n  " << iView << ") " << pView -> GetName ();
	}
	G4cout << "\nChoose by specifying integer parameter.";
	G4cout << G4endl;
      }
      else {
	G4VViewer* pView = vList[iSelect];
	fpVMan -> SetCurrentViewer (pView);
	if (fpVMan -> GetVerboseLevel () > 0) {
	  G4cout << "Current view now " << pView -> GetName () << G4endl;
	  if (fpVMan -> GetVerboseLevel () > 1) {
	    fpVMan -> PrintCurrentView ();
	  }
	}
      }
    }
  }
}
