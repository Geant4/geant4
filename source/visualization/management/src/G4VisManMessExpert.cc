// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessExpert.cc,v 1.6 2001-02-23 15:43:30 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Expert sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4Polyline.hh"
#include "G4Square.hh"
#include "G4Circle.hh"
#include "G4Text.hh"
#include "G4Polymarker.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "Randomize.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

void G4VisManMessenger::AddCommandExpert () {

  G4UIcommand* command;
  G4UIparameter* param;

  ////////////////////////////////////////////////////  /vis~/expert/...  ////
  //vis \hline
  //vis /vis~/expert/ &&
  //vis ...menu of expert commands. \\%
  command = new G4UIcommand ("/vis~/expert/", this);
  command -> SetGuidance
    (
     "...menu of expert commands."
     );
  fCommandList.push_back (command);

  ////////////////////////////////////////  /vis~/expert/draw_circle  ////
  //expert \hline
  //expert /vis~/expert/draw\_circle & $x$, $y$, $z$, $size$ &
  //expert Specify coordinates and world size. \\%
  command = new G4UIcommand ("/vis~/expert/draw_circle", this);
  command -> SetGuidance ("Specify coordinates x, y, z and world size");
  param   =  new G4UIparameter ("x", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("size", 'd', true);
  param   -> SetDefaultValue (10.);
  command -> SetParameter (param);
  fCommandList.push_back (command);

  ////////////////////////////////////////  /vis~/expert/draw_line  ////
  //expert \hline
  //expert /vis~/expert/draw\_line & $x_0$, $y_0$, $z_0$, $x_1$, $y_1$, $z_1$ &
  //expert Specify coordinates of start and end. \\%
  command = new G4UIcommand ("/vis~/expert/draw_line", this);
  command -> SetGuidance ("Specify coordinates x0, y0, z0,  x1, y1 and z1");
  param   =  new G4UIparameter ("x0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("x1", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y1", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z1", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  fCommandList.push_back (command);

  ////////////////////////////////////  /vis~/expert/draw_marks_and_show  ////
  //expert \hline
  //expert /vis~/expert/ draw\_marks\_and\_show & $x$, $y$, $z$, $size$ &
  //expert Specify coordinates and world size. \\%
  command = new G4UIcommand ("/vis~/expert/draw_marks_and_show", this);
  command -> SetGuidance ("Specify coordinates x, y, z and world size");
  param   =  new G4UIparameter ("x", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("size", 'd', true);
  param   -> SetDefaultValue (10.);
  command -> SetParameter (param);
  fCommandList.push_back (command);

  /***************************************
  ////////////////////////////////////  /vis~/expert/draw_physical_volume  ////
  //expert \hline
  //expert /vis~/expert/ draw\_physical\_volume &&
  //expert A test of
  //expert {\tt G4VisManager::Draw (const G4VPhysicalVolume\&},... \\%
  command = new G4UIcommand ("/vis~/expert/draw_physical_volume", this);
  command -> SetGuidance
    (
     "A test of G4VisManager::Draw (const G4VPhysicalVolume&,..."
     );
  fCommandList.push_back (command);
  **************************************/

  ////////////////////////////////////////  /vis~/expert/draw_polymarkers  ////
  //expert \hline
  //expert /vis~/expert/draw\_polymarkers & $N_\mathrm{spirals}$ &
  //expert Draws {\tt G4Polymarker}s in the form of coloured spirals. \\%
  command = new G4UIcommand ("/vis~/expert/draw_polymarkers", this);
  command -> SetGuidance ("Takes 2 parameters: no. of random sets, type");
  param   =  new G4UIparameter ("no. of polymarker sets", 'i', true);
  param   -> SetDefaultValue (1);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("marker type", 'i', true);
  param   -> SetGuidance ("0/1/2/3 = line/dots/circles/squares");
  param   -> SetDefaultValue (0);
  command -> SetParameter (param);
  fCommandList.push_back (command);

  ////////////////////////////////////////  /vis~/expert/draw_spiral  ////
  //expert \hline
  //expert /vis~/expert/draw\_spiral &&
  //expert Draws a {\tt G4Polyline} in the form of a spiral. \\%
  command = new G4UIcommand ("/vis~/expert/draw_spiral", this);
  command -> SetGuidance
    (
     "Draws a G4Polyline in the form of a spiral."
     );
  fCommandList.push_back (command);

  ////////////////////////////////////////  /vis~/expert/draw_spirals  ////
  //expert \hline
  //expert /vis~/expert/draw\_spirals & $N_\mathrm{spirals}$ &
  //expert Draws {\tt G4Polyline}s in the form of coloured spirals. \\%
  command = new G4UIcommand ("/vis~/expert/draw_spirals", this);
  command -> SetGuidance
    (
     "Draws G4Polylines in the form of coloured spirals."
     );
  param   =  new G4UIparameter ("no. of spirals", 'i', true);
  param   -> SetDefaultValue (1);
  command -> SetParameter (param);
  fCommandList.push_back (command);

  ////////////////////////////////////////  /vis~/expert/draw_square  ////
  //expert \hline
  //expert /vis~/expert/draw\_square & $x$, $y$, $z$, $size$ &
  //expert Specify coordinates and world size. \\%
  command = new G4UIcommand ("/vis~/expert/draw_square", this);
  command -> SetGuidance ("Specify coordinates x, y, z and world size");
  param   =  new G4UIparameter ("x", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("size", 'd', true);
  param   -> SetDefaultValue (10.);
  command -> SetParameter (param);
  fCommandList.push_back (command);

}

void G4VisManMessenger::DoCommandExpert (const G4String& commandPath,
					 G4String& newValues) {

  ///////////////////////////////////////  /vis~/expert/draw_circle  ////
  if (commandPath == "/vis~/expert/draw_circle") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      G4double x, y, z, worldSize;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString); is >> x >> y >> z >> worldSize;
      G4Circle circle = G4Circle(G4Point3D(x,y,z));
      circle.SetWorldSize (worldSize);
      G4Colour c(0.0, 1.0, 1.0 );
      G4VisAttributes a(c);
      circle.SetVisAttributes(&a);
      fpVMan->Draw (circle);
    }
  }

  ///////////////////////////////////////  /vis~/expert/draw_line  ////
  if (commandPath == "/vis~/expert/draw_line") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      G4double x0, y0, z0;
      G4double x1, y1, z1;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString); is >> x0 >> y0 >> z0 >> x1 >> y1 >> z1 ;

      G4Polyline line;
      G4Colour c(1.0, 0.0, 0.0);
      G4VisAttributes a(c);
      line.SetVisAttributes(&a);
      line.append (G4Point3D (x0,y0,z0));
      line.append (G4Point3D (x1,y1,z1));
      fpVMan -> Draw (line);
      line.clear ();
    }
  }

  ////////////////////////////////////  /vis~/expert/draw_marks_and_show  ////
  if (commandPath == "/vis~/expert/draw_marks_and_show") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");

    const G4int     num_each_marks = 100 ;
    const G4double  mark_step      = 1000.0 ;
    G4int             i                   ;
    G4double xx, yy, zz ;
	

    if (ViewValid ()) {
      G4double x, y, z, worldSize;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString); is >> x >> y >> z >> worldSize;

       // Draw squares
      xx = x ; yy = y; zz = z ;
      for( i = 0 ; i < num_each_marks ; i++ ) {
         G4Square Square = G4Square(G4Point3D(xx,yy,zz));
         Square.SetWorldSize (worldSize);
	 G4Colour c(1.0, 0.0, 0.0 );
	 G4VisAttributes a(c);
	 Square.SetVisAttributes(&a);
         fpVMan->Draw (Square);

	 xx += mark_step ;
         yy  = y ;    zz = z ;         
      } // for

       // Draw circles
      xx = x ; yy = y; zz = z ;
      for( i = 0 ; i < num_each_marks ; i++ ) {
         G4Circle circle = G4Circle (G4Point3D(xx,yy,zz));
         circle.SetWorldSize (worldSize);
	 G4Colour c(0.0, 1.0, 1.0 );
	 G4VisAttributes a(c);
	 circle.SetVisAttributes(&a);
         fpVMan->Draw (circle);

	 yy += mark_step ;
         xx  = x ;    zz = z ;         
      } // for

        // show
      fpVMan -> Show (); 

    } // if
  }

  /****************************************
  //////////////////////////////////  /vis~/expert/draw_physical_volume  ////
  if (commandPath == "/vis~/expert/draw_physical_volume") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {

      // Given a physical volume...
      G4VPhysicalVolume* pPhysVol =
	fpVMan -> GetCurrentScene ().GetPhysicalVolume ();

      // Given its transformation (normally, this would be known from
      // the hits or digis information)...
      HepRotation m;
      Hep3Vector axis (G4UniformRand(), G4UniformRand(), G4UniformRand());
      m.rotate (twopi * G4UniformRand(), axis);
      G4double scale = 1000. * G4UniformRand();
      Hep3Vector v (scale * G4UniformRand(),
		    scale * G4UniformRand(),
		    scale * G4UniformRand());
      G4Transform3D transform (m, v);

      // Change attributes - begin snippet.
      // pPhysVol is the pointer to the physical volume.
      // transform is the required G4Transform3D in the world
      //   coordinate system.

      // Find the corresponding logical volume...
      G4LogicalVolume* pLV  = pPhysVol -> GetLogicalVolume ();

      // Copy its visualization attributes into local variable...
      const G4VisAttributes* pVA = pLV -> GetVisAttributes ();
      G4VisAttributes attribs;     // Default vis attributes.
      if (pVA) attribs = *pVA;     // Set to those of logical volume if any.

      // Modify the visualization attributes...
      G4double red   = G4UniformRand();
      G4double green = G4UniformRand();
      G4double blue  = G4UniformRand();
      G4Colour c (red, green, blue);
      attribs.SetColour (c);
      attribs.SetForceSolid (true);

      // Draw it...
      fpVMan -> Draw (*pPhysVol, attribs, transform);
      // Change attributes - end snippet.
    }
  }
  ***********************************/

  ///////////////////////////////////////  /vis~/expert/draw_polymarkers  ////
  if (commandPath == "/vis~/expert/draw_polymarkers") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      int nPolyMarkers, type;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ;is >> nPolyMarkers >> type;
      G4Polymarker::MarkerType realType;
      switch (type) {
      default:
      case 0: realType = G4Polymarker::line; break;
      case 1: realType = G4Polymarker::dots; break;
      case 2: realType = G4Polymarker::circles; break;
      case 3: realType = G4Polymarker::squares; break;
      }
      G4int i;
      G4double a, fa, z, dz, t, dt = M_PI / 20.;
      for (int iPolyMarker = 0; iPolyMarker < nPolyMarkers; iPolyMarker++) {
	G4double x0 = 2000. * (G4UniformRand() - 0.5);
	G4double y0 = 2000. * (G4UniformRand() - 0.5);
	G4double z0 = 2000. * (G4UniformRand() - 0.5);
	a = 1000. * G4UniformRand();
	fa = 0.99 + 0.01 * G4UniformRand();
	z = 500. * (G4UniformRand() - 0.5);
	dz = 50. * (G4UniformRand() - 0.5);
	G4Polymarker polymarker;
	G4int nPoints = G4int (100 * G4UniformRand());
	G4double tStart = twopi * G4UniformRand();
	for (i = 0, t = tStart; i < nPoints; t += dt, a *= fa, z += dz, i++) {
	  polymarker.append (G4Point3D (x0 + a * cos (t),
					y0 + a * sin (t), z0 + z));
	}
	G4double red   = G4UniformRand();
	G4double green = G4UniformRand();
	G4double blue  = G4UniformRand();
	G4Colour c (red, green, blue);
	G4VisAttributes a(c);
	polymarker.SetVisAttributes(&a);
	polymarker.SetMarkerType (realType);
	fpVMan -> Draw (polymarker);
      }
    }
  }

  ///////////////////////////////////////  /vis~/expert/draw_spiral  ////
  if (commandPath == "/vis~/expert/draw_spiral") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      G4Polyline line;
      G4double t = 0., dt = M_PI / 50.;
      G4double a = 2000., fa = 0.999;
      G4double z = -500., dz = 1.;
      for (int i = 0; i < 1000; t += dt, a *= fa, z += dz, i++) {
        line.append (G4Point3D (a * cos (t), a * sin (t), z));
      }
      fpVMan -> Draw (line);
    }
  }

  ///////////////////////////////////////  /vis~/expert/draw_spirals  ////
  if (commandPath == "/vis~/expert/draw_spirals") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      int nSpirals;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> nSpirals;
      G4int i;
      G4double a, fa, z, dz, t, dt = M_PI / 50.;
      G4Polyline line;
      for (int iSpiral = 0; iSpiral < nSpirals; iSpiral++) {
	G4double x0 = 2000. * (G4UniformRand() - 0.5);
	G4double y0 = 2000. * (G4UniformRand() - 0.5);
	G4double z0 = 2000. * (G4UniformRand() - 0.5);
	a = 1000. * G4UniformRand();
	fa = 0.9985 + 0.001 * G4UniformRand();
	z = 500. * (G4UniformRand() - 0.5);
	dz = 5. * (G4UniformRand() - 0.5);
	G4int nPoints = G4int (1000 * G4UniformRand());
	G4double tStart = twopi * G4UniformRand();
	for (i = 0, t = tStart; i < nPoints; t += dt, a *= fa, z += dz, i++) {
	  line.append (G4Point3D (x0 + a * cos (t), y0 + a * sin (t), z0 + z));
	}
	G4Colour c (G4UniformRand(), G4UniformRand(), G4UniformRand());
	G4VisAttributes a(c);
	line.SetVisAttributes(&a);
	fpVMan -> Draw (line);
	line.clear ();
      }
    }
  }

  ///////////////////////////////////////  /vis~/expert/draw_square  ////
  if (commandPath == "/vis~/expert/draw_square") {
    G4VisManager::PrintCommandDeprecation
      ("This command will no longer be maintained.");
    if (ViewValid ()) {
      G4double x, y, z, worldSize;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString); is >> x >> y >> z >> worldSize;
      G4Square Square = G4Square (G4Point3D (x, y, z));
      // Set attributes - begin snippet.
      Square.SetWorldSize (worldSize);
      G4Colour c (1.0, 0.0, 0.0 );
      G4VisAttributes a(c);
      Square.SetVisAttributes(&a);
      // Set attributes - end snippet.
      fpVMan->Draw (Square);
    }
  }

}
