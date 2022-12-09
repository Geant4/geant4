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

// /vis/viewer/set commands - John Allison  16th May 2000

#include "G4VisCommandsViewerSet.hh"

#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UnitsTable.hh"
#include "G4VisManager.hh"
#include "G4Polyhedron.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>
#include <iomanip>

#define G4warn G4cout

G4VisCommandsViewerSet::G4VisCommandsViewerSet ():
fLightsVector    (G4ThreeVector(1.,1.,1.)),
fUpVector        (G4ThreeVector(0.,1.,0.)),
fViewpointVector (G4ThreeVector(0.,0.,1.))
{
  G4bool omitable;
  G4UIparameter* parameter;

  fpCommandAll = new G4UIcmdWithAString ("/vis/viewer/set/all",this);
  fpCommandAll->SetGuidance
  ("Copies view parameters.");
  fpCommandAll->SetGuidance
  ("Copies ALL view parameters (except the autoRefresh status) from"
   "\nfrom-viewer to current viewer. You may need \"/vis/viewer/rebuild\".");
  fpCommandAll->SetGuidance
  ("Note: to copy only the camera-specific parameters use"
   "\n\"/vis/viewer/copyfrom\".");
  fpCommandAll->SetParameterName ("from-viewer-name",omitable = false);

  fpCommandAutoRefresh = new G4UIcmdWithABool
  ("/vis/viewer/set/autoRefresh",this);
  fpCommandAutoRefresh->SetGuidance("Sets auto-refresh.");
  fpCommandAutoRefresh->SetGuidance
  ("If true, view is automatically refreshed after a change of"
   "\nview parameters.");
  fpCommandAutoRefresh->SetParameterName("auto-refresh",omitable = true);
  fpCommandAutoRefresh->SetDefaultValue(true);

  fpCommandAuxEdge = new G4UIcmdWithABool
  ("/vis/viewer/set/auxiliaryEdge",this);
  fpCommandAuxEdge->SetGuidance("Sets visibility of auxiliary edges");
  fpCommandAuxEdge->SetGuidance
  ("Auxiliary edges, i.e., those that are part of a curved surface,"
   "\nsometimes called soft edges, become visible/invisible.");
  fpCommandAuxEdge->SetParameterName("edge",omitable = true);
  fpCommandAuxEdge->SetDefaultValue(true);

  fpCommandBackground = new G4UIcommand
  ("/vis/viewer/set/background",this);
  fpCommandBackground->SetGuidance
  ("Set background colour and transparency (default black and opaque).");
  fpCommandBackground->SetGuidance(ConvertToColourGuidance());
  parameter = new G4UIparameter("red_or_string", 's', omitable = true);
  parameter -> SetDefaultValue ("0.");
  fpCommandBackground -> SetParameter (parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommandBackground -> SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommandBackground -> SetParameter (parameter);
  parameter = new G4UIparameter ("opacity", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandBackground -> SetParameter (parameter);

  fpCommandCulling = new G4UIcommand("/vis/viewer/set/culling",this);
  fpCommandCulling->SetGuidance ("Set culling options.");
  fpCommandCulling->SetGuidance
  ("\"global\": enables/disables all other culling options.");
  fpCommandCulling->SetGuidance
  ("\"coveredDaughters\": culls, i.e., eliminates, volumes that would not"
   "\nbe seen because covered by ancestor volumes in surface drawing mode,"
   "\nand then only if the ancestors are visible and opaque, and then only"
   "\nif no sections or cutaways are in operation.  Intended solely to"
   "\nimprove the speed of rendering visible volumes.");
  fpCommandCulling->SetGuidance
  ("\"invisible\": culls objects with the invisible attribute set.");
  fpCommandCulling->SetGuidance
  ("\"density\": culls volumes with density lower than threshold.  Useful"
   "\nfor eliminating \"container volumes\" with no physical correspondence,"
   "\nwhose material is usually air.  If this is selected, provide threshold"
   "\ndensity and unit (e.g., g/cm3, mg/cm3 or kg/m3)."
   );
  parameter = new G4UIparameter("culling-option",'s',omitable = false);
  parameter->SetParameterCandidates
  ("global coveredDaughters invisible density");
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("action",'b',omitable = true);
  parameter->SetDefaultValue(1);
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("density-threshold",'d',omitable = true);
  parameter->SetDefaultValue("0.01");
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("unit",'s',omitable = true);
  // parameter->SetParameterCandidates ("g/cm3 mg/cm3 kg/m3");
  // Instead of the above, SetNewValue accepts *any* density unit.
  parameter->SetDefaultValue("g/cm3");
  fpCommandCulling->SetParameter(parameter);

  fpCommandCutawayMode =
  new G4UIcmdWithAString ("/vis/viewer/set/cutawayMode", this);
  fpCommandCutawayMode->SetGuidance
  ("Sets cutaway mode - add (union) or multiply (intersection).");
  fpCommandCutawayMode->SetParameterName ("cutaway-mode",omitable = false);
  fpCommandCutawayMode->SetCandidates ("add union multiply intersection");
  fpCommandCutawayMode->SetDefaultValue("union");

  fpCommandDefaultColour = new G4UIcommand
  ("/vis/viewer/set/defaultColour",this);
  fpCommandDefaultColour->SetGuidance
  ("Set defaultColour colour and transparency (default white and opaque).");
  fpCommandDefaultColour->SetGuidance(ConvertToColourGuidance());
  parameter = new G4UIparameter("red_or_string", 's', omitable = true);
  parameter -> SetDefaultValue ("1.");
  fpCommandDefaultColour -> SetParameter (parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandDefaultColour -> SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandDefaultColour -> SetParameter (parameter);
  parameter = new G4UIparameter ("opacity", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandDefaultColour -> SetParameter (parameter);

  fpCommandDefaultTextColour = new G4UIcommand
  ("/vis/viewer/set/defaultTextColour",this);
  fpCommandDefaultTextColour->SetGuidance
  ("Set defaultTextColour colour and transparency (default blue and opaque).");
  fpCommandDefaultTextColour->SetGuidance(ConvertToColourGuidance());
  parameter = new G4UIparameter("red_or_string", 's', omitable = true);
  parameter -> SetDefaultValue ("0.");
  fpCommandDefaultTextColour -> SetParameter (parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommandDefaultTextColour -> SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandDefaultTextColour -> SetParameter (parameter);
  parameter = new G4UIparameter ("opacity", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandDefaultTextColour -> SetParameter (parameter);

  fpCommandEdge = new G4UIcmdWithABool("/vis/viewer/set/edge",this);
  fpCommandEdge->SetGuidance
  ("Edges become visible/invisible in surface mode.");
  fpCommandEdge->SetParameterName("edge",omitable = true);
  fpCommandEdge->SetDefaultValue(true);

  fpCommandExplodeFactor = new G4UIcommand
  ("/vis/viewer/set/explodeFactor", this);
  fpCommandExplodeFactor->SetGuidance
  ("Moves top-level drawn volumes by this factor from this centre.");
  parameter = new G4UIparameter("explodeFactor", 'd', omitable=true);
  parameter->SetParameterRange("explodeFactor>=1.");
  parameter->SetDefaultValue(1.);
  fpCommandExplodeFactor->SetParameter(parameter);
  parameter = new G4UIparameter("x",'d',omitable = true);
  parameter->SetDefaultValue  (0);
  parameter->SetGuidance      ("Coordinate of explode centre.");
  fpCommandExplodeFactor->SetParameter(parameter);
  parameter = new G4UIparameter("y",'d',omitable = true);
  parameter->SetDefaultValue  (0);
  parameter->SetGuidance      ("Coordinate of explode centre.");
  fpCommandExplodeFactor->SetParameter(parameter);
  parameter = new G4UIparameter("z",'d',omitable = true);
  parameter->SetDefaultValue  (0);
  parameter->SetGuidance      ("Coordinate of explode centre.");
  fpCommandExplodeFactor->SetParameter(parameter);
  parameter = new G4UIparameter("unit",'s',omitable = true);
  parameter->SetDefaultValue  ("m");
  parameter->SetGuidance      ("Unit of explode centre.");
  fpCommandExplodeFactor->SetParameter(parameter);

  fpCommandGlobalLineWidthScale = new G4UIcmdWithADouble
  ("/vis/viewer/set/globalLineWidthScale", this);
  fpCommandGlobalLineWidthScale->SetGuidance
  ("Multiplies line widths by this factor.");
  fpCommandGlobalLineWidthScale->
  SetParameterName("scale-factor", omitable=true);
  fpCommandGlobalLineWidthScale->SetDefaultValue(1.);

  fpCommandGlobalMarkerScale = new G4UIcmdWithADouble
  ("/vis/viewer/set/globalMarkerScale", this);
  fpCommandGlobalMarkerScale->SetGuidance
  ("Multiplies marker sizes by this factor.");
  fpCommandGlobalMarkerScale->
  SetParameterName("scale-factor", omitable=true);
  fpCommandGlobalMarkerScale->SetDefaultValue(1.);

  fpCommandHiddenEdge =
  new G4UIcmdWithABool("/vis/viewer/set/hiddenEdge",this);
  fpCommandHiddenEdge->SetGuidance
  ("Edges become hidden/seen in wireframe or surface mode.");
  fpCommandHiddenEdge->SetParameterName("hidden-edge",omitable = true);
  fpCommandHiddenEdge->SetDefaultValue(true);

  fpCommandHiddenMarker =
  new G4UIcmdWithABool("/vis/viewer/set/hiddenMarker",this);
  fpCommandHiddenMarker->SetGuidance
  ("If true, closer objects hide markers. Otherwise, markers always show.");
  fpCommandHiddenMarker->SetParameterName("hidden-marker",omitable = true);
  fpCommandHiddenMarker->SetDefaultValue(true);

  fpCommandLightsMove = new G4UIcmdWithAString
  ("/vis/viewer/set/lightsMove",this);
  fpCommandLightsMove->SetGuidance
  ("Lights move with camera or with object");
  fpCommandLightsMove->SetParameterName("lightsMove",omitable = false);
  fpCommandLightsMove->SetCandidates
  ("cam camera with-camera obj object with-object");

  fpCommandLightsThetaPhi = new G4UIcommand
  ("/vis/viewer/set/lightsThetaPhi", this);
  fpCommandLightsThetaPhi->SetGuidance
  ("Set direction from target to lights.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetDefaultValue(60.);
  fpCommandLightsThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetDefaultValue(45.);
  fpCommandLightsThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandLightsThetaPhi -> SetParameter (parameter);

  fpCommandLightsVector = new G4UIcommand
  ("/vis/viewer/set/lightsVector", this);
  fpCommandLightsVector->SetGuidance
  ("Set direction from target to lights.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetDefaultValue (1);
  fpCommandLightsVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetDefaultValue (1);
  fpCommandLightsVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetDefaultValue (1);
  fpCommandLightsVector -> SetParameter (parameter);

  fpCommandLineSegments = new G4UIcmdWithAnInteger
  ("/vis/viewer/set/lineSegmentsPerCircle",this);
  fpCommandLineSegments->SetGuidance
  ("Set number of sides per circle for polygon/polyhedron drawing.");
  fpCommandLineSegments->SetGuidance
  ("Refers to graphical representation of objects with curved lines/surfaces.");
  fpCommandLineSegments->SetParameterName("line-segments",omitable = true);
  fpCommandLineSegments->SetDefaultValue(24);

  fpCommandLineWidth = new G4UIcmdWithoutParameter
  ("/vis/viewer/set/lineWidth",this);
  fpCommandLineWidth->SetGuidance
  ("Use \"/vis/viewer/set/globalLineWidthScale\" instead."
   "\nFor trajectories use \"/vis/modeling/trajectories/*/default/setLineWidth\"."
   "\nFor volumes use \"/vis/geometry/set/lineWidth\".");

  fpCommandNumberOfCloudPoints = new G4UIcmdWithAnInteger
  ("/vis/viewer/set/numberOfCloudPoints",this);
  fpCommandNumberOfCloudPoints->SetGuidance
  ("Set number of points to be used for cloud representation of volumes.");
  fpCommandNumberOfCloudPoints->SetParameterName("points",omitable = true);
  fpCommandNumberOfCloudPoints->SetDefaultValue(10000);

  fpCommandPicking = new G4UIcmdWithABool
  ("/vis/viewer/set/picking",this);
  fpCommandPicking->SetGuidance("Sets picking, if available.");
  fpCommandPicking->SetGuidance
  ("If true, view is set up for picking, if available.");
  fpCommandPicking->SetGuidance
  ("You may need to issue \"/vis/viewer/update\".");
  fpCommandPicking->SetGuidance
  ("For required actions, watch for instructions for viewer.");
  fpCommandPicking->SetParameterName("picking",omitable = true);
  fpCommandPicking->SetDefaultValue(true);

  fpCommandProjection = new G4UIcommand("/vis/viewer/set/projection",this);
  fpCommandProjection->SetGuidance
  ("Set projection style - o[rthogonal] or p[erspective]."
   "\nIf p[erspective], also set field half angle.");
  parameter = new G4UIparameter("projection",'s',omitable = true);
  parameter->SetParameterCandidates("o orthogonal p perspective");
  parameter->SetDefaultValue("orthogonal");
  fpCommandProjection->SetParameter(parameter);
  parameter = new G4UIparameter("field-half-angle",'d',omitable = true);
  parameter->SetDefaultValue(30.);
  //parameter->SetCurrentAsDefault(true);
  fpCommandProjection->SetParameter(parameter);
  parameter = new G4UIparameter("unit",'s',omitable = true);
  parameter->SetDefaultValue("deg");
  //parameter->SetCurrentAsDefault(true);
  fpCommandProjection->SetParameter(parameter);

  fpCommandRotationStyle = new G4UIcmdWithAString
  ("/vis/viewer/set/rotationStyle",this);
  fpCommandRotationStyle->SetGuidance
  ("Set style of rotation - constrainUpDirection or freeRotation.");
  fpCommandRotationStyle->SetGuidance
  ("constrainUpDirection: conventional HEP view.");
  fpCommandRotationStyle->SetGuidance
  ("freeRotation: Google-like rotation, using mouse-grab.");
  fpCommandRotationStyle->SetParameterName ("style",omitable = false);
  fpCommandRotationStyle->SetCandidates("constrainUpDirection freeRotation");

  fpCommandSectionPlane = new G4UIcommand("/vis/viewer/set/sectionPlane",this);
  fpCommandSectionPlane -> SetGuidance
  ("Set plane for drawing section (DCUT).");
  fpCommandSectionPlane -> SetGuidance
  ("E.g., for a y-z plane at x = 1 cm:"
   "\n\"/vis/viewer/set/sectionPlane on 1 0 0 cm 1 0 0\"."
   "\nTo turn off: /vis/viewer/set/sectionPlane off");
  parameter  =  new G4UIparameter("Selector",'c',true);
  parameter  -> SetDefaultValue  ("on");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("x",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("y",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("z",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("unit",'s',omitable = true);
  parameter  -> SetDefaultValue  ("m");
  parameter  -> SetGuidance      ("Unit of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("nx",'d',omitable = true);
  parameter  -> SetDefaultValue  (1);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("ny",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("nz",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);

  fpCommandSpecialMeshRendering = new G4UIcmdWithABool
  ("/vis/viewer/set/specialMeshRendering",this);
  fpCommandSpecialMeshRendering -> SetGuidance
  ("Request special rendering of volumes (meshes) that use G4VParameterisation.");
  fpCommandSpecialMeshRendering->SetParameterName("render",omitable = true);
  fpCommandSpecialMeshRendering->SetDefaultValue(true);

  fpCommandSpecialMeshRenderingOption = new G4UIcmdWithAString
  ("/vis/viewer/set/specialMeshRenderingOption",this);
  fpCommandSpecialMeshRenderingOption->SetGuidance
  ("Set special mesh rendering option - \"dots\" or \"surfaces\".");
  fpCommandSpecialMeshRenderingOption->SetParameterName ("option",omitable = true);
  fpCommandSpecialMeshRenderingOption->SetCandidates("dots surfaces");
  fpCommandSpecialMeshRenderingOption->SetDefaultValue("dots");

  fpCommandSpecialMeshVolumes = new G4UIcommand
  ("/vis/viewer/set/specialMeshVolumes",this);
  fpCommandSpecialMeshVolumes -> SetGuidance
  ("Specify the volumes for special rendering. No arguments resets the list"
   "\nand is interpreted to mean \"all found meshes\".");
  fpCommandSpecialMeshVolumes->SetGuidance
  ("Please provide a list of space-separated physical volume names and copy"
   "\nnumber pairs. Negative copy number means \"all volumes of that name\".");
  parameter = new G4UIparameter("volumes",'s',omitable = true);
  parameter->SetGuidance
  ("List of physical volume names and copy number pairs");
  fpCommandSpecialMeshVolumes->SetParameter(parameter);

  fpCommandStyle = new G4UIcmdWithAString ("/vis/viewer/set/style",this);
  fpCommandStyle->SetGuidance
  ("Set style of drawing - w[ireframe] or s[urface] or c[loud].");
  fpCommandStyle->SetGuidance
  ("(Hidden line drawing is controlled by \"/vis/viewer/set/hiddenEdge\".)");
  fpCommandStyle->SetParameterName ("style",omitable = false);
  fpCommandStyle->SetCandidates("w wireframe s surface c cloud");

  fpCommandTargetPoint = new G4UIcmdWith3VectorAndUnit
  ("/vis/viewer/set/targetPoint", this);
  fpCommandTargetPoint->SetGuidance
  ("Set target point.");
  fpCommandTargetPoint->SetGuidance
  ("This sets the \"Current Target Point\" relative to the \"Standard");
  fpCommandTargetPoint->SetGuidance
  ("Target Point\" so that the actual target point is as requested.");
  fpCommandTargetPoint->SetGuidance
  ("(See G4ViewParameters.hh for an explanation of target points.)");
  fpCommandTargetPoint->SetParameterName("x", "y", "z", omitable = false);
  fpCommandTargetPoint->SetUnitCategory("Length");

  fpCommandUpThetaPhi = new G4UIcommand
  ("/vis/viewer/set/upThetaPhi", this);
  fpCommandUpThetaPhi -> SetGuidance ("Set up vector.");
  fpCommandUpThetaPhi -> SetGuidance
  ("Viewer will attempt always to show this direction upwards.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetDefaultValue (90.);
  fpCommandUpThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetDefaultValue (90.);
  fpCommandUpThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandUpThetaPhi -> SetParameter (parameter);

  fpCommandUpVector = new G4UIcommand
  ("/vis/viewer/set/upVector", this);
  fpCommandUpVector -> SetGuidance ("Set up vector.");
  fpCommandUpVector -> SetGuidance
  ("Viewer will attempt always to show this direction upwards.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommandUpVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandUpVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommandUpVector -> SetParameter (parameter);

  fpCommandViewpointThetaPhi = new G4UIcommand
  ("/vis/viewer/set/viewpointThetaPhi", this);
  fpCommandViewpointThetaPhi -> SetGuidance
  ("Set direction from target to camera.");
  fpCommandViewpointThetaPhi -> SetGuidance
  ("Also changes lightpoint direction if lights are set to move with camera.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetDefaultValue (60.);
  fpCommandViewpointThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetDefaultValue (45.);
  fpCommandViewpointThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandViewpointThetaPhi -> SetParameter (parameter);

  fpCommandViewpointVector = new G4UIcommand
  ("/vis/viewer/set/viewpointVector", this);
  fpCommandViewpointVector -> SetGuidance
  ("Set direction from target to camera.");
  fpCommandViewpointVector -> SetGuidance
  ("Also changes lightpoint direction if lights are set to move with camera.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandViewpointVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandViewpointVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetDefaultValue (1.);
  fpCommandViewpointVector -> SetParameter (parameter);

  fpTimeWindowDirectory = new G4UIdirectory ("/vis/viewer/set/timeWindow/");
  fpTimeWindowDirectory -> SetGuidance ("Set time window parameters of current viewer.");
  G4String timeWindowGuidance =
  "For these commands use"
  "\n  /vis/scene/add/trajectories rich"
  "\n  /vis/modeling/trajectories/drawByCharge-0/default/setTimeSliceInterval 0.01 ns"
  "\nthen typically"
  "\n  /vis/viewer/set/timeWindow/displayLightFront true 0 0 -50 cm -0.5 ns"
  "\n  /vis/viewer/set/timeWindow/displayHeadTime true"
  "\n  /vis/viewer/set/timeWindow/fadeFactor 1"
  "\n  /run/beamOn # or several until you get a good event or events"
  "\n  /vis/viewer/set/timeWindow/startTime 0 ns 1 ns"
  "\n  /vis/viewer/save"
  "\n  /vis/viewer/set/timeWindow/startTime 1 ns 1 ns"
  "\nthen zoom, pan etc to a view of interest and"
  "\n  /vis/viewer/save"
  "\nthen repeat with next start time, another view and a save, then try"
  "\n  /vis/viewer/interpolate";

  fpCommandTimeWindowDisplayHeadTime =
  new G4UIcommand("/vis/viewer/set/timeWindow/displayHeadTime", this);
  fpCommandTimeWindowDisplayHeadTime->SetGuidance
  ("Display head time of range in 2D text.");
  fpCommandTimeWindowDisplayHeadTime->SetGuidance(timeWindowGuidance);
  parameter = new G4UIparameter ("displayHeadTime", 'b', omitable = false);
  parameter->SetDefaultValue(false);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenX", 'd', omitable = true);
  parameter->SetGuidance("-1 < screenX < 1");
  parameter->SetParameterRange("screenX >= -1. && screenX <= 1.");
  parameter->SetDefaultValue(-0.9);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenY", 'd', omitable = true);
  parameter->SetGuidance("-1 < screenY < 1");
  parameter->SetParameterRange("screenY >= -1. && screenY <= 1.");
  parameter->SetDefaultValue(-0.9);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenSize", 'd', omitable = true);
  parameter->SetDefaultValue(24.);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetParameterRange("red >= 0. && red <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetParameterRange("green >= 0. && green <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetParameterRange("blue >= 0. && blue <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandTimeWindowDisplayHeadTime->SetParameter(parameter);

  fpCommandTimeWindowDisplayLightFront =
  new G4UIcommand("/vis/viewer/set/timeWindow/displayLightFront", this);
  fpCommandTimeWindowDisplayLightFront->SetGuidance
  ("Display the light front at head time.");
  fpCommandTimeWindowDisplayLightFront->SetGuidance
  ("Tip: The trajectories can appear of jump ahead of the light front"
   "\nbecause their time range overlaps the viewer's time range.  To"
   "\naverage out this discrete time effect, advance the light front by"
   "\nhalf the trajectories interval. E.g., if the trajectory time slice"
   "\ninterval is 0.01 ns:"
   "\n  /vis/viewer/set/timeWindow/displayLightFront true -90 0 0 mm -0.005 ns"
   "\nTo prevent them beating the light front at all:"
   "\n  /vis/viewer/set/timeWindow/displayLightFront true -90 0 0 mm -0.01 ns");
  fpCommandTimeWindowDisplayLightFront->SetGuidance(timeWindowGuidance);
  parameter = new G4UIparameter ("displayLightFront", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originX", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originY", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originZ", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("space_unit", 's', omitable = true);
  parameter->SetDefaultValue("m");
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originT", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("time_unit", 's', omitable = true);
  parameter->SetDefaultValue("s");
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetParameterRange("red >= 0. && red <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetParameterRange("green >= 0. && green <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetParameterRange("blue >= 0. && blue <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandTimeWindowDisplayLightFront->SetParameter(parameter);

  fpCommandTimeWindowEndTime =
  new G4UIcommand("/vis/viewer/set/timeWindow/endTime", this);
  fpCommandTimeWindowEndTime->SetGuidance("Set end and range of track time.");
  fpCommandTimeWindowEndTime->SetGuidance(timeWindowGuidance);
  parameter = new G4UIparameter ("end-time", 'd', omitable = false);
  parameter->SetDefaultValue(G4VisAttributes::fVeryLongTime);
  fpCommandTimeWindowEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("end-time-unit", 's', omitable = false);
  parameter->SetDefaultValue("ns");
  fpCommandTimeWindowEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range", 'd', omitable = true);
  parameter->SetDefaultValue(-1.);
  fpCommandTimeWindowEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range-unit", 's', omitable = true);
  parameter->SetDefaultValue("ns");
  fpCommandTimeWindowEndTime->SetParameter(parameter);

  fpCommandTimeWindowFadeFactor =
  new G4UIcmdWithADouble("/vis/viewer/set/timeWindow/fadeFactor", this);
  fpCommandTimeWindowFadeFactor->SetGuidance
  ("0: no fade; 1: maximum fade with time window.");
  fpCommandTimeWindowFadeFactor->SetGuidance(timeWindowGuidance);
  fpCommandTimeWindowFadeFactor->SetParameterName("fade_factor", omitable = false);
  fpCommandTimeWindowFadeFactor->SetRange("fade_factor>=0.&&fade_factor<=1.");
  fpCommandTimeWindowFadeFactor->SetDefaultValue(0.);

  fpCommandTimeWindowStartTime =
  new G4UIcommand("/vis/viewer/set/timeWindow/startTime", this);
  fpCommandTimeWindowStartTime->SetGuidance("Set start and range of track time.");
  fpCommandTimeWindowStartTime->SetGuidance(timeWindowGuidance);
  parameter = new G4UIparameter ("start-time", 'd', omitable = false);
  parameter->SetDefaultValue(-G4VisAttributes::fVeryLongTime);
  fpCommandTimeWindowStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("start-time-unit", 's', omitable = false);
  parameter->SetDefaultValue("ns");
  fpCommandTimeWindowStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range", 'd', omitable = true);
  parameter->SetDefaultValue(-1.);
  fpCommandTimeWindowStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range-unit", 's', omitable = true);
  parameter->SetDefaultValue("ns");
  fpCommandTimeWindowStartTime->SetParameter(parameter);
}

G4VisCommandsViewerSet::~G4VisCommandsViewerSet() {
  delete fpCommandTimeWindowStartTime;
  delete fpCommandTimeWindowFadeFactor;
  delete fpCommandTimeWindowEndTime;
  delete fpCommandTimeWindowDisplayLightFront;
  delete fpCommandTimeWindowDisplayHeadTime;
  delete fpTimeWindowDirectory;
  delete fpCommandViewpointVector;
  delete fpCommandViewpointThetaPhi;
  delete fpCommandUpVector;
  delete fpCommandUpThetaPhi;
  delete fpCommandTargetPoint;
  delete fpCommandStyle;
  delete fpCommandSpecialMeshVolumes;
  delete fpCommandSpecialMeshRenderingOption;
  delete fpCommandSpecialMeshRendering;
  delete fpCommandSectionPlane;
  delete fpCommandRotationStyle;
  delete fpCommandProjection;
  delete fpCommandPicking;
  delete fpCommandNumberOfCloudPoints;
  delete fpCommandLineWidth;
  delete fpCommandLineSegments;
  delete fpCommandLightsVector;
  delete fpCommandLightsThetaPhi;
  delete fpCommandLightsMove;
  delete fpCommandHiddenMarker;
  delete fpCommandHiddenEdge;
  delete fpCommandGlobalMarkerScale;
  delete fpCommandGlobalLineWidthScale;
  delete fpCommandExplodeFactor;
  delete fpCommandEdge;
  delete fpCommandDefaultTextColour;
  delete fpCommandDefaultColour;
  delete fpCommandCutawayMode;
  delete fpCommandCulling;
  delete fpCommandBackground;
  delete fpCommandAuxEdge;
  delete fpCommandAutoRefresh;
  delete fpCommandAll;
}

G4String G4VisCommandsViewerSet::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandsViewerSet::SetNewValue
(G4UIcommand* command,G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: G4VisCommandsViewerSet::SetNewValue: no current viewer."
      << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandAll) {
    G4VViewer* fromViewer = fpVisManager->GetViewer(newValue);
    if (!fromViewer) {
      if (verbosity >= G4VisManager::errors) {
        G4warn <<
        "ERROR: G4VisCommandsViewerSet::SetNewValue: all:"
        "\n  unrecognised from-viewer."
        << G4endl;
      }
      return;
    }
    if (fromViewer == currentViewer) {
      if (verbosity >= G4VisManager::warnings) {
        G4warn <<
        "WARNING: G4VisCommandsViewerSet::SetNewValue: all:"
        "\n  from-viewer and current viewer are identical."
        << G4endl;
      }
      return;
    }
    // Copy view parameters except for autoRefresh and background...
    auto keepAutoRefresh = vp.IsAutoRefresh();
    auto keepBackground  = vp.GetBackgroundColour();
    vp = fromViewer->GetViewParameters();
    vp.SetAutoRefresh(keepAutoRefresh);
    vp.SetBackgroundColour(keepBackground);
    // Concatenate any private vis attributes modifiers...
    const std::vector<G4ModelingParameters::VisAttributesModifier>*
    privateVAMs = fromViewer->GetPrivateVisAttributesModifiers();
    if (privateVAMs) {
      std::vector<G4ModelingParameters::VisAttributesModifier>::const_iterator i;
      for (i = privateVAMs->begin(); i != privateVAMs->end(); ++i) {
        vp.AddVisAttributesModifier(*i);
      }
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "View parameters of viewer \"" << currentViewer->GetName()
      << "\"\n  set to those of viewer \"" << fromViewer->GetName()
      << "\"."
      << G4endl;
    }
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "You may need \"/vis/viewer/rebuild\"."
      << G4endl;
    }
  }

  else if (command == fpCommandAutoRefresh) {
    G4bool autoRefresh = G4UIcommand::ConvertToBool(newValue);
    const G4ViewParameters& defaultVP =
    currentViewer->GetDefaultViewParameters();
    if (autoRefresh && !defaultVP.IsAutoRefresh()) {
      if (verbosity >= G4VisManager::warnings) {
        G4warn
        << "WARNING: "
        << currentViewer->GetName() << " is NOT auto-refesh by default"
        << "\n  so cannot be set to auto-refresh."
        << G4endl;
      }
      return;
    }
    vp.SetAutoRefresh(autoRefresh);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Views will ";
      if (!vp.IsAutoRefresh()) G4cout << "not ";
      G4cout << "be automatically refreshed after a change of view parameters."
      << G4endl;
    }
    if (!vp.IsAutoRefresh()) {
      currentViewer->SetViewParameters(vp);
      return;  // Avoid a refresh if auto-refresh has been set to off...
    }  // ...otherwise take normal action.
  }

  else if (command == fpCommandAuxEdge) {
    vp.SetAuxEdgeVisible(G4UIcommand::ConvertToBool(newValue));
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Auxiliary edges will ";
      if (!vp.IsAuxEdgeVisible()) G4cout << "not ";
      G4cout << "be visible." << G4endl;
    }
  }

  else if (command == fpCommandBackground) {
    G4String redOrString;
    G4double green, blue, opacity;
    std::istringstream iss(newValue);
    iss >> redOrString >> green >> blue >> opacity;
    G4Colour colour(0.,0.,0.);  // Default black and opaque.
    ConvertToColour(colour, redOrString, green, blue, opacity);
    vp.SetBackgroundColour(colour);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Background colour "
      << vp.GetBackgroundColour()
      << " requested."
      << G4endl;
    }
  }

  else if (command == fpCommandCulling) {
    G4String cullingOption, stringFlag, unit;
    G4double density;
    std::istringstream is (newValue);
    is >> cullingOption >> stringFlag >> density >> unit;
    G4bool boolFlag = G4UIcommand::ConvertToBool(stringFlag);
    if (cullingOption == "global") {
      vp.SetCulling(boolFlag);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout <<
        "G4VisCommandsViewerSet::SetNewValue: culling: global culling flag"
        " set to " << G4UIcommand::ConvertToString(boolFlag) <<
        ".\n  Does not change specific culling flags."
        << G4endl;
      }
    }
    else if (cullingOption == "coveredDaughters") {
      vp.SetCullingCovered(boolFlag);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout <<
        "G4VisCommandsViewerSet::SetNewValue: culling: culling covered"
        "\n  daughters flag set to "
        << G4UIcommand::ConvertToString(boolFlag) <<
        ".  Daughters covered by opaque mothers"
        "\n  will be culled, i.e., not drawn, if this flag is true."
        "\n  Note: this is only effective in surface drawing style,"
        "\n  and then only if the volumes are visible and opaque, and then"
        "\n  only if no sections or cutaways are in operation."
        << G4endl;
      }
    }
    else if (cullingOption == "invisible") {
      vp.SetCullingInvisible(boolFlag);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout <<
        "G4VisCommandsViewerSet::SetNewValue: culling: culling invisible"
        "\n  flag set to "
        << boolFlag << G4UIcommand::ConvertToString(boolFlag) <<
        ".  Volumes marked invisible will be culled,"
        "\n  i.e., not drawn, if this flag is true."
        << G4endl;
      }
    }
    else if (cullingOption == "density") {
      const G4String where =
      "G4VisCommandsViewerSet::SetNewValue: culling: culling by density";
      if (boolFlag) {
        G4double valueOfUnit;
        // "Volumic Mass" is Michel's phrase for "Density"
        if (ProvideValueOfUnit(where,unit,"Volumic Mass",valueOfUnit)) {
          // Successful outcome of unit search
          vp.SetDensityCulling(boolFlag);
          density *= valueOfUnit;
          vp.SetVisibleDensity(density);
        } else {
          // Unsuccessful outcome of unit search. Flag and density unchanged.
          density = vp.GetVisibleDensity();
        }
      } else {  // Reset flag but density unchanged.
        vp.SetDensityCulling(boolFlag);
      }
      if (verbosity >= G4VisManager::confirmations) {
        G4cout
        << where
        << "\n  flag set to "
        << std::boolalpha << vp.IsDensityCulling()
        << ".  Volumes with density less than "
        << G4BestUnit(density,"Volumic Mass")
        << "\n  will be culled, i.e., not drawn, if this flag is true."
        << G4endl;
      }
    }
    else {
      if (verbosity >= G4VisManager::errors) {
        G4warn <<
        "ERROR: G4VisCommandsViewerSet::SetNewValue: culling:"
        "\n  option not recognised."
        << G4endl;
      }
    }
  }

  else if (command == fpCommandCutawayMode) {
    if (newValue == "add" || newValue == "union")
      vp.SetCutawayMode(G4ViewParameters::cutawayUnion);
    if (newValue == "multiply" || newValue == "intersection")
      vp.SetCutawayMode(G4ViewParameters::cutawayIntersection);

    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Cutaway mode set to ";
      if (vp.GetCutawayMode() == G4ViewParameters::cutawayUnion)
        G4cout << "cutawayUnion";
      if (vp.GetCutawayMode() == G4ViewParameters::cutawayIntersection)
        G4cout << "cutawayIntersection";
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandDefaultColour) {
    G4String redOrString;
    G4double green, blue, opacity;
    std::istringstream iss(newValue);
    iss >> redOrString >> green >> blue >> opacity;
    G4Colour colour(1.,1.,1.);  // Default white and opaque.
    ConvertToColour(colour, redOrString, green, blue, opacity);
    G4VisAttributes va = vp.GetDefaultVisAttributes();
    va.SetColour(colour);
    vp.SetDefaultVisAttributes(va);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Default colour "
      << vp.GetDefaultVisAttributes()->GetColour()
      << " requested."
      << G4endl;
    }
  }

  else if (command == fpCommandDefaultTextColour) {
    G4String redOrString;
    G4double green, blue, opacity;
    std::istringstream iss(newValue);
    iss >> redOrString >> green >> blue >> opacity;
    G4Colour colour(1.,1.,1.);  // Default white and opaque.
    ConvertToColour(colour, redOrString, green, blue, opacity);
    G4VisAttributes va = vp.GetDefaultTextVisAttributes();
    va.SetColour(colour);
    vp.SetDefaultTextVisAttributes(va);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Default colour "
      << vp.GetDefaultTextVisAttributes()->GetColour()
      << " requested."
      << G4endl;
    }
  }

  else if (command == fpCommandEdge) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    if (G4UIcommand::ConvertToBool(newValue)) {  // true
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          break;
        case G4ViewParameters::hlr:
          break;
        case G4ViewParameters::hsr:
          vp.SetDrawingStyle(G4ViewParameters::hlhsr);
          break;
        case G4ViewParameters::hlhsr:
          break;
        case G4ViewParameters::cloud:
          break;
      }
    }
    else {  // false
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          break;
        case G4ViewParameters::hlr:
          break;
        case G4ViewParameters::hsr:
          break;
        case G4ViewParameters::hlhsr:
          vp.SetDrawingStyle(G4ViewParameters::hsr);
          break;
        case G4ViewParameters::cloud:
          break;
      }
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
      << "\" set to " << vp.GetDrawingStyle()
      << G4endl;
    }
  }

  else if (command == fpCommandExplodeFactor) {
    G4double explodeFactor, x, y, z;
    G4String unitString;
    std::istringstream is (newValue);
    is >> explodeFactor >> x >> y >> z >> unitString;
    G4double unit = G4UIcommand::ValueOf(unitString);
    vp.SetExplodeFactor(explodeFactor);
    vp.SetExplodeCentre(G4Point3D(x * unit, y * unit, z * unit));
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Explode factor changed to " << vp.GetExplodeFactor()
      << " from centre " << vp.GetExplodeCentre()
      << G4endl;
    }
  }

  else if (command == fpCommandGlobalLineWidthScale) {
    G4double globalLineWidthScale
    = fpCommandGlobalLineWidthScale->GetNewDoubleValue(newValue);
    vp.SetGlobalLineWidthScale(globalLineWidthScale);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Global Line Width Scale changed to "
      << vp.GetGlobalLineWidthScale() << G4endl;
    }
  }

  else if (command == fpCommandGlobalMarkerScale) {
    G4double globalMarkerScale
    = fpCommandGlobalMarkerScale->GetNewDoubleValue(newValue);
    vp.SetGlobalMarkerScale(globalMarkerScale);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Global Marker Scale changed to "
      << vp.GetGlobalMarkerScale() << G4endl;
    }
  }

  else if (command == fpCommandHiddenEdge) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    if (G4UIcommand::ConvertToBool(newValue)) {
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          vp.SetDrawingStyle(G4ViewParameters::hlr);
          break;
        case G4ViewParameters::hlr:
          break;
        case G4ViewParameters::hsr:
          vp.SetDrawingStyle(G4ViewParameters::hlhsr);
          break;
        case G4ViewParameters::hlhsr:
          break;
        case G4ViewParameters::cloud:
          break;
      }
    }
    else {
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          break;
        case G4ViewParameters::hlr:
          vp.SetDrawingStyle(G4ViewParameters::wireframe);
          break;
        case G4ViewParameters::hsr:
          break;
        case G4ViewParameters::hlhsr:
          vp.SetDrawingStyle(G4ViewParameters::hsr);
          break;
        case G4ViewParameters::cloud:
          break;
      }
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
      << "\" set to " << vp.GetDrawingStyle()
      << G4endl;
    }
  }

  else if (command == fpCommandHiddenMarker) {
    G4bool hidden = G4UIcommand::ConvertToBool(newValue);
    if (hidden) vp.SetMarkerHidden();
    else vp.SetMarkerNotHidden();
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Markers will ";
      if (vp.IsMarkerNotHidden()) G4cout << "not ";
      G4cout << "be hidden under solid objects." << G4endl;
    }
  }

  else if (command == fpCommandLightsMove) {
    if (newValue.find("cam") != G4String::npos)
      vp.SetLightsMoveWithCamera(true);
    else if(newValue.find("obj") != G4String::npos)
      vp.SetLightsMoveWithCamera(false);
    else {
      if (verbosity >= G4VisManager::errors) {
        G4warn << "ERROR: \"" << newValue << "\" not recognised."
        "  Looking for \"cam\" or \"obj\" in string." << G4endl;
      }
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Lights move with ";
      if (vp.GetLightsMoveWithCamera())
        G4cout << "camera (object appears to rotate).";
      else G4cout << "object (the viewer appears to be moving).";
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandLightsThetaPhi) {
    G4double theta, phi;
    if (ConvertToDoublePair(newValue, theta, phi)) {
      G4double x = std::sin (theta) * std::cos (phi);
      G4double y = std::sin (theta) * std::sin (phi);
      G4double z = std::cos (theta);
      fLightsVector = G4ThreeVector (x, y, z);
      vp.SetLightpointDirection(fLightsVector);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "Lights direction set to "
        << vp.GetLightpointDirection() << G4endl;
      }
    }
  }

  else if (command == fpCommandLightsVector) {
    fLightsVector = G4UIcommand::ConvertTo3Vector(newValue);
    vp.SetLightpointDirection(fLightsVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Lights direction set to "
      << vp.GetLightpointDirection() << G4endl;
    }
  }

  else if (command == fpCommandLineWidth) {
    if (verbosity >= G4VisManager::errors) {
      // A do-nothing command
      G4warn << command->GetGuidanceLine(0) << G4endl;
    }
  }

  else if (command == fpCommandLineSegments) {
    G4int nSides = G4UIcommand::ConvertToInt(newValue);
    nSides = vp.SetNoOfSides(nSides);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
      "Number of line segments per circle in polygon approximation is "
      << nSides << G4endl;
    }
  }

  else if (command == fpCommandNumberOfCloudPoints) {
    G4int nPoints = G4UIcommand::ConvertToInt(newValue);
    nPoints = vp.SetNumberOfCloudPoints(nPoints);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
      "Number of points to be used in cloud representation of volumes is "
      << nPoints << G4endl;
    }
  }

  else if (command == fpCommandPicking) {
    vp.SetPicking(G4UIcommand::ConvertToBool(newValue));
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Picking ";
      if (vp.IsPicking()) G4cout << "requested.";
      else G4cout << "inhibited.";
      G4cout << G4endl;
    }
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "You may need to issue \"/vis/viewer/update\"."
      << G4endl;
    }
  }

  else if (command == fpCommandProjection) {
    G4double fieldHalfAngle;
    const size_t iPos0 = 0;
    if (newValue[iPos0] == 'o') {  // "orthogonal"
      fieldHalfAngle = 0.;
    }
    else if (newValue[iPos0] == 'p') {  // "perspective"
      G4String dummy;
      G4String unit;
      std::istringstream is (newValue);
      is >> dummy >> fieldHalfAngle >> unit;
      fieldHalfAngle *= G4UIcommand::ValueOf(unit);
      if (fieldHalfAngle > 89.5 * deg || fieldHalfAngle <= 0.0) {
        if (verbosity >= G4VisManager::errors) {
          G4warn <<
          "ERROR: Field half angle should be 0 < angle <= 89.5 degrees.";
          G4warn << G4endl;
        }
        return;
      }
    }
    else {
      if (verbosity >= G4VisManager::errors) {
        G4warn << "ERROR: \"" << newValue << "\" not recognised."
        "  Looking for 'o' or 'p' first character." << G4endl;
      }
      return;
    }
    vp.SetFieldHalfAngle(fieldHalfAngle);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Projection style of viewer \"" << currentViewer->GetName()
      << "\" set to ";
      if (fieldHalfAngle == 0.) {
        G4cout << "orthogonal.";
      }
      else {
        G4cout << "perspective\n  with half angle " << fieldHalfAngle / deg
        << " degrees.";
      }
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandRotationStyle) {
    G4ViewParameters::RotationStyle style;
    if (newValue == "constrainUpDirection")
      style = G4ViewParameters::constrainUpDirection;
    else if (newValue == "freeRotation")
      style = G4ViewParameters::freeRotation;
    else {
      if (verbosity >= G4VisManager::errors) {
	G4warn << "ERROR: \"" << newValue << "\" not recognised." << G4endl;
      }
      return;
    }
    vp.SetRotationStyle(style);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Rotation style of viewer \"" << currentViewer->GetName()
      << "\" set to " << vp.GetRotationStyle()
      << G4endl;
    }
  }

  else if (command == fpCommandSectionPlane) {
    G4String choice, unit;
    G4double x, y, z, nx, ny, nz, F = 1.;
    std::istringstream is (newValue);
    is >> choice >> x >> y >> z >> unit >> nx >> ny >> nz;
    G4int iSelector = -1;
    if (G4StrUtil::icompare(choice, "off") == 0 ||
        !G4UIcommand::ConvertToBool(choice)) iSelector = 0;
    if (G4StrUtil::icompare(choice, "on") == 0 ||
        G4UIcommand::ConvertToBool(choice)) iSelector = 1;
    if (iSelector < 0) {
      if (verbosity >= G4VisManager::errors) {
        G4warn << "Choice not recognised (on/true or off/false)." << G4endl;
      }
      goto write_result;
    }
    // iSelector can only be 0 or 1
    switch (iSelector) {
      case 0:
        vp.UnsetSectionPlane();
        break;
      case 1:
        F = G4UIcommand::ValueOf(unit);
        x *= F; y *= F; z *= F;
        if (nx == 0. && ny == 0. && nz == 0.) {
          if (verbosity >= G4VisManager::errors) {
            G4warn << "Null normal." << G4endl;
          }
          break;;
        }
        // Make sure normal is normalised
        const G4Normal3D& normal = G4Normal3D(nx,ny,nz).unit();
        vp.SetSectionPlane(G4Plane3D(normal, G4Point3D(x,y,z)));
        vp.SetViewpointDirection(normal);
        break;
    }
  write_result:
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Section drawing is: ";
      if (vp.IsSection ()) G4cout << "on";
      else                 G4cout << "off";
      G4cout << ".\nSection plane is now: "
      << vp.GetSectionPlane ();
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandSpecialMeshRendering) {
    vp.SetSpecialMeshRendering(G4UIcommand::ConvertToBool(newValue));
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Special mesh rendering";
      if (vp.IsSpecialMeshRendering()) {
        G4cout << " requested. Current option is \""
        << vp.GetSpecialMeshRenderingOption() << "\" for ";
        if (vp.GetSpecialMeshVolumes().empty()) {
          G4cout << "any mesh.";
        } else{
          G4cout << "selected volumes:";
          for (const auto& pvNameCopyNo: vp.GetSpecialMeshVolumes()) {
            G4cout << "\n  " << pvNameCopyNo.GetName();
            if (pvNameCopyNo.GetCopyNo() >= 0) G4cout << ':' << pvNameCopyNo.GetCopyNo();
          }
        }
      }
      else G4cout << ": off.";
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandSpecialMeshRenderingOption) {
    G4ViewParameters::SMROption option = G4ViewParameters::meshAsDots;
    if (newValue == "surfaces") {
      option = G4ViewParameters::meshAsSurfaces;
    }
    vp.SetSpecialMeshRenderingOption(option);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Special mesh rendering option set to \""
      << vp.GetSpecialMeshRenderingOption() << "\"."
      << G4endl;
    }
  }

  else if (command == fpCommandSpecialMeshVolumes) {
    std::vector<G4ModelingParameters::PVNameCopyNo> requestedMeshes;
    if (newValue.empty()) {
      vp.SetSpecialMeshVolumes(requestedMeshes);  // Empty list
    } else {
      // Algorithm from Josuttis p.476.
      G4String::size_type iBegin, iEnd;
      iBegin = newValue.find_first_not_of(' ');
      while (iBegin != G4String::npos) {
	iEnd = newValue.find_first_of(' ',iBegin);
	if (iEnd == G4String::npos) {
	  iEnd = newValue.length();
	}
	G4String name(newValue.substr(iBegin,iEnd-iBegin));
	iBegin = newValue.find_first_not_of(' ',iEnd);
	if (iBegin == G4String::npos) {
	  if (verbosity >= G4VisManager::warnings) {
	    G4warn <<
	    "WARNING: G4VisCommandsViewerSet::SetNewValue: /vis/viewer/set/specialMeshVolumes"
	    "\n  A pair not found.  (There should be an even number of parameters.)"
	    "\n  Command ignored."
	    << G4endl;
	    return;
	  }
	}
	iEnd = newValue.find_first_of(' ',iBegin);
	if (iEnd == G4String::npos) {
	  iEnd = newValue.length();
	}
	G4int copyNo;
	std::istringstream iss(newValue.substr(iBegin,iEnd-iBegin));
	if (!(iss >> copyNo)) {
	  if (verbosity >= G4VisManager::warnings) {
	    G4warn <<
	    "WARNING: G4VisCommandsViewerSet::SetNewValue: /vis/viewer/set/specialMeshVolumes"
	    "\n  Error reading copy number - it was not numeric?"
	    "\n  Command ignored."
	    << G4endl;
	    return;
	  }
	}
	requestedMeshes.push_back(G4ModelingParameters::PVNameCopyNo(name,copyNo));
	iBegin = newValue.find_first_not_of(' ',iEnd);
      }
      vp.SetSpecialMeshVolumes(requestedMeshes);
    }
    if (verbosity >= G4VisManager::confirmations) {
      if (vp.GetSpecialMeshVolumes().empty()) {
	G4cout <<
	"Special mesh list empty, which means \"all meshes\"."
	<< G4endl;
      } else {
	G4cout << "Selected special mesh volumes are:";
	for (const auto& pvNameCopyNo: vp.GetSpecialMeshVolumes()) {
	  G4cout << "\n  " << pvNameCopyNo.GetName();
	  if (pvNameCopyNo.GetCopyNo() >= 0) G4cout << ':' << pvNameCopyNo.GetCopyNo();
	}
	G4cout << G4endl;
      }
    }
  }

  else if (command == fpCommandStyle) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    const size_t iPos0 = 0;
    if (newValue[iPos0] == 'w') {  // "wireframe"
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          break;
        case G4ViewParameters::hlr:
          break;
        case G4ViewParameters::hsr:
          vp.SetDrawingStyle(G4ViewParameters::wireframe);
          break;
        case G4ViewParameters::hlhsr:
          vp.SetDrawingStyle(G4ViewParameters::hlr);
          break;
        case G4ViewParameters::cloud:
          vp.SetDrawingStyle(G4ViewParameters::wireframe);
          break;
      }
    }
    else if (newValue[iPos0] == 's') {  // "surface"
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          vp.SetDrawingStyle(G4ViewParameters::hsr);
          break;
        case G4ViewParameters::hlr:
          vp.SetDrawingStyle(G4ViewParameters::hlhsr);
          break;
        case G4ViewParameters::hsr:
          break;
        case G4ViewParameters::hlhsr:
          break;
        case G4ViewParameters::cloud:
          vp.SetDrawingStyle(G4ViewParameters::hsr);
          break;
     }
    }
    else if (newValue[iPos0] == 'c') {  // "cloud"
      switch (existingStyle) {
        case G4ViewParameters::wireframe:
          vp.SetDrawingStyle(G4ViewParameters::cloud);
          break;
        case G4ViewParameters::hlr:
          vp.SetDrawingStyle(G4ViewParameters::cloud);
          break;
        case G4ViewParameters::hsr:
          vp.SetDrawingStyle(G4ViewParameters::cloud);
          break;
        case G4ViewParameters::hlhsr:
          vp.SetDrawingStyle(G4ViewParameters::cloud);
          break;
        case G4ViewParameters::cloud:
          break;
      }
    }
    else {
      if (verbosity >= G4VisManager::errors) {
        G4warn << "ERROR: \"" << newValue << "\" not recognised."
        "  Looking for 'w' or 's' or 'c' first character." << G4endl;
      }
      return;
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
      << "\" set to " << vp.GetDrawingStyle()
      << G4endl;
    }
  }

  else if (command == fpCommandTargetPoint) {
    G4ThreeVector targetPoint =
    fpCommandTargetPoint->GetNew3VectorValue(newValue);
    const G4Point3D& standardTargetPoint =
    currentViewer->GetSceneHandler()->GetScene()->GetStandardTargetPoint();
    vp.SetCurrentTargetPoint(targetPoint - standardTargetPoint);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Target point set to "
      << fpCommandTargetPoint->ConvertToStringWithBestUnit
      (targetPoint)
      << "\n\"Current Target Point\" set to  "
      << fpCommandTargetPoint->ConvertToStringWithBestUnit
      (vp.GetCurrentTargetPoint())
      << "\n\"Standard Target Point\" is "
      << fpCommandTargetPoint->ConvertToStringWithBestUnit
      (standardTargetPoint)
      << G4endl;
    }
  }

  else if (command == fpCommandUpThetaPhi) {
    G4double theta, phi;
    if (ConvertToDoublePair(newValue, theta, phi)) {
      G4double x = std::sin (theta) * std::cos (phi);
      G4double y = std::sin (theta) * std::sin (phi);
      G4double z = std::cos (theta);
      fUpVector = G4ThreeVector (x, y, z);
      vp.SetUpVector(fUpVector);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "Up direction set to " << vp.GetUpVector() << G4endl;
      }
    }
  }

  else if (command == fpCommandUpVector) {
    fUpVector = G4UIcommand::ConvertTo3Vector(newValue).unit();
    vp.SetUpVector(fUpVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Up direction set to " << vp.GetUpVector() << G4endl;
    }
  }

  else if (command == fpCommandViewpointThetaPhi) {
    G4double theta, phi;
    if (ConvertToDoublePair(newValue, theta, phi)) {
      G4double x = std::sin (theta) * std::cos (phi);
      G4double y = std::sin (theta) * std::sin (phi);
      G4double z = std::cos (theta);
      fViewpointVector = G4ThreeVector (x, y, z);
      vp.SetViewAndLights(fViewpointVector);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "Viewpoint direction set to "
        << vp.GetViewpointDirection() << G4endl;
        if (vp.GetLightsMoveWithCamera ()) {
          G4cout << "Lightpoint direction set to "
          << vp.GetActualLightpointDirection () << G4endl;
        }
      }
    }
  }

  else if (command == fpCommandViewpointVector) {
    G4ThreeVector viewpointVector = G4UIcommand::ConvertTo3Vector(newValue);
    if (viewpointVector.mag2() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4warn << "ERROR: Null viewpoint vector. No action taken." << G4endl;
      }
    } else {
      fViewpointVector = viewpointVector.unit();
      vp.SetViewAndLights(fViewpointVector);
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "Viewpoint direction set to "
        << vp.GetViewpointDirection() << G4endl;
        if (vp.GetLightsMoveWithCamera ()) {
          G4cout << "Lightpoint direction set to "
          << vp.GetActualLightpointDirection () << G4endl;
        }
      }
    }
  }

  else if (command == fpCommandTimeWindowDisplayHeadTime)
  {
    G4String display;
    G4double screenX, screenY, screenSize, red, green, blue;
    std::istringstream iss(newValue);
    iss >> display >> screenX >> screenY
    >> screenSize >> red >> green >> blue;
    vp.SetDisplayHeadTime(command->ConvertToBool(display));
    vp.SetDisplayHeadTimeX(screenX);
    vp.SetDisplayHeadTimeY(screenY);
    vp.SetDisplayHeadTimeSize(screenSize);
    vp.SetDisplayHeadTimeRed(red);
    vp.SetDisplayHeadTimeGreen(green);
    vp.SetDisplayHeadTimeBlue(blue);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Display head time flag set: "
      << vp
      << G4endl;
    }
  }

  else if (command == fpCommandTimeWindowDisplayLightFront)
  {
    G4String display, originX, originY, originZ, unitS, originT, unitT;
    G4double red, green, blue;
    std::istringstream iss(newValue);
    iss >> display
    >> originX >> originY >> originZ >> unitS
    >> originT >> unitT
    >> red >> green >> blue;
    vp.SetDisplayLightFront(command->ConvertToBool(display));
    vp.SetDisplayLightFrontX
    (command->ConvertToDimensionedDouble(G4String(originX + ' ' + unitS)));
    vp.SetDisplayLightFrontY
    (command->ConvertToDimensionedDouble(G4String(originY + ' ' + unitS)));
    vp.SetDisplayLightFrontZ
    (command->ConvertToDimensionedDouble(G4String(originZ + ' ' + unitS)));
    vp.SetDisplayLightFrontT
    (command->ConvertToDimensionedDouble(G4String(originT + ' ' + unitT)));
    vp.SetDisplayLightFrontRed(red);
    vp.SetDisplayLightFrontGreen(green);
    vp.SetDisplayLightFrontBlue(blue);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Display light front flag set: "
      << vp
      << G4endl;
    }
  }

  else if (command == fpCommandTimeWindowEndTime)
  {
    G4String end_time_string, end_time_unit,
    time_range_string, time_range_unit;
    std::istringstream iss(newValue);
    iss >> end_time_string >> end_time_unit
    >> time_range_string >> time_range_unit;
    vp.SetEndTime
    (command->ConvertToDimensionedDouble
     (G4String(end_time_string + ' ' + end_time_unit)));
    G4double timeRange = command->ConvertToDimensionedDouble
    (G4String(time_range_string + ' ' + time_range_unit));
    if (timeRange > 0.) {
      vp.SetStartTime
      (vp.GetEndTime() - timeRange);
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout
      << "Time window start time: " << vp.GetStartTime()/ns << " ns"
      << ", time window end time: " << vp.GetEndTime()/ns << " ns";
      if (timeRange > 0.) {
        G4cout << "\n  (time range: " << timeRange/ns << " ns)";
      }
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandTimeWindowFadeFactor) {
    vp.SetFadeFactor(command->ConvertToDouble(newValue));
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Time window fade factor changed to " << vp.GetFadeFactor()
      << G4endl;
    }
  }

  else if (command == fpCommandTimeWindowStartTime)
  {
    G4String start_time_string, start_time_unit,
    time_range_string, time_range_unit;
    std::istringstream iss(newValue);
    iss >> start_time_string >> start_time_unit
    >> time_range_string >> time_range_unit;
    vp.SetStartTime
    (command->ConvertToDimensionedDouble
     (G4String(start_time_string + ' ' + start_time_unit)));
    G4double timeRange = command->ConvertToDimensionedDouble
    (G4String(time_range_string + ' ' + time_range_unit));
    if (timeRange > 0.) {
      vp.SetEndTime
      (vp.GetStartTime() + timeRange);
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout
      << "Time window start time: " << vp.GetStartTime()/ns << " ns"
      << ", time window end time: " << vp.GetEndTime()/ns << " ns";
      if (timeRange > 0.) {
        G4cout << "\n  (time range: " << timeRange/ns << " ns)";
      }
      G4cout << G4endl;
    }
  }

  else {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: G4VisCommandsViewerSet::SetNewValue: unrecognised command."
      << G4endl;
    }
    return;
  }

  SetViewParameters(currentViewer,vp);
}
