// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayView.cc,v 2.9 1998/11/06 13:41:56 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing view.

#include "G4RayView.hh"

#include "G4RayScene.hh"

#include "G4ParticleGun.hh"
#include "G4TransportationManager.hh"
#include "geomdefs.hh"

#include "G4VSolid.hh"

#include <math.h>

#ifdef GNU_GCC
  template class RWTValHashDictionary <G4RayCoordinate, G4RayHitColour>;
  template class RWTValHashDictionary <G4RayHitColour, G4int>;
#endif

G4RayView::G4RayView (G4RayScene& scene, const G4String& name):
G4VView (scene, -1, name),
fScene (scene),
fpWorld (0)
{
 //RUN MANAGER
  fG4RayRunManager         = new G4RunManager;

  //G4RaySteppingAction
   fG4RaySteppingAction  = new G4RaySteppingAction;

  //Gun
   fG4RayGenerator         = new G4RayPrimaryGeneratorAction;

   fG4RayRunManager -> SetUserAction ( fG4RayGenerator );
   fG4RayRunManager -> SetUserAction ( fG4RaySteppingAction);

   fG4RayGun                = fG4RayGenerator -> GetParticleGun(); 
}

G4RayView::~G4RayView(){

  delete fG4RayRunManager;
  // Note: G4RunManager deletes userPrimaryGeneratorAction, what we
  // have called fG4RayGenerator).  And G4SteppingManager deletes
  // fUserSteppingAction, what we have called fG4RaySteppingAction

}

void G4RayView::SetView () {

  fpWorld = G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();

  // We need to know the world so we can be sure to start traceing inside.
  
  if (!fpWorld) {
    G4Exception ("G4RayView::SetView: no world established yet!");
  }
  
  // VECTORS

  const G4Vector3D& viewvector = fVP.GetViewpointDirection().unit();
  // view vector on the line from the target point to the eye-- view
  // point direction

  const G4Vector3D& upvector = fVP.GetUpVector().unit();   // the up vector of the volume 

  const G4Vector3D& lightDirection = fVP.GetActualLightpointDirection().unit();
  // vector of incident light

 
  // DOUBLES
  G4double radius = fScene.GetSceneData ().GetExtent ().GetExtentRadius ();
  //radius of the bounding sphere.

  G4double cameradistance = fVP.GetCameraDistance (radius);
  // from the center of the bounding sphere to the camera.

  G4double johnsneardistance = fVP.GetNearDistance (cameradistance, radius);

  G4double halfheight = fVP.GetFrontHalfHeight (johnsneardistance, radius);
  // half the size of the front view plane (screen)

  // xprime ,yprime, zprime define the screen coordinate system.
  G4Vector3D zprime       = viewvector.unit();
  G4Vector3D xprime       = (upvector.cross(viewvector)).unit();
  G4Vector3D yprime       = (zprime.cross(xprime)).unit();
  // Make xprime and yprime have lengths equal to halfheight so that
  // they can be multplied by a screen coordinate in the range -1 to
  // +1.
  fXprime = halfheight * xprime;
  fYprime = halfheight * yprime;
 

  // POINTS
  const G4Point3D& centeroftarget = fVP.GetCurrentTargetPoint ();// the center of target
 

  //fCenterofscreen = centeroftarget + radius * viewvector.unit();
  fCamerapoint  = centeroftarget + cameradistance * viewvector.unit();
  fCenterofscreen = fCamerapoint - johnsneardistance * viewvector.unit();

  G4cout << " The viewvector is : " << viewvector<<endl;
  G4cout << " The lightDirection is : " << lightDirection << endl;
  G4cout << " The cameraDistance is: " << cameradistance << endl;
  G4cout << "  fCamerapoint is  : " << fCamerapoint << endl;
  G4cout << "  fCenterofscreen is : " << fCenterofscreen << endl;
  G4cout << "  centeroftarget is : " << centeroftarget << endl;

}

G4ViewParameters G4RayView::fImageViewParameters;
G4double         G4RayView::fInitialResolution      = 27.;
G4double         G4RayView::fRequestedResolution    = 27.;
G4double         G4RayView::fCurrentResolution      = 27.;
G4float          G4RayView::fMaxMinColourSeparation = 0.;
G4int            G4RayView::fImageWidth             = 0;
G4int            G4RayView::fImageHeight            = 0;
G4int            G4RayView::fMaxImageMapEntries     = 0;
G4bool           G4RayView::fImageMapRejigged       = true;
G4bool           G4RayView::fTriggerRescan          = true;
G4int            G4RayView::fWaitLimit              = 200;
G4int            G4RayView::fWaitCount              = 200;

RWTValOrderedVector <G4RayColourCell>
G4RayView::fRayHitColourMap (256);

RWTValHashDictionary <G4RayCoordinate, G4RayHitColour>
G4RayView::fRayHits (G4RayView::RayCoordHash, 1000);

RWTValHashDictionary <G4RayHitColour, G4int>
G4RayView::fRayHitColourIndex (G4RayView::RayColourHash, 1000);

void G4RayView::InitialiseImage () {
  fRayHitColourMap.clear ();
  fRayHits.clear ();
  fMaxMinColourSeparation = 0.;
  fInitialResolution = 27;
  fCurrentResolution = fInitialResolution;
  fRequestedResolution = fInitialResolution;
  G4cout << "Requested resolution reset to " << fRequestedResolution << endl;
  fWaitLimit = 200;
  G4cout << "Wait limit reset to: " << fWaitLimit
       << "\n(e.g, display refreshed every this many rays traced.)"
       << endl;
  TriggerRedraw ();
}

void G4RayView::TriggerRedraw () {
  fImageMapRejigged = true;
  fTriggerRescan = true;
  fWaitCount = fWaitLimit;
}

G4bool G4RayView::Scan (int& i, int& j, int& k,
			float& red, float& green, float& blue,
			int& iPixVal) {

  static G4bool Done = true;
  static G4bool nearlyDone;
  static G4double cellSize;
  static G4int requestedNoOfRays;
  static G4int doneCount;
  static G4int u, v;


  if (Done) {  // Assume called anew.
    if (fRequestedResolution < 1.0) {
      G4cout <<
	"Over-sampling not implemented - resetting requested resolution to 1."
	<< endl;
      fRequestedResolution = 1.;
    }
    Done = false;
    nearlyDone = false;
    requestedNoOfRays = (int) 
      ((((fImageWidth - 1) / fRequestedResolution) + 1) *
       (((fImageHeight - 1) / fRequestedResolution) + 1));
    doneCount = 1;
    cellSize = fInitialResolution;
    u = v = 0;
    if (fCurrentResolution >= 3 * fRequestedResolution) {
      G4cout << "Resizing Ray Hits Dictionaries - " << requestedNoOfRays
	   << " buckets." << endl;
      fRayHits.resize (requestedNoOfRays);
      fRayHitColourIndex.resize (requestedNoOfRays);
      G4cout << "Continuing..." << endl;
    }
  }

  if (fTriggerRescan) {
    fTriggerRescan = false;
    cellSize = 3 * cellSize;
    u = v = 0;
  }

  if (fRayHitColourMap.entries ()  >= fMaxImageMapEntries) {
    static int count = 0;
    if (count++ < 10 || count % 10 == 0) {
      G4cout << "G4RayView::Scan: compacting colour map - count = "
	   << count << endl;
    }
    CompactColour (false);
    if (fRayHitColourMap.entries () >= fMaxImageMapEntries) {
      G4Exception ("G4RayView::Scan: compaction failure!");
    }
    fImageMapRejigged = true;
  }

  i = u;
  j = v;
  k = int (cellSize);

  G4double x = 2. * (G4double (i + cellSize / 2) / fImageWidth - 0.5);
  G4double y = 2. * (0.5 - G4double (j + cellSize / 2) / fImageHeight);
  if (x >  1.0) x =  1.0;
  if (y < -1.0) y = -1.0;

  G4RayCoordinate rayCoord (x, y);
  G4RayHitColour rayHitColour;
  if (fRayHits.findValue (rayCoord, rayHitColour)) {
    iPixVal = fRayHitColourIndex [rayHitColour];  // Old pixval.
    red   = rayHitColour.GetRed   ();
    green = rayHitColour.GetGreen ();
    blue  = rayHitColour.GetBlue  ();
    //G4cout << "G4RayView::Scan: " << x << ", " << y << " already traced."
    //  "  Pixval = " << iPixVal
    //     << endl;
  }
  else {
    // Initiate a new ray.
    Trace (x, y, red, green, blue);
    doneCount++;
    int ii;
    for (ii = 0; ii < fRayHitColourMap.entries (); ii++) {
      G4RayColourCell& colourCell = fRayHitColourMap [ii];
      if (ColourSeparation
	  (red, green, blue,
	   colourCell.red, colourCell.green, colourCell.blue) < 0.01) {
	break;  // ...if there's a close-match.
      }
    }
    iPixVal = ii; // If no close-match found, ii = no. of entries, which
                  // will be no. of entries - 1 after append.
    if (iPixVal == fRayHitColourMap.entries ())
      fRayHitColourMap.append (G4RayColourCell (red, green, blue, 1));
    rayHitColour = G4RayHitColour (red, green, blue);
    fRayHits.insertKeyAndValue (rayCoord, rayHitColour);
    fRayHitColourIndex.insertKeyAndValue (rayHitColour, iPixVal);
  }

  //G4cout << "G4RayView::Scan: Draw x,y,r,g,b,i,j,k,iPixVal: "
  //     << x << ' ' << y << ' ' << red << ' ' << green << ' ' << blue << ' '
  //     << i << ' ' << j << ' ' << k << ' ' << iPixVal 
  //     << endl;

  if (!(doneCount % fWaitLimit) && !nearlyDone) {
    doneCount = 1;
    G4cout << "G4RayView::Scan: "
	 << (int) ((100. * (fImageWidth * v + u)) /
		   (fImageWidth * fImageHeight))
	 << "% Done at cell-size "
	 << cellSize
	 << endl;
  }

  u += int (cellSize);

  if (u >= fImageWidth) {
    u = 0;
    v += int (cellSize);
  }

  if (v >= fImageHeight) {
    u = v = 0;
    cellSize = cellSize / 3;  // Next size down.
  }

  if ((cellSize < fRequestedResolution
      || 3 * cellSize < 1.)
      // The above line can be removed when over-sampling implemented.
      && !nearlyDone
      ) {
    G4cout << "G4RayView::Scan: final redraw." << endl;
    nearlyDone = true;
    cellSize = 3 * cellSize;
    u = v = 0;
    TriggerRedraw ();
  }

  if ((cellSize < fRequestedResolution
      || 3 * cellSize < 1.)
      // The above line can be removed when over-sampling implemented.
      && nearlyDone
      ) {
    fCurrentResolution = 3 * cellSize;
    G4cout << "G4RayView::Scan: requested resolution = "
	 << fRequestedResolution
	 << ", achieved\nresolution = "
	 << fCurrentResolution << " (always a Power of 3).";
    G4cout << "\nNo. of ray hits          = " << fRayHits.entries ();
    G4cout << "\nNo. of ray Hit colours   = " << fRayHitColourIndex.entries ();
    G4cout << "\nNo. of colourmap entries = " << fRayHitColourMap.entries ();
    fWaitLimit = 3 * fWaitLimit;
    G4cout << "\nWait limit increased to: " << fWaitLimit
	 << "\n(e.g, display refreshed every this many rays traced.)";
    if (fRequestedResolution < 1.) {
      G4cout << "\nOver-sampling (resolution < 1.) not implemented.";
    }
    G4cout << "\nYou can continue as follows:"
      "\n/vis~/ray/resolution " << cellSize
	 << "\n/vis~/draw/current"
      "\n(It may take some time - 10 Times longer than so far!!!"
      "\nMoreover, it will use about another "
	 << 30 * (((fImageWidth - 1) / cellSize) + 1) *
      (((fImageHeight - 1) / cellSize) + 1) * 1.e-6
	 << " MB of memory to Store the image.)";
    G4cout << endl;
    Done = true;
  }

  return !Done;
}

void G4RayView::CompactColour (G4bool Print) {

  G4int i;
  G4int j;

  RWTValOrderedVector <G4RayColourCell>& cmap = fRayHitColourMap;

  G4double minColourSeparation = 2.;
  // i.e, bigger than possible in the colour cube.
  for (i = 0; i < cmap.entries (); i++) {
    G4RayColourCell& celli = fRayHitColourMap [i];
    for (j = i + 1; j < cmap.entries (); j++) {
      G4RayColourCell& cellj = fRayHitColourMap [j];
      float colourSeparation = ColourSeparation
	(celli.red, celli.green, celli.blue,
	 cellj.red, cellj.green, cellj.blue);
      if (minColourSeparation > colourSeparation) {
	minColourSeparation = colourSeparation;
      }
    }
  }

  if (fMaxMinColourSeparation < minColourSeparation) {
    fMaxMinColourSeparation = minColourSeparation;
    G4cout << "G4RayView::CompactColour: minColourSeparation has peaked at "
	 << minColourSeparation << endl;
  }

  if (minColourSeparation > 0.125) {
    G4cerr << "G4RayView::CompactColour: WARNING: minimum colour separation"
      "\nis greater than 0.125, the separation of R/G/B 3/3/2!" << endl;
  }

  G4int foundi;
  G4int foundij [2];
  RWTValOrderedVector <G4int> foundjs;
  G4int jj;

  // Process all pairs with colour separation < compaction * minimum.
  const float compaction = 1.5;
  do {
    foundi = -1;
    for (i = 0; i < cmap.entries (); i++) {
      G4RayColourCell& celli = fRayHitColourMap [i];
      for (j = i + 1; j < cmap.entries (); j++) {
	G4RayColourCell& cellj = fRayHitColourMap [j];
	float colourSeparation = ColourSeparation
	  (celli.red, celli.green, celli.blue,
	   cellj.red, cellj.green, cellj.blue);
	if (colourSeparation < compaction * minColourSeparation) {
	  if (foundi < 0) {
	    foundi = i;
	    foundjs.clear ();
	    foundjs.append (i);  // First entry is the i value.
	  }
	  foundjs.append (j);
	}
      }
      if (foundi >= 0) break;
    }

    if (foundi >= 0) {

      if (Print) {
	G4cout << "G4RayView::CompactColour: foundjs:";
	for (jj = 0; jj < foundjs.entries (); jj++) {
	  G4cout << ' ' << foundjs [jj];
	}
	G4int iii = 0;
	G4cout << "\nG4RayView::CompactColour: (colour, index) before:";
	fRayHitColourIndex.applyToKeyAndValue (PrintPixVals, (void *) &iii);
      }

      // Change all ray hits with pix val j to pix val i.
      fRayHitColourIndex.applyToKeyAndValue (ChangePixVals, (void *) &foundjs);
      // Decrement pix vals of all ray hits with pix val j + 1  and and above.
      for (jj = foundjs.entries (); jj-- > 1;) { // ...in reverse order so
	// that decrements are applied in the correct order.
	G4int jjj = foundjs [jj];
	fRayHitColourIndex.applyToKeyAndValue
	  (DecrementPixVals, (void *) &jjj);
      }
      
      // Recalculate colour map and remove redundant entries.
      G4RayColourCell& celli = fRayHitColourMap [foundi];
      for (jj = foundjs.entries (); jj-- > 1;) { // ...in reverse order so
	// that we can safely remove colour cells.
	G4int jjj = foundjs [jj];
	G4RayColourCell& cellj = fRayHitColourMap [jjj];
	// Recompute colour of cell.
	G4int N = (celli.nHits + cellj.nHits);
	celli.red   = (celli.nHits * celli.red
		       + cellj.nHits * cellj.red)   / N;
	celli.green = (celli.nHits * celli.green
		       + cellj.nHits * cellj.green) / N;
	celli.blue  = (celli.nHits * celli.blue
		       + cellj.nHits * cellj.blue)  / N;
	celli.nHits = N;
	// Finally, remove the redundant cells.
	if (Print) G4cout << "Removing " << jjj << endl;
	fRayHitColourMap.removeAt (jjj);
      }

      if (Print) {
	G4int iii = 0;
	G4cout << "\nG4RayView::CompactColour: (colour, index) after:";
	fRayHitColourIndex.applyToKeyAndValue (PrintPixVals, (void *) &iii);
	G4cout << endl;
      }

    }
  } while (foundi >=0);

  if (Print) {
    G4cout << "G4RayView::CompactColour: map compressed to "
	 << fRayHitColourMap.entries ()
	 << endl;
  }
}

float G4RayView::ColourSeparation (float ri, float gi, float bi,
				   float rj, float gj, float bj) {
  float sep = abs (rj -ri);
  float sepg = abs (gj - gi);
  if (sep < sepg) sep = sepg;
  float sepb = abs (bj - bi) / 2.;
  if (sep < sepb) sep = sepb;
  return sep;
}

void G4RayView::ChangePixVals (const G4RayHitColour& colour,
			       G4int& index, void* pfjs) {

  RWTValOrderedVector <G4int>* pfoundjs = (RWTValOrderedVector <G4int>*) pfjs;
  G4int i = (*pfoundjs) [0];

  // Change all ray Hit colours with index j to index i.
  for (int jj = 1; jj < (*pfoundjs).entries (); jj++) {
    G4int j = (*pfoundjs) [jj];
    if (index == j) {
      //G4cout << "G4RayView::ChangePixVal: ray Hit colour index changed from "
      //     << index;
      index = i;
      //G4cout << " to " << index << endl;
    }
  }
}

void G4RayView::DecrementPixVals (const G4RayHitColour& colour,
				  G4int& index, void* pjjj) {
  G4int jjj = * (G4int*) pjjj;
  // Decrement indices all ray Hit colours with index jjj + 1 and above.
  if (index > jjj) {
    //G4cout << index  << " decremented." << endl;
    index--;
  }
}

void G4RayView::PrintPixVals (const G4RayHitColour& colour,
				  G4int& index, void* pjjj) {
  G4cout << ' '
    //<< '(' << colour.colour << ','
       << index
    //<< ')'
    ;
}

void G4RayView::Trace (G4double x, G4double y,
		       float& red, float& green, float& blue) {

  // Colour
  G4Colour colour;           // The colour of the object
  G4double Ip             = 0.6 ; // incident light intensity
  G4double Id             = 0.9 ;  // reflected light intensity
  G4double Ia             = 0.2;  // ambient light intensity
  G4double I;                      // total intensity
  G4double intrinsicred;           // intrinsic colour of the object
  G4double intrinsicgreen;         // given in red green 
  G4double intrinsicblue;          // and blue
  G4double outsideDistance;


  G4Point3D        newPoint;

  G4double Dot;               // Dot product of the normal vector of surface
                            // with the normalized vector of incident light 
  G4double reflectionDot;     // Dot product of the normal vector of reflection
                            // with the normal vector of angle of viewing
  G4Vector3D  rayDirection;      //unit vector, the ray direction.
  G4Point3D  rayPosition;    // position of ray on the front view plane

  G4Vector3D roughrayDirection; // the vector from camera to 
                                //the point on screen
  G4Vector3D reflightDirection; // vector of reflection light
  G4Vector3D screenvector; // vector on the front view plane pointing from
                         // the center of screen to a given point on screen 

  const G4Vector3D& viewvector = fVP.GetViewpointDirection().unit();
  // view vector on the line from the target point to the eye-- view
  // point direction
 
  const G4Vector3D& lightDirection = fVP.GetActualLightpointDirection().unit();
  // vector of incident light

 

  EInside    whereisit;  
 
  screenvector    = fXprime * x  + fYprime * y;
  rayPosition     = fCenterofscreen +  screenvector;

  if (fVP.GetFieldHalfAngle() == 0 ) {
     rayDirection = - viewvector.unit();
  }
  else {  
  roughrayDirection = rayPosition - fCamerapoint;
  rayDirection      = roughrayDirection.unit();
  }
   
  // Default - black.
  I = 0;

  // Check if rayPosition is in the world.
  outsideDistance = kInfinity;
  whereisit = fpWorld -> GetLogicalVolume() -> GetSolid()->Inside(rayPosition);
  if (whereisit == kInside) {
    outsideDistance = 0.;
  }
  else {  // It's outside the world, so move it inside.
    outsideDistance = fpWorld -> GetLogicalVolume() ->
      GetSolid()-> DistanceToIn(rayPosition,rayDirection);  
    if ( outsideDistance != kInfinity) {
      rayPosition = rayPosition + outsideDistance * rayDirection;
    }
  }

  if (outsideDistance != kInfinity) {

    fG4RayGun -> SetParticlePosition(rayPosition);
    fG4RayGun -> SetParticleMomentumDirection(rayDirection);

    fG4RayRunManager -> BeamOn(1);

    if (fG4RaySteppingAction -> WeHaveContact ()) {

      // Store these in the scene so we can use culling algorithms.
      // Not figured out how to do this yet.
      // Certainly we don't need to Store these in the scene now
      // culling is Done in the modeling category.
      // Rethink this.
      //fScene.SetWorkingPhysicalVolume
      //(fG4RaySteppingAction -> GetPhysicalVolume ());
      //fScene.SetWorkingLogicalVolume
      //(fG4RaySteppingAction -> GetLogicalVolume ());
      //fScene.SetWorkingSolid (fG4RaySteppingAction -> GetSolid ());
      const G4VisAttributes* pVA =
	fG4RaySteppingAction -> GetVisAttributes();
      pVA = GetApplicableVisAttributes (pVA);
      //fScene.SetWorkingVisAttributes (pVA);

      /**** if (!fScene.IsThisCulled ()) ****/ {

	G4Normal3D normal = fG4RaySteppingAction -> GetSurfaceNormal ();

	//  Return colour of the object
	colour = pVA -> GetColour();
	intrinsicred   = colour.GetRed();
	intrinsicgreen = colour.GetGreen();
	intrinsicblue  = colour.GetBlue();
	
	//G4cout << "The intrnsic red colour is: " << red << endl;     
	//G4cout << "The intrinsic green colour is: " << green << endl;     
	//G4cout << "The intrinsic blue colour is: " << blue << endl;     
	
	// Calculate Intencity.
	Dot               = normal.unit() * lightDirection.unit();
	reflightDirection = 2.0 * Dot * normal.unit() - lightDirection.unit();
	reflectionDot     = viewvector.unit() * reflightDirection.unit();
        
	if (Dot < 0) Dot = 0;
	if (reflectionDot < 0) reflectionDot = 0;
	I = Ia + Ip * Dot ; //+ Id * pow(reflectionDot,90);
	
      }
    }
 }
  /*************
 else {
   
    outsideDistance = fpWorld -> GetLogicalVolume() -> GetSolid()-> DistanceToIn(rayPosition,rayDirection);  
    if ( outsideDistance == kInfinity) {I=0;}
                                        
    else{
       newPoint = rayPosition + outsideDistance * rayDirection;
     
       fG4RayGun -> SetParticlePosition(newPoint);
       fG4RayGun -> SetParticleMomentumDirection(rayDirection);

       fG4RayRunManager -> BeamOn(1);

       if (fG4RaySteppingAction -> WeHaveContact ()) {

	 // Store these in the scene so we can use culling algorithms.
	 fScene.SetWorkingPhysicalVolume
	   (fG4RaySteppingAction -> GetPhysicalVolume ());
	 fScene.SetWorkingLogicalVolume
	   (fG4RaySteppingAction -> GetLogicalVolume ());
	 fScene.SetWorkingSolid (fG4RaySteppingAction -> GetSolid ());
	 const G4VisAttributes* pVA =
	   fG4RaySteppingAction -> GetVisAttributes();
	 pVA = GetApplicableVisAttributes (pVA);
	 fScene.SetWorkingVisAttributes (pVA);

	 if (!fScene.IsThisCulled ()) {

	   G4Normal3D normal = fG4RaySteppingAction -> GetSurfaceNormal ();

	   //  Return colour of the object
	   colour = pVA -> GetColour();
	   intrinsicred   = colour.GetRed();
	   intrinsicgreen = colour.GetGreen();
	   intrinsicblue  = colour.GetBlue();

	   //G4cout << "The intrnsic red colour is: " << red << endl;     
	   //G4cout << "The intrinsic green colour is: " << green << endl;     
	   //G4cout << "The intrinsic blue colour is: " << blue << endl;     

	   // Calculate Intencity.
	   Dot               = normal.unit() * lightDirection.unit();
	   reflightDirection =
	     2.0 * Dot * normal.unit() - lightDirection.unit();
	   reflectionDot     = viewvector.unit() * reflightDirection.unit();
        
	   if (Dot < 0) Dot = 0;
	   if (reflectionDot < 0) reflectionDot = 0;
	   I = Ia + Ip * Dot ; //+ Id * pow(reflectionDot,90);
	 }
       }
    }
 }
 **********************/

  //  G4cout << " Intensity : " << I << endl;
 
  red   = intrinsicred * I ;
  green = intrinsicgreen * I ;
  blue  = intrinsicblue * I;
  if (I > 0.) {
    G4double specular = Id * pow(reflectionDot,10);
    red   += specular;
    green += specular;
    blue  += specular;
  }

  /****************************************************************

  ////////////////////////////////////////////////////////
  // A very testing pallette...
  ////////////////////////////////////////////////////////

  if (x < -0.8) {
    red   = (y + 1.) / 2.;
    green = 0;
    blue  = 0;
  }
  else if (x < -0.6) {
    red   = 0;
    green = (y + 1.) / 2.;
    blue  = 0;
  }
  else if (x < -0.4) {
    red   = 0;
    green = 0;
    blue  = (y + 1.) / 2.;
  }
  else if (x < -0.2) {
    red   = 0;
    green = (y + 1.) / 2.;
    blue  = (y + 1.) / 2.;
  }
  else if (x < 0.) {
    red   = (y + 1.) / 2.;
    green = 0;
    blue  = (y + 1.) / 2.;
  }
  else if (x < 0.2) {
    red   = (y + 1.) / 2.;
    green = (y + 1.) / 2.;
    blue  = 0.;
  }
  else if (x < 0.4) {
    red   = (y + 1.) / 2.;
    green = (y + 1.) / 2.;
    blue  = (y + 1.) / 2.;
  }
  else {
    red   = pow (sin (y), 2);
    green = pow (cos (y), 2);
    blue  = x * x;
  }

  ***********************************/
  /**********************************

  ////////////////////////////////////////////////////////
  // A shiny red sphere...
  ////////////////////////////////////////////////////////

  const G4double size (0.5);
  const G4Point3D centre;
  const G4double ambient (0.2);
  const G4Vector3D viewpointDirection (0.0, 0.0, 1.0);
  const G4Vector3D lightpointDirection = G4Vector3D (1.0, 1.0, 1.0).unit ();
  G4double p = sqrt (pow (x - centre.x (), 2) + pow (y - centre.y (), 2));
  if (p > 0.5) {
    red = green = blue = 0.;
  }
  else {
    G4double z = sqrt (pow (size,2) - pow (p, 2));
    G4Vector3D radiusVector (x, y, z);
    G4double brightness = radiusVector.unit () * lightpointDirection;
    if (brightness < 0.0) brightness = 0.0;
    G4Vector3D reflightDirection =
      2.0 * brightness * radiusVector.unit () - lightpointDirection;
    G4double reflectionDot =
      viewpointDirection.unit() * reflightDirection.unit();
    if (reflectionDot < 0) reflectionDot = 0;
    G4double specular = pow (reflectionDot, 10);
    red = ambient + (1.0 - ambient) * brightness;
    red += specular;
    green += specular;
    blue += specular;
  }

  **************************/

  ////////////////////////////////////////////////////////
  // Always check...
  ////////////////////////////////////////////////////////

  if (red   < 0.) red   = 0.;
  if (green < 0.) green = 0.;
  if (blue  < 0.) blue  = 0.;
  if (red   > 1.) red   = 1.;
  if (green > 1.) green = 1.;
  if (blue  > 1.) blue  = 1.;
}
