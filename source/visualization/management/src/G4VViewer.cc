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
//
// 
// John Allison  27th March 1996
// Abstract interface class for graphics views.

#include "G4VViewer.hh"

#include "G4Timer.hh"

#include "G4ios.hh"
#include <sstream>

#include "G4VisManager.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4UImanager.hh"

G4VViewer::G4VViewer (G4VSceneHandler& sceneHandler,
		      G4int id, const G4String& name):
fSceneHandler (sceneHandler),
fViewId (id),
//fModified (true),
fNeedKernelVisit (true)
{
  if (name == "") {
    std::ostringstream ost;
    ost << fSceneHandler.GetName () << '-' << fViewId;
    fName = ost.str();
  }
  else {
    fName = name;
  }
  fShortName = fName.substr(0, fName.find (' '));
  G4StrUtil::strip(fShortName);

  fVP = G4VisManager::GetInstance()->GetDefaultViewParameters();
  fDefaultVP = fVP;
}

G4VViewer::~G4VViewer () {
  fSceneHandler.RemoveViewerFromList(this);
}

void G4VViewer::SetName (const G4String& name) {
  fName = name;
  fShortName = fName.substr(0, fName.find (' '));
  G4StrUtil::strip(fShortName);
}

void G4VViewer::NeedKernelVisit () {

  fNeedKernelVisit = true;

  // At one time I thought we'd better notify all viewers.  But I guess
  // each viewer can take care of itself, so the following code is
  // redundant (but keep it commented out for now).   (John Allison)
  // Notify all viewers that a kernel visit is required.
  // const G4ViewerList& viewerList = fSceneHandler.GetViewerList ();
  // G4ViewerListConstIterator i;
  // for (i = viewerList.begin(); i != viewerList.end(); i++) {
  //   (*i) -> SetNeedKernelVisit ();
  // }
  // ??...but, there's a problem in OpenGL Stored which seems to
  // require *all* viewers to revisit the kernel, so...
  //  const G4ViewerList& viewerList = fSceneHandler.GetViewerList ();
  //  G4ViewerListConstIterator i;
  //  for (i = viewerList.begin(); i != viewerList.end(); i++) {
  //    (*i) -> SetNeedKernelVisit (true);
  //  }
  // Feb 2005 - commented out.  Let's fix OpenGL if necessary.
}

void G4VViewer::FinishView () {}

void G4VViewer::ShowView () {}

void G4VViewer::ProcessView ()
{
  // If the scene has changed, or if the concrete viewer has decided
  // that it necessary to visit the kernel, perhaps because the view
  // parameters have changed significantly (this should be done in the
  // concrete viewer's DrawView)...
  if (fNeedKernelVisit) {
    // Reset flag.  This must be done before ProcessScene to prevent
    // recursive calls when recomputing transients...
    G4Timer timer;
    timer.Start();
    fNeedKernelVisit = false;
    fSceneHandler.ClearStore ();
    fSceneHandler.ProcessScene ();
    timer.Stop();
    fKernelVisitElapsedTimeSeconds = timer.GetRealElapsed();
  }
}

void G4VViewer::SetViewParameters (const G4ViewParameters& vp) {
  fVP = vp;
}

void G4VViewer::SetTouchable
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath)
{
  // Set the touchable for /vis/touchable/set/... commands.
  std::ostringstream oss;
  const auto& pvStore = G4PhysicalVolumeStore::GetInstance();
  for (const auto& pvNodeId: fullPath) {
    const auto& pv = pvNodeId.GetPhysicalVolume();
    auto iterator = find(pvStore->cbegin(),pvStore->cend(),pv);
    if (iterator == pvStore->cend()) {
      G4ExceptionDescription ed;
      ed << "Volume no longer in physical volume store.";
      G4Exception("G4VViewer::SetTouchable", "visman0501", JustWarning, ed);
    } else {
      oss
      << ' ' << pvNodeId.GetPhysicalVolume()->GetName()
      << ' ' << pvNodeId.GetCopyNo();
    }
  }
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/set/touchable" + oss.str());
}

void G4VViewer::TouchableSetVisibility
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
 G4bool visibiity)
{
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  // The following is equivalent to
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/touchable/set/visibility ...");
  // (assuming the touchable has already been set), but avoids view rebuild.

  // Instantiate a working copy of a G4VisAttributes object...
  G4VisAttributes workingVisAtts;
  // and set the visibility.
  workingVisAtts.SetVisibility(visibiity);

  fVP.AddVisAttributesModifier
  (G4ModelingParameters::VisAttributesModifier
   (workingVisAtts,
    G4ModelingParameters::VASVisibility,
    G4PhysicalVolumeModel::GetPVNameCopyNoPath(fullPath)));
  // G4ModelingParameters::VASVisibility (VAS = Vis Attribute Signifier)
  // signifies that it is the visibility that should be picked out
  // and merged with the touchable's normal vis attributes.
}

void G4VViewer::TouchableSetColour
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
 const G4Colour& colour)
{
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  // The following is equivalent to
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/touchable/set/colour ...");
  // (assuming the touchable has already been set), but avoids view rebuild.

  // Instantiate a working copy of a G4VisAttributes object...
  G4VisAttributes workingVisAtts;
  // and set the colour.
  workingVisAtts.SetColour(colour);

  fVP.AddVisAttributesModifier
  (G4ModelingParameters::VisAttributesModifier
   (workingVisAtts,
    G4ModelingParameters::VASColour,
    G4PhysicalVolumeModel::GetPVNameCopyNoPath(fullPath)));
  // G4ModelingParameters::VASColour (VAS = Vis Attribute Signifier)
  // signifies that it is the colour that should be picked out
  // and merged with the touchable's normal vis attributes.
}

std::vector <G4ThreeVector> G4VViewer::ComputeFlyThrough(G4Vector3D* /*aVect*/)
{
    enum CurveType {
        Bezier,
        G4SplineTest};
    
    // Choose a curve type (for testing)
//    int myCurveType = Bezier;

    // number if step points
    G4int stepPoints = 500;

    
    G4Spline spline;

    
    // At the moment we don't use the aVect parameters, but build it here :
    // Good step points for exampleB5
    spline.AddSplinePoint(G4Vector3D(0,1000,-14000));
    spline.AddSplinePoint(G4Vector3D(0,1000,0));
    spline.AddSplinePoint(G4Vector3D(-4000,1000,4000));

    
    std::vector <G4ThreeVector> viewVect;

//    if(myCurveType == Bezier) {

        
        // Draw the spline
        
        for (G4int i = 0; i < stepPoints; ++i) {
            G4float t = (G4float)i / (G4float)stepPoints;
            G4Vector3D cameraPosition = spline.GetInterpolatedSplinePoint(t);
            //        G4Vector3D targetPoint = spline.GetInterpolatedSplinePoint(t);
            
            //        viewParam->SetViewAndLights(G4ThreeVector (cameraPosition.x(), cameraPosition.y(), cameraPosition.z()));
            //        viewParam->SetCurrentTargetPoint(targetPoint);
            G4cout << "FLY CR("<< i << "):" << cameraPosition << G4endl;
            viewVect.push_back(G4ThreeVector (cameraPosition.x(), cameraPosition.y(), cameraPosition.z()));
        }
        
//    } else if (myCurveType == G4SplineTest) {
        /*
         This method is a inspire from a Bezier curve. The problem of the Bezier curve is that the path does not go straight between two waypoints.
         This method add "stay straight" parameter which could be between 0 and 1 where the pass will follow exactly the line between the waypoints
         Ex : stay straight = 50%
         m1 = 3*(P1+P0)/2
         
         Ex : stay straight = 0%
         m1 = (P1+P0)/2
         
         P1
         /  \
         /    \
         a--x--b
         /  °  °  \
         / °      ° \
         m1           m2
         /              \
         /                \
         /                  \
         /                    \
         P0                     P2
         
         */
//        G4Vector3D a;
//        G4Vector3D b;
//        G4Vector3D m1;
//        G4Vector3D m2;
//        G4Vector3D P0;
//        G4Vector3D P1;
//        G4Vector3D P2;
//        G4double stayStraight = 0;
//        G4double bezierSpeed = 0.4; // Spend 40% time in bezier curve (time between m1-m2 is 40% of time between P0-P1)
//        
//        G4Vector3D firstPoint;
//        G4Vector3D lastPoint;
//        
//        float nbBezierSteps = (stepPoints * bezierSpeed*(1-stayStraight)) * (2./spline.GetNumPoints());
//        float nbFirstSteps = ((stepPoints/2-nbBezierSteps/2) /(1+stayStraight)) * (2./spline.GetNumPoints());
//        
//        // First points
//        firstPoint = spline.GetPoint(0);
//        lastPoint = (firstPoint + spline.GetPoint(1))/2;
//        
//        for( float j=0; j<1; j+= 1/nbFirstSteps) {
//            G4ThreeVector pt = firstPoint + (lastPoint - firstPoint) * j;
//            viewVect.push_back(pt);
//            G4cout << "FLY Bezier A1("<< viewVect.size()<< "):" << pt << G4endl;
//        }
//        
//        for (int i = 0; i < spline.GetNumPoints()-2; i++) {
//            P0 = spline.GetPoint(i);
//            P1 = spline.GetPoint(i+1);
//            P2 = spline.GetPoint(i+2);
//            
//            m1 = P1 - (P1-P0)*(1-stayStraight)/2;
//            m2 = P1 + (P2-P1)*(1-stayStraight)/2;
//            
//            // We have to get straight path from (middile of P0-P1) to (middile of P0-P1 + (dist P0-P1) * stayStraight/2)
//            if (stayStraight >0) {
//                
//                firstPoint = (P0 + P1)/2;
//                lastPoint = (P0 + P1)/2 + (P1-P0)*stayStraight/2;
//                
//                for( float j=0; j<1; j+= 1/(nbFirstSteps*stayStraight)) {
//                    G4ThreeVector pt = firstPoint + (lastPoint - firstPoint)* j;
//                    viewVect.push_back(pt);
//                    G4cout << "FLY Bezier A2("<< viewVect.size()<< "):" << pt << G4endl;
//                }
//            }
//            // Compute Bezier curve
//            for( float delta = 0 ; delta < 1 ; delta += 1/nbBezierSteps)
//            {
//                // The Green Line
//                a = m1 + ( (P1 - m1) * delta );
//                b = P1 + ( (m2 - P1) * delta );
//                
//                // Final point
//                G4ThreeVector pt = a + ((b-a) * delta );
//                viewVect.push_back(pt);
//                G4cout << "FLY Bezier("<< viewVect.size()<< "):" << pt << G4endl;
//            }
//            
//            // We have to get straight path
//            if (stayStraight >0) {
//                firstPoint = (P1 + P2)/2 - (P2-P1)*stayStraight/2;
//                lastPoint = (P1 + P2)/2;
//                
//                for( float j=0; j<1; j+= 1/(nbFirstSteps*stayStraight)) {
//                    G4ThreeVector pt = firstPoint + (lastPoint - firstPoint)* j;
//                    viewVect.push_back(pt);
//                    G4cout << "FLY Bezier B1("<< viewVect.size()<< "):" << pt << G4endl;
//                }
//            }
//        }
//        
//        // last points
//        firstPoint = spline.GetPoint(spline.GetNumPoints()-2);
//        lastPoint = spline.GetPoint(spline.GetNumPoints()-1);
//        for( float j=1; j>0; j-= 1/nbFirstSteps) {
//            G4ThreeVector pt = lastPoint - ((lastPoint-firstPoint)*((1-stayStraight)/2) * j );
//            viewVect.push_back(pt);
//            G4cout << "FLY Bezier B2("<< viewVect.size()<< "):" << pt << G4endl;
//        }
//    }
    return viewVect;
}


#ifdef G4MULTITHREADED

void G4VViewer::DoneWithMasterThread () {
  // G4cout << "G4VViewer::DoneWithMasterThread" << G4endl;
}

void G4VViewer::MovingToMasterThread () {
  // G4cout << "G4VViewer::MovingToMasterThread" << G4endl;
}

void G4VViewer::SwitchToVisSubThread () {
  // G4cout << "G4VViewer::SwitchToVisSubThread" << G4endl;
}

void G4VViewer::DoneWithVisSubThread () {
  // G4cout << "G4VViewer::DoneWithVisSubThread" << G4endl;
}

void G4VViewer::MovingToVisSubThread () {
  // G4cout << "G4VViewer::MovingToVisSubThread" << G4endl;
}

void G4VViewer::SwitchToMasterThread () {
  // G4cout << "G4VViewer::SwitchToMasterThread" << G4endl;
}

#endif

std::ostream& operator << (std::ostream& os, const G4VViewer& v) {
  os << "View " << v.fName << ":\n";
  os << v.fVP;
  return os;
}


// ===== G4Spline class =====

G4VViewer::G4Spline::G4Spline()
: vp(), delta_t(0)
{
}


G4VViewer::G4Spline::~G4Spline()
{}

// Solve the Catmull-Rom parametric equation for a given time(t) and vector quadruple (p1,p2,p3,p4)
G4Vector3D G4VViewer::G4Spline::CatmullRom_Eq(G4float t, const G4Vector3D& p1, const G4Vector3D& p2, const G4Vector3D& p3, const G4Vector3D& p4)
{
    G4float t2 = t * t;
    G4float t3 = t2 * t;
    
    G4float b1 = .5 * (  -t3 + 2*t2 - t);
    G4float b2 = .5 * ( 3*t3 - 5*t2 + 2);
    G4float b3 = .5 * (-3*t3 + 4*t2 + t);
    G4float b4 = .5 * (   t3 -   t2    );
    
    return (p1*b1 + p2*b2 + p3*b3 + p4*b4);
}

void G4VViewer::G4Spline::AddSplinePoint(const G4Vector3D& v)
{
    vp.push_back(v);
    delta_t = (G4float)1 / (G4float)vp.size();
}


G4Vector3D G4VViewer::G4Spline::GetPoint(G4int a)
{
    return vp[a];
}

G4int G4VViewer::G4Spline::GetNumPoints()
{
    return (G4int)vp.size();
}

G4Vector3D G4VViewer::G4Spline::GetInterpolatedSplinePoint(G4float t)
{
    // Find out in which interval we are on the spline
    G4int p = (G4int)(t / delta_t);
    // Compute local control point indices
#define BOUNDS(pp) { if (pp < 0) pp = 0; else if (pp >= (G4int)vp.size()-1) pp = (G4int)vp.size() - 1; }
    G4int p0 = p - 1;     BOUNDS(p0);
    G4int p1 = p;         BOUNDS(p1);
    G4int p2 = p + 1;     BOUNDS(p2);
    G4int p3 = p + 2;     BOUNDS(p3);
    // Relative (local) time
    G4float lt = (t - delta_t*p) / delta_t;
    // Interpolate
    return CatmullRom_Eq(lt, vp[p0], vp[p1], vp[p2], vp[p3]);
}
