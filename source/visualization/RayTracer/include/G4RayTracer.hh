
#ifndef G4RayTracer_H
#define G4RayTracer_H 1
// 
//
// G4RayTracer
//

#if defined (G4VIS_BUILD_RAYTRACER_DRIVER) || defined (G4VIS_USE_RAYTRACER)

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VGraphicsSystem.hh"

class G4Event;
class G4EventManager;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4RTTrackingAction;
class G4RTSteppingAction;
class G4Colour;
class G4RTMessenger;
class G4RayShooter;
class G4VFigureFileMaker;
class G4RayTrajectoryPoint;
class G4VisAttributes;


class G4RayTracer : public G4VGraphicsSystem
{
  public:
    G4RayTracer(G4VFigureFileMaker* figMaker = 0);
    ~G4RayTracer();

  public:
    void Trace(G4String fileName);

  public:
    virtual G4VSceneHandler* CreateSceneHandler (const G4String& ) {;}
    virtual G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& ) {;}

  private:
    G4bool CreateBitMap();
    // Event loop
    void CreateFigureFile(G4String fileName);
    // Create figure file after an event loop
    G4bool GenerateColour(G4Event* anEvent);
    // Calcurate RGB for one trajectory
    void StoreUserActions();
    void RestoreUserActions();
    // Store and restore user action classes if defined

    G4Colour GetSurfaceColour(G4RayTrajectoryPoint* point);
    G4Colour GetMixedColour(G4Colour surfCol,G4Colour transCol,G4double weight=0.5);
    G4Colour Attenuate(G4RayTrajectoryPoint* point, G4Colour sourceCol);
    G4bool ValidColour(const G4VisAttributes* visAtt);

  public:
    inline void SetFigureFileMaker(G4VFigureFileMaker* figMaker)
    { theFigMaker = figMaker; }

  private:
    G4RayShooter * theRayShooter;
    G4VFigureFileMaker * theFigMaker;
    G4RTMessenger * theMessenger;

    G4EventManager * theEventManager;

    G4UserEventAction * theUserEventAction;
    G4UserStackingAction * theUserStackingAction;
    G4UserTrackingAction * theUserTrackingAction;
    G4UserSteppingAction * theUserSteppingAction;

    G4UserEventAction * theRayTracerEventAction;
    G4UserStackingAction * theRayTracerStackingAction;
    G4RTTrackingAction * theRayTracerTrackingAction;
    G4RTSteppingAction * theRayTracerSteppingAction;

    unsigned char* colorR;
    unsigned char* colorG;
    unsigned char* colorB;

    G4int nColumn;
    G4int nRow;

    G4ThreeVector eyePosition;
    G4ThreeVector targetPosition;
    G4ThreeVector eyeDirection;
    G4ThreeVector lightDirection;
    G4double headAngle;
    G4double viewSpan; // Angle per 100 pixels
    G4double attenuationLength;

    G4bool distortionOn;

  public:
    inline void SetNColumn(G4int val) { nColumn = val; }
    inline G4int GetNColumn() const { return nColumn; }
    inline void SetNRow(G4int val) { nRow = val; }
    inline G4int GetNRow() const { return nRow; }
    inline void SetEyePosition(G4ThreeVector val) { eyePosition = val; }
    inline G4ThreeVector GetEyePosition() const { return eyePosition; }
    inline void SetTargetPosition(G4ThreeVector val) { targetPosition = val; }
    inline G4ThreeVector GetTargetPosition() const { return targetPosition; }
    inline void SetLightDirection(G4ThreeVector val) { lightDirection = val.unit(); }
    inline G4ThreeVector GetLightDirection() const { return lightDirection; }
    inline void SetHeadAngle(G4double val) { headAngle = val; }
    inline G4double GetHeadAngle() const { return headAngle; }
    inline void SetViewSpan(G4double val) { viewSpan = val; }
    inline G4double GetViewSpan() const { return viewSpan; }
    inline void SetAttenuationLength(G4double val) { attenuationLength = val; }
    inline G4double GetAttenuationLength() const { return attenuationLength; }
    inline void SetDistortion(G4bool val) { distortionOn = val; }
    inline G4bool GetDistortion() const { return distortionOn; }
};

#endif
#endif


