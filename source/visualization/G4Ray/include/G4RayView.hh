// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayView.hh,v 2.2 1998/11/06 13:41:51 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing view.

#include "G4RunManager.hh"
#include "G4RaySteppingAction.hh"
#include "G4RayPrimaryGeneratorAction.hh"


#ifndef G4RAYVIEW_HH
#define G4RAYVIEW_HH

#include "G4VView.hh"
#include "G4Colour.hh"

#include "G4VPhysicalVolume.hh"

#include <rw/tvordvec.h>
#include <rw/tvhdict.h>

class G4RayScene;
class G4RaySteppingAction;

class G4RayColourCell {
public:
  G4RayColourCell (): red (0), green (0), blue (0), nHits (0) {}
  G4RayColourCell (float rr, float gg, float bb, G4int ii):
    red (rr), green (gg), blue (bb), nHits (ii) {}
  float red, green, blue;
  G4int nHits;  // No. of ray hits using this colour.
  G4bool operator == (const G4RayColourCell& h) const {
    return this == &h;
  }
};

class G4RayCoordinate {
public:
  G4RayCoordinate (): coord (0) {}
  G4RayCoordinate (G4double x, G4double y):
    coord (((unsigned int) (0x7fff * (x + 1.0))) | 
	     ((unsigned int) (0x7fff * (1.0 - y)) << 16)) {}
  // Assumes -1.0 <= x, y <= 1.0.  Assigns 16 bits each.
  unsigned int coord;
  G4bool operator == (const G4RayCoordinate& c) const {
    return coord == c.coord;
  }
};

class G4RayHitColour {
public:
  G4RayHitColour (): colour (0) {}
  G4RayHitColour (float r, float g, float b):
    colour ((((unsigned int) (0x3ff * r)) << 20) |
	      (((unsigned int) (0x3ff * g)) << 10) |
	      ((unsigned int) (0x3ff * b))) {}
  // Assumes 0.0 <= r, g, b <= 1.0.  Assigns 10 bits each.
  unsigned int colour;
  G4bool operator == (const G4RayHitColour& h) const {
    return colour == h.colour;
  }
  float GetRed   () {return (float ((colour >> 20)        )) / 0x3ff;}
  float GetGreen () {return (float ((colour >> 10) & 0x3ff)) / 0x3ff;}
  float GetBlue  () {return (float ((colour      ) & 0x3ff)) / 0x3ff;}
};

class G4RayView: virtual public G4VView {
  
public:

  G4RayView (G4RayScene& scene, const G4String& name = "");

  ~G4RayView();

  G4double GetResolution () const;

  void SetResolution (G4double resolution);

protected:

  void SetView ();

  void InitialiseImage ();

  void TriggerRedraw ();

  G4bool Scan (int& i, int& j, int& k,
	       float& red, float& green, float& blue,
	       int& iPixVal);
  // Iterate over window, return false when finished.  Return i, j,
  // and k which define a coloured Square of side k whose upper left
  // corner is (i, j).  Note origin is at upper left corner of window
  // - X and PostScript convention.  iPixVal is the pixel value for
  // this colour, i.e., the key in fRayHitColourIndex.

  // The following static members define an image.  They are static so
  // that (a) there is only ever one of them (they can take a lot of
  // memory) and (b) any object of the class or its descendants can
  // access it.  It is *not* cleared by ClearView ().  It is only
  // cleared by a change of view parameters and some other crucial
  // parameters.  It can be progressively refined by multiple calls to
  // DrawView.
  static G4ViewParameters fImageViewParameters;
  static G4double         fInitialResolution;
  static G4double         fRequestedResolution;
  static G4double         fCurrentResolution;
  static G4float          fMaxMinColourSeparation;
  static G4int            fImageWidth, fImageHeight;
  static G4int            fMaxImageMapEntries;
  static G4bool           fImageMapRejigged;
  static G4bool           fTriggerRescan;
  static G4int            fWaitLimit;  // Waits between refresh.
  static G4int            fWaitCount;
  static RWTValOrderedVector <G4RayColourCell> fRayHitColourMap;
  static RWTValHashDictionary <G4RayCoordinate, G4RayHitColour> fRayHits;
  static unsigned int     RayCoordHash (const G4RayCoordinate& rayCoord);
  static RWTValHashDictionary <G4RayHitColour, G4int> fRayHitColourIndex;
  static unsigned int     RayColourHash (const G4RayHitColour& rayColour);

private:
  void Trace (G4double x, G4double y, float& red, float& green, float& blue);
  // Trace the ray from (x, y) in screen coordinates (-1 <= x,y <= 1)
  // and fill colour.

  void CompactColour (G4bool Print);

  float ColourSeparation (float ri, float gi, float bi,
			  float rj, float gj, float bj);

  static void ChangePixVals (const G4RayHitColour& colour,
			     G4int& index, void* ij);

  static void DecrementPixVals (const G4RayHitColour& colour,
				G4int& index, void* ij);

  static void PrintPixVals (const G4RayHitColour& colour,
			    G4int& index, void* ij);

  G4RayScene&                 fScene; // Graphics Scene for this view.
  
  G4Vector3D                  fXprime;
  G4Vector3D                  fYprime;

  G4Point3D                   fCenterofscreen;
  G4Point3D                   fCamerapoint;
  
  G4RunManager*               fG4RayRunManager ;
  G4RaySteppingAction*        fG4RaySteppingAction;
  G4RayPrimaryGeneratorAction*   fG4RayGenerator;
 
  G4ParticleGun *             fG4RayGun ;

  G4VPhysicalVolume * fpWorld;

};

#include "G4RayView.icc"

#endif
