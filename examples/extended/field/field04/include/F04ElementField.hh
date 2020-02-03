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
/// \file field/field04/include/F04ElementField.hh
/// \brief Definition of the F04ElementField class
//

#ifndef F04ElementField_h
#define F04ElementField_h 1

#include "globals.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

//  class F04ElementField - interface for the EM field of one element

//  This is the interface class used by GlobalField to compute the field
//  value at a given point[].

//  An element that represents an element with an EM field will
//  derive a class from this one and implement the computation for the
//  element. The Construct() function will add the derived object into
//  GlobalField.

class F04ElementField 
{

  private:

    F04ElementField& operator=(const F04ElementField&);

  public:

    ///  Constructor.
    F04ElementField(const G4ThreeVector, G4LogicalVolume*);

    /// the actual implementation constructs the F04ElementField
    void Construct();

    ///  Destructor.
    virtual ~F04ElementField() {}

    /// SetMaxStep(G4double) sets the max. step size
    void SetMaxStep(G4double stp)
    {
      fMaxStep = stp;
      fUserLimits->SetMaxAllowedStep(fMaxStep);
      fVolume->SetUserLimits(fUserLimits);
    }

    /// GetMaxStep() returns the max. step size
    G4double GetMaxStep() { return fMaxStep; }

    /// SetColor(G4String) sets the color
    void SetColor(G4String c)
    {
      fColor = c;
      fVolume->SetVisAttributes(GetVisAttribute(fColor));
    }

    /// GetColor() returns the color
    G4String GetColor() { return fColor; }

    ///  GetVisAttribute() returns the appropriate G4VisAttributes.
    static G4VisAttributes* GetVisAttribute(G4String color);

    ///  SetGlobalPoint() ensures that the point is within the global
    ///  bounding box of this ElementField's global coordinates.
    ///  Normally called 8 times for the corners of the local bounding
    ///  box, after a local->global coordinate transform.
    ///  If never called, the global bounding box is infinite.
    ///  BEWARE: if called only once, the bounding box is just a point.
    void SetGlobalPoint(const G4double point[4])
    {
      if(fMinX == -DBL_MAX || fMinX > point[0]) fMinX = point[0];
      if(fMinY == -DBL_MAX || fMinY > point[1]) fMinY = point[1];
      if(fMinZ == -DBL_MAX || fMinZ > point[2]) fMinZ = point[2];
      if(fMaxX ==  DBL_MAX || fMaxX < point[0]) fMaxX = point[0];
      if(fMaxY ==  DBL_MAX || fMaxY < point[1]) fMaxY = point[1];
      if(fMaxZ ==  DBL_MAX || fMaxZ < point[2]) fMaxZ = point[2];
    }

    ///  IsInBoundingBox() returns true if the point is within the
    ///  global bounding box - global coordinates.
    bool IsInBoundingBox(const G4double point[4]) const
    {
      if(point[2] < fMinZ || point[2] > fMaxZ) return false;
      if(point[0] < fMinX || point[0] > fMaxX) return false;
      if(point[1] < fMinY || point[1] > fMaxY) return false;
      return true;
    }

    ///  AddFieldValue() will add the field value for this element to field[].
    ///  Implementations must be sure to verify that point[] is within
    ///  the field region, and do nothing if not.
    ///  point[] is in global coordinates and geant4 units; x,y,z,t.
    ///  field[] is in geant4 units; Bx,By,Bz,Ex,Ey,Ez.
    ///  For efficiency, the caller may (but need not) call
    ///  IsInBoundingBox(point), and only call this function if that
    ///  returns true.
    virtual void
        AddFieldValue(const G4double point[4], G4double field[6]) const = 0;

    virtual G4double GetLength() = 0;
    virtual G4double GetWidth()  = 0;
    virtual G4double GetHeight() = 0;

  protected:

    G4LogicalVolume* fVolume;

    G4AffineTransform fGlobal2local;

//    F04ElementField(const F04ElementField&);

  private:

    static G4ThreadLocal G4Navigator* fNavigator;

    G4String fColor;

    G4ThreeVector fCenter;
    G4double fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ;

    G4double fMaxStep;
    G4UserLimits* fUserLimits;

};

#endif
