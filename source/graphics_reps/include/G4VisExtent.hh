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
// $Id$
//
// 
// A.Walkden 28/11/95

// Class Description:
// G4VisExtent defines a bounding box in a visualisable object's local
// coordinate system which includes the object.
// WARNING: it also attempts to support the concept of a bounding
// sphere.  (This is used extensively in the G4 Visualisation System
// to calculate camera parameters, etc.)  Be aware that this involves
// loss of information.  Given a bounding box, one can calculate the
// bounding sphere; inverting this will produce a cube which *might*
// *not* include the object.  E.g., a long thin object of length l
// will have a bounding sphere of diameter l; the corresponding cube
// will have side l/std::sqrt(3) so that the bounding sphere diameter is
// still l.  Thus the long thin object will stick out of the cube.
// So, if you once use the concept of bounding sphere you must stick
// with it and abandon the concept of bounding box.
// Class Description - End:

#ifndef G4VISEXTENT_HH
#define G4VISEXTENT_HH

#include "globals.hh"
#include "G4Point3D.hh"

class G4VisExtent
{
public: // With description

  G4VisExtent (G4double xmin = 0., G4double xmax = 0., 
               G4double ymin = 0., G4double ymax = 0., 
               G4double zmin = 0., G4double zmax = 0.);
  G4VisExtent (const G4Point3D& centre, G4double radius);
  ~G4VisExtent ();
  G4bool operator != (const G4VisExtent& e) const;
  G4double  GetXmin         () const;
  G4double  GetXmax         () const;
  G4double  GetYmin         () const;
  G4double  GetYmax         () const;
  G4double  GetZmin         () const;
  G4double  GetZmax         () const;
  const G4Point3D& GetExtentCentre () const;
  const G4Point3D& GetExtentCenter () const;
  G4double  GetExtentRadius () const;
  void SetXmin (G4double xmin);
  void SetXmax (G4double xmax);
  void SetYmin (G4double ymin);
  void SetYmax (G4double ymax);
  void SetZmin (G4double zmin);
  void SetZmax (G4double zmax);
  
  void SetExtent(size_t axis/*0 --> 2*/, G4bool min_max, G4double value);
  G4double GetExtent(size_t axis/*0 --> 2*/, G4bool min_max);
  
  friend std::ostream& operator << (std::ostream& os, const G4VisExtent& e);
  static const G4VisExtent NullExtent;

private:
  G4double fExtent[3][2];  
  mutable G4bool fRadiusCached, fCentreCached;
  mutable G4double fRadius;
  mutable G4Point3D fCentre;
};

inline G4double G4VisExtent::GetXmin () const { return fExtent[0][0]; }
inline G4double G4VisExtent::GetXmax () const { return fExtent[0][1]; }
inline G4double G4VisExtent::GetYmin () const { return fExtent[1][0]; }
inline G4double G4VisExtent::GetYmax () const { return fExtent[1][1]; }
inline G4double G4VisExtent::GetZmin () const { return fExtent[2][0]; }
inline G4double G4VisExtent::GetZmax () const { return fExtent[2][1]; }

inline const G4Point3D& G4VisExtent::GetExtentCenter () const {
  return GetExtentCentre ();
}

inline void G4VisExtent::SetXmin (G4double xmin)
{fExtent[0][0] = xmin; fRadiusCached = false; fCentreCached = false;}
inline void G4VisExtent::SetXmax (G4double xmax)
{fExtent[0][1] = xmax; fRadiusCached = false; fCentreCached = false;}
inline void G4VisExtent::SetYmin (G4double ymin)
{fExtent[1][0] = ymin; fRadiusCached = false; fCentreCached = false;}
inline void G4VisExtent::SetYmax (G4double ymax)
{fExtent[1][1] = ymax; fRadiusCached = false; fCentreCached = false;}
inline void G4VisExtent::SetZmin (G4double zmin)
{fExtent[2][0] = zmin; fRadiusCached = false; fCentreCached = false;}
inline void G4VisExtent::SetZmax (G4double zmax)
{fExtent[2][1] = zmax; fRadiusCached = false; fCentreCached = false;}

inline void G4VisExtent::SetExtent(size_t axis, G4bool min_max, G4double value)
{
    fExtent[axis][min_max] = value;
}    
inline G4double G4VisExtent::GetExtent(size_t axis, G4bool min_max)
{
    return fExtent[axis][min_max];
}
#endif
