//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisExtent.hh,v 1.6 2001-07-24 21:39:43 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// will have side l/sqrt(3) so that the bounding sphere diameter is
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
  G4Point3D GetExtentCentre () const;
  G4Point3D GetExtentCenter () const;
  G4double  GetExtentRadius () const;
  void SetXmin (G4double xmin);
  void SetXmax (G4double xmax);
  void SetYmin (G4double ymin);
  void SetYmax (G4double ymax);
  void SetZmin (G4double zmin);
  void SetZmax (G4double zmax);
  friend G4std::ostream& operator << (G4std::ostream& os, const G4VisExtent& e);

private:
  G4double fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;
};

inline G4double G4VisExtent::GetXmin () const { return fXmin; }
inline G4double G4VisExtent::GetXmax () const { return fXmax; }
inline G4double G4VisExtent::GetYmin () const { return fYmin; }
inline G4double G4VisExtent::GetYmax () const { return fYmax; }
inline G4double G4VisExtent::GetZmin () const { return fZmin; }
inline G4double G4VisExtent::GetZmax () const { return fZmax; }

inline G4Point3D G4VisExtent::GetExtentCenter () const {
  return GetExtentCentre ();
}

inline void G4VisExtent::SetXmin (G4double xmin) {fXmin = xmin;}
inline void G4VisExtent::SetXmax (G4double xmax) {fXmax = xmax;}
inline void G4VisExtent::SetYmin (G4double ymin) {fYmin = ymin;}
inline void G4VisExtent::SetYmax (G4double ymax) {fYmax = ymax;}
inline void G4VisExtent::SetZmin (G4double zmin) {fZmin = zmin;}
inline void G4VisExtent::SetZmax (G4double zmax) {fZmax = zmax;}

#endif
