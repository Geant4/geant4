#include "ExP02GeoTree.hh"
//
#include "ExP02DetectorConstruction.hh"
//
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Transform3D.h"
//
#include "G4RunManager.hh"
#include "G4VisExtent.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VSolid.hh"
#include "G4VCSGfaceted.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Polyhedron.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyPhiFace.hh"
#include "G4IntersectingCone.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Cons.hh"
#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPVParameterisation.hh"
#include "G4VUserRegionInformation.hh"
#include "G4UserLimits.hh"
//
std::vector<CLHEP::Hep3Vector> a;
std::vector<G4VCSGface*> b;
std::vector<G4PolyconeSideRZ> c;
std::vector<G4PolyPhiFaceEdge> d;
std::vector<G4PolyPhiFaceVertex> e;
//
#undef __G4String
