// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: locg4templates.hh,v 1.1 1999-01-08 16:31:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Include file used by applications
// to solve templates on some system.
//
// OSF1/cxx : don't forget to compile the file that contains this include 
//  with -define_templates
//
//    09/10/96 : creation. Guy Barrand
//    30/06/97 : HP_aCC directed template declarations. P.Mora de Freitas
//
#ifndef locg4templates_hh
#define locg4templates_hh

#ifdef G4_SOLVE_TEMPLATES

#include <rw/tpordvec.h>
#include "G4Allocator.hh"
#include "G4VSolid.hh"
static RWTPtrOrderedVector<G4VSolid>          dummy1;
#include "G4VPhysicalVolume.hh"
static RWTPtrOrderedVector<G4VPhysicalVolume> dummy2;
#include "G4LogicalVolume.hh"
static RWTPtrOrderedVector<G4LogicalVolume>   dummy3;
#include "G4SmartVoxelProxy.hh"
static RWTPtrOrderedVector<G4SmartVoxelProxy> dummy5;
#include "G4SmartVoxelNode.hh"
static RWTPtrOrderedVector<G4SmartVoxelNode>  dummy6;
#include "G4ThreeVector.hh"
static RWTValOrderedVector<G4ThreeVector>     dummy20;
static RWTValOrderedVector<G4int>             dummy21;
#include "G4OrderedTable.hh"
static RWTPtrOrderedVector<G4ValVector>       dummy24;
static RWTValVector<double>                   dummy26;
#include "G4Point3D.hh"
static RWTValOrderedVector<G4Point3D>            dummy41;
#include "G4VStateDependent.hh"
static RWTPtrOrderedVector<G4VStateDependent>    dummy48;
#include "G4AffineTransform.hh"
static RWTValVector<G4AffineTransform>           dummy49;
static RWTValVector<EVolume>                     dummy50;
#include "G4Transform3D.hh"
static RWTValVector<G4Transform3D>               dummy53;
static RWTValOrderedVector<G4Transform3D>        dummy54;


// To have a clean shared lib link with option -Wl,-v :
void dummyUse () {
 RWCString                dummy35;
 dummy24.clearAndDestroy  ();

// dummy26.boundsCheck     (0);
 double                   dummy40;
 dummy40                  = dummy26[0];


 RWBoolean                b;
}

//Solve visualization templates.
#ifdef G4_SOLVE_VIS_TEMPLATES

#include <G4RayView.hh>
static RWTValVector<G4RayView::ColourCell>                                     vis_1;
inline static unsigned vis_hash_fun_1 (const G4RayView::RayCoordinate&) {return (unsigned) 1;}
static RWTValHashDictionary<G4RayView::RayCoordinate, G4RayView::RayHitColour> vis_2(vis_hash_fun_1);
inline static unsigned vis_hash_fun_2 (const G4RayView::RayHitColour&)  {return (unsigned) 1;}
static RWTValHashDictionary<G4RayView::RayHitColour, int>                      vis_3(vis_hash_fun_2);

#endif //G4_SOLVE_VIS_TEMPLATES

#endif //G4_SOLVE_TEMPLATES

//
// Directed instantiations to solve templates
// on HP aCC compiler.
//
//
#ifdef HP_aCC
#include <rw/tpordvec.h>
// #include "G4PhysicsVector.hh"
// template class RWTPtrOrderedVector<G4PhysicsVector>;
#include "G4VSolid.hh"
template class RWTPtrOrderedVector<G4VSolid>;
#include "globals.hh"
#include <rw/tvhdict.h>
inline static unsigned dummyHashFun1 (const unsigned long&) {return (unsigned) 1;}
static RWTValHashDictionary<unsigned long,int> dummy_aCC_2047 (dummyHashFun1);
#include <CLHEP/Vector/ThreeVector.h>
template class RWTValOrderedVector<Hep3Vector>;
static RWTValOrderedVector<Hep3Vector> dummy_aCC_1011;
#include "G4VPhysicalVolume.hh"
template class RWTPtrOrderedVector<G4VPhysicalVolume>;
#include "G4LogicalVolume.hh"
template class RWTPtrOrderedVector<G4LogicalVolume>;
#include "G4SmartVoxelNode.hh"
template class RWTPtrOrderedVector<G4SmartVoxelNode>;
#include "G4SmartVoxelProxy.hh"
template class RWTPtrOrderedVector<G4SmartVoxelProxy>;
template class RWTValOrderedVector<G4int>;
#include "G4OrderedTable.hh"
template class RWTValOrderedVector<G4double>;
static RWTPtrOrderedVector<G4ValVector> dummy_aCC_1010;
template class RWTPtrOrderedVector<G4ValVector>;
template class RWTValVector<double>;
static RWTValVector<double> dummy_aCC_1014;
#include "G4Point3D.hh"
template class RWTValOrderedVector<G4Point3D>;
template class RWTValVector<G4Point3D>;
#include "G4Plane3D.hh"
template class RWTValOrderedVector<G4Plane3D>;
template class RWTValVector<G4Plane3D>;
#include "G4VStateDependent.hh"
template class RWTPtrOrderedVector<G4VStateDependent>;
#include "G4AffineTransform.hh"
template class RWTValVector<G4AffineTransform>;
static RWTValVector<G4AffineTransform> dummy_aCC_1015;
template class RWTValVector<EVolume>;
static RWTValVector<EVolume> dummy_aCC_1016;
#include "G4Transform3D.hh"
template class RWTValOrderedVector<G4Transform3D>;
template class RWTValOrderedVector<RWCString>;
template class RWTValVector<int>;
static RWTValVector<int> dummy_aCC_1017;
template class RWTValVector<EAxis>;
static RWTValVector<EAxis> dummy_aCC_1018;

#ifdef G4VISMANAGER_HH
#include "G4VScene.hh"
template class RWTPtrOrderedVector<G4VScene>;
#include "G4VView.hh"
template class RWTPtrOrderedVector<G4VView>;
#include "G4VGraphicsSystem.hh"
template class RWTPtrOrderedVector<G4VGraphicsSystem>;
//#include "G4RayView.hh"
//template class RWTValHashDictionary<G4RayView::RayHitColour,int>;
//template class RWTValHashDictionary<G4RayView::RayCoordinate,G4RayView::RayHitColour>;
//template class RWTValOrderedVector<G4RayView::ColourCell>;
#endif // G4VISMANAGER_HH
#endif //HP_aCC

#endif //locg4templates_hh

