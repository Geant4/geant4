#ifndef G4GEOMETRYCREATOR_HH
#define G4GEOMETRYCREATOR_HH
//#include "G4OrderedTable.hh"
#include "G4STEPEntity.hh"
#include "G4NISTStepReader.hh"
#include "G4Curve.hh"
#include "G4PlacedSolid.hh"
//#include "G4Placement.hh"
#include "G4Axis2Placement3D.hh"

//typedef RWTPtrOrderedVector<G4Curve> G4CurveVector;
//typedef RWTPtrOrderedVector<G4CurveVector> G4BoundaryVector;
typedef RWTPtrOrderedVector<G4PlacedSolid> G4PlacedSolidVector;
//typedef RWTPtrOrderedVector<G4GeometryCreator> G4CreatorVector;
typedef RWTPtrOrderedVector<G4Surface> G4SurfaceVector;
typedef RWTPtrOrderedVector<G4BREPSolid> G4SolidVector;
//typedef RWTPtrOrderedVector<G4Placement> G4PlacementVector;

//#include "G4STEPEntity.hh"
//erator ==(class RWTPtrOrderedVector<G4Curve> a, class RWTPtrOrderedVector<G4Curve> b){}
#include "globals.hh"
#include "instmgr.h"
#include "STEPentity.h"
#include "STEPaggregate.h"

#include "G4Curve.hh"
//#include "G4Surface.hh"
#include "G4FPlane.hh"
#include "G4FConicalSurface.hh"
#include "G4FCylindricalSurface.hh"
#include "G4CylindricalSurface.hh"
#include "G4ToroidalSurface.hh"
#include "G4SphericalSurface.hh"
#include "G4PlacedSolid.hh"

class G4GeometryCreator
{
  
public:
  G4GeometryCreator(){}
  virtual ~G4GeometryCreator(){};

  virtual void CreateG4Geometry(STEPentity&)=0;
  virtual void CreateSTEPGeometry(void* =0)=0;

  virtual void* GetCreatedObject(){return createdObject;}
  
  virtual G4String Name()=0;
  virtual G4bool operator==(const G4GeometryCreator&){return 0;}
  virtual STEPattribute* GetNamedAttribute(G4String&,STEPentity&);
  virtual STEPentity* GetNamedEntity(G4String&,STEPentity&);  
  G4int GetNextId(){objectId +=10 ; return objectId; }
  
  static G4int objectId;
  static InstMgr instanceManager;
  void* createdObject;

};

#endif



