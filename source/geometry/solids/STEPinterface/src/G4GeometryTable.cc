#include "G4GeometryTable.hh"
#include "includes"


G4GeometryTable G4GeometryTable::gt;
G4CreatorVector G4GeometryTable::RegisteredObjects;

G4GeometryTable::G4GeometryTable()
{
}

G4GeometryTable::~G4GeometryTable()
{

}

G4bool G4GeometryTable::ExistsInTable(G4String& objectName)
{
  G4bool existsInTable = false;

  for(G4int a=0;a<RegisteredObjects.length();a++)
  {
    G4GeometryCreator* gcTmp =  RegisteredObjects(a);
    G4String rName = gcTmp->Name();
    if(rName == objectName) existsInTable = true;
  }

  return existsInTable;
}

void G4GeometryTable::RegisterObject(G4GeometryCreator* gc)
{
  G4String newName = gc->Name();
  
  if(!ExistsInTable(newName))
    RegisteredObjects.append(gc);
}

G4GeometryCreator* G4GeometryTable::GetObject(G4String objectName)
{
  for(G4int a=0; a<RegisteredObjects.length(); a++)
  {
    G4GeometryCreator* gcTmp =  RegisteredObjects(a);
    G4String rName = gcTmp->Name();
    
    if(rName == objectName)
      return gcTmp;
  }
  
  return (G4GeometryCreator*)0;
}


void* G4GeometryTable::CreateObject(STEPentity& Ent)
{
  Ent.ResetAttributes();
  G4String name = Ent.EntityName();
  G4GeometryCreator* gctmp = GetObject(name);

  if(gctmp)
  {
    gctmp->CreateG4Geometry(Ent);
    void* obj = gctmp->GetCreatedObject();

//#define G4_STEPINTERFACE_DEBUG 1 
#ifdef G4_STEPINTERFACE_DEBUG
    if (obj == 0) 
      G4cout << "G4 object of type " << gctmp->Name() 
	     << " not necessary" << endl;
    else 
      G4cout << "Created G4 object of type " << gctmp->Name() << endl;
#endif

    return obj;
  }
  else
  {
    G4String err = "Geometry creator for entity " + name +  " not found.";
    G4cout << err << endl;
    //      G4Exception(err);
    return 0;
  } 
}


void* G4GeometryTable::CreateSTEPObject(void* G4obj, G4String& objName)
{

  G4GeometryCreator* gctmp = GetObject(objName);

  if(gctmp)
    {
      gctmp->CreateSTEPGeometry(G4obj);
      void* obj = gctmp->GetCreatedObject();
      return obj;
    }
  else
    {
      G4String err = "\nSTEP creator for entity " + objName +  " not found.\n";
      G4cout << err;
      //      G4Exception(err);
      return 0;
    }

}


void G4GeometryTable::PrintObjectNames()
{
  G4String name;
  G4GeometryCreator* gc=0;
  G4int a=0;

  G4int objectCount = RegisteredObjects.length();
  if(!objectCount)
    {
      G4cout << "No registered objects in table." << endl;
      return;
    }

  G4cout << "Objects in table:" << endl;  
  
  for(a=0;a<objectCount;a++)
    {
      gc = RegisteredObjects[a];
      name = gc->Name();
      G4cout << name << endl;
    }

}


G4PointCreator a;
G4CartesianPointCreator b;
G4PointOnCurveCreator c;
G4PointOnSurfaceCreator d;
G4PointReplicaCreator e;
G4VertexPointCreator f;

G4CurveCreator c1;
G4LineCreator c2;
G4CircleCreator c3;
G4EllipseCreator c4;
G4HyperbolaCreator c5;
G4ParabolaCreator c6;

G4EdgeCurveCreator c7;
G4BSplineCurveWithKnotsCreator bsc1;

G4EdgeLoopCreator l1;
G4OrientedEdgeCreator e1;
G4FaceBoundCreator b1;
G4FaceOuterBoundCreator b2 ;

G4DirectionCreator dir;
G4VectorCreator vec;

G4Axis1PlacementCreator ax1;
G4Axis2Placement2dCreator ax2;
G4Axis2Placement3dCreator ax3;
G4Axis2PlacementCreator ax4;
G4Axis2PlacementsCreator ax5;

G4AdvancedFaceCreator s2;
G4PlaneCreator s3;
G4CylindricalSurfaceCreator s4;
G4FaceSurfaceCreator fsc1;
G4ToroidalSurfaceCreator s5;
G4SphericalSurfaceCreator s6;
G4ConicalSurfaceCreator s7;
G4BSplineSurfaceCreator s8;
G4BoundedSurfaceCreator s9;
G4RationalBSplineSurfaceCreator s10;
G4BSplineSurfaceWithKnotsCreator s11;

G4ClosedShellCreator sc1;
G4OpenShellCreator sc2;

G4AdvancedBrepShapeRepresentationCreator shape1;
G4ManifoldSolidBrepCreator ms1;

G4ConnectedFaceSetCreator cfsc1;

G4ItemDefinedTransformationCreator trans1;

G4RepresentationRelationshipCreator reprel1;
G4RepresentationRelationshipWithTransformationCreator repretrans;

G4ContextDependentShapeRepresentationCreator condep1;

G4ProductDefinitionShapeCreator pdsc;
G4ShapeRepresentationRelationshipCreator srrc;
G4ShapeRepresentationCreator srcs;

G4ShapeDefinitionRepresentationCreator sdr;
G4GeometricRepresentationContextCreator grc;







