#include "includes"
#include "G4AssemblyCreator.hh"

int main(int argc,char** argv) 
{
  
  G4PointCreator::GetInstance();  
  G4CartesianPointCreator::GetInstance();
  G4PointOnCurveCreator::GetInstance();
  G4PointOnSurfaceCreator::GetInstance();  
  G4PointReplicaCreator::GetInstance();
  G4VertexPointCreator::GetInstance();
  
  G4CurveCreator::GetInstance();
  G4AssemblyCreator::GetInstance();
  G4ContextDependentShapeRepresentationCreator::GetInstance();
  G4ShapeDefinitionRepresentationCreator::GetInstance();
  
  char *filename;
  if ( argc == 1 ) {
    filename = "stepfiles/G4rod_place_asm.stp";
  } else {
    filename = argv[1];
  }
  G4AssemblyCreator MyAC(filename);

  //  G4AssemblyCreator MyAC("stepfiles/CMSTracker.Step");
  //G4AssemblyCreator MyAC("stepfiles/coil1.Step");    
  MyAC.ReadStepFile();
  STEPentity *ent=0;
  MyAC.CreateG4Geometry(*ent);


  
    /*
       // G4Placement* place = new G4Placement();
    //MyAC.CreateSTEPGeometry(place);    

    G4ThreeVec o(0,0,0);
    G4ThreeVec a(1,0,0);
    G4ThreeVec d(0,1,0);
    //    G4FCylindricalSurface *fCyl = new G4FCylindricalSurface(o,a,100,50);
    //    MyAC.CreateSTEPGeometry(fCyl);

    G4ToroidalSurface* tor = new G4ToroidalSurface(o,a,d,10,100);
    MyAC.CreateSTEPGeometry(tor);
    */
    //  G4GeometryTable::PrintObjectNames();

  
}

