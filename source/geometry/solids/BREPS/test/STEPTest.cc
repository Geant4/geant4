//////////////////////////////////////////////////////////////////////////
// $Id: STEPTest.cc,v 1.5 2000-08-28 08:58:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4Timer.hh"

#include "G4AssemblyCreator.hh"
#include "G4Assembly.hh"
#include "instmgr.h"
#include "globals.hh"


G4int main()
{

  G4cout<<"//////////////////////////////////////////////////////////////\n";
  G4cout<<"// STEP assembly creation\n\n\n";

  // Define the STEP file, 
  // and the STEP file reader used is G4NISTStepReader 
  G4cout<<" Choose the step file to be used : \n";
  G4cout<<"    (1) G4cyl_v.stp      (volume)  \n";
  G4cout<<"    (2) G4sphere_v.stp   (volume)  \n";
  G4cout<<"    (3) G4cone_v.stp     (volume)  \n";
  G4cout<<"    (4) G4boite_v.stp    (volume)  \n";
  G4cout<<"    (5) G4assembly_v.stp (volume)  \n";
  G4cout<<"    (6) G4rod_solid.stp  (ProE)    \n";
  G4cout<<"    (7) G4rod_solid2.stp (SldWrks) \n";
  G4cout<<"    (8) G4rod_place_asm.stp        \n";
  G4cout<<"    (9) G4spline1.stp                :";

  G4String stepfile;
  G4int stepf;
  G4cin>>stepf;
  
  if(stepf == 1)
    stepfile = "stepfiles/G4cyl_v.stp";   
  else if(stepf == 2)
    stepfile = "stepfiles/G4sphere_v.stp"; 
  else if(stepf == 3)
    stepfile = "stepfiles/G4cone_v.stp"; 
  else if(stepf == 4)
    stepfile = "stepfiles/G4boite_v.stp";
  else if(stepf == 5)
    stepfile = "stepfiles/G4assembly_v.stp";
  else if(stepf == 6)
    stepfile = "stepfiles/G4rod_solid.stp";
  else if(stepf == 7)
    stepfile = "stepfiles/G4rod_solid2.stp";
  else if(stepf == 8)
    stepfile = "stepfiles/G4rod_place_asm.stp";
  else if(stepf == 9)
    stepfile = "stepfiles/G4spline1.stp";

  G4AssemblyCreator MyAC(stepfile);


  // Read the STEP file and put the instances before 
  // into G4NISTStepReader.InstanceList (which is a InstMgr) and after 
  // into G4GeometryCreator.instanceManager (which is a InstMgr)
  MyAC.ReadStepFile();

  G4int Advbrepshapes = MyAC.instanceManager.
    EntityKeywordCount("Advanced_Brep_Shape_Representation");

  G4int Condepshapes = MyAC.instanceManager.
    EntityKeywordCount("Context_Dependent_Shape_Representation");
  
  G4int Shapedefreps = MyAC.instanceManager.
    EntityKeywordCount("Shape_Definition_Representation");   


  G4cout<<"\n Advanced_Brep_Shape_Representation = "<<Advbrepshapes<<G4endl;
  G4cout<<"\n Context_Dependent_Shape_Representation = "<<Condepshapes<<G4endl;
  G4cout<<"\n Shape_Definition_Representation = "<<Shapedefreps<<G4endl;
  
  // Create the STEP entities into psv (which is a G4PlacedSolidVector*)
  // and after copy into G4GeometryCreator.createdObject (which is a void*)
  STEPentity* ent =0;
  MyAC.CreateG4Geometry(*ent);

  G4cout<<"\n\n STEP assembly created\n\n";


  G4cout<<"///////////////////////////////////////////////////////////////\n";
  G4cout<<"// Test the STEP assembly \n\n\n";

  // Chek that reader output entities equal to file


  // Get the created objects into tmp (which is a void*)
  void *tmp =  MyAC.GetCreatedObject();

  // Set placed vectors into assembly.placedVec (which is a G4PlacedVector)
  G4Assembly* assembly = new G4Assembly();
  assembly->SetPlacedVector(*(G4PlacedVector*)tmp);

  
  G4PlacedSolid* ps = 0;
  G4int solids = assembly->GetNumberOfSolids();
  G4cout<<"Number of solids : "<<solids<<G4endl;

  // Check that BREP solids & surfaces build 
  // the solid specified by the reader output
  
  G4Point3D Pt[4];
  Pt[0] = G4Point3D(   99,   99,  100);
  Pt[1] = G4Point3D( 1000, 1000,  150);
  Pt[2] = G4Point3D( -150, -150,-5000);
  Pt[3] = G4Point3D(    1,    1,   50);
  
  G4Vector3D Dir[4];
  Dir[0] = G4Vector3D(    1,    0,    0);
  Dir[1] = G4Vector3D(   -1,   -1,    0);
  Dir[2] = G4Vector3D(    0,    0,    1);
  Dir[3] = G4Vector3D(    1,    1,   -1);
  
  EInside in[4];
  G4double dist[4][3];

  for(G4int c=0; c<solids; c++)
  {
    G4cout<<"\n=============================================";
    G4cout<<"\n\n Get placed solid "<<c+1<<"/"<<solids;
    ps = assembly->GetPlacedSolid(c);

    G4cout<<"\n\n Get solid";
    G4BREPSolid* sol = (G4BREPSolid*)ps->GetSolid();
    G4cout<<"\n  solid bbox is";
    G4cout<<"\n     bbox min :"
	  <<" x="<<sol->GetBBox()->GetBoxMin().x()
	  <<" y="<<sol->GetBBox()->GetBoxMin().y()
	  <<" z="<<sol->GetBBox()->GetBoxMin().z();
    G4cout<<"\n     bbox max :"
	  <<" x="<<sol->GetBBox()->GetBoxMax().x()
	  <<" y="<<sol->GetBBox()->GetBoxMax().y()
	  <<" z="<<sol->GetBBox()->GetBoxMax().z()<<G4endl;

    G4cout<<"\n Get Translation"<<G4endl;
    G4ThreeVector* tr = ps->GetTranslation();
    G4cout<<"   x="<<tr->x()<<" y="<<tr->y()<<" z="<<tr->z();
  
    G4cout<<"\n\n Get Rotation"<<G4endl;
    HepRotation* hr = ps->GetRotation();
    G4cout<<"   xx="<<hr->xx()<<" xy="<<hr->xy()<<" xz="<<hr->xz()<<G4endl;
    G4cout<<"   yx="<<hr->yx()<<" yy="<<hr->yy()<<" yz="<<hr->yz()<<G4endl;
    G4cout<<"   zx="<<hr->zx()<<" zy="<<hr->zy()<<" zz="<<hr->zz()<<G4endl;


    // -> Check methods :
    //  - Inside
    //  - DistanceToIn
    //  - DistanceToOut
    
    G4cout<<G4endl<<G4endl;
    
    for (G4int a = 0; a<4; a++)
    {
      in[a] = sol->Inside(Pt[a]);
      
      G4cout<<"----------------------------"<<G4endl;
      
      G4cout<<"x="<<Pt[a].x()
	    <<"  y="<<Pt[a].y()<<"  z="<<Pt[a].z();
      
      if( in[a] == kInside )
      {
	G4cout <<" is inside"<<G4endl;
	
	dist[a][1] = sol->DistanceToOut(Pt[a]);
	G4cout<<"\nDistance to out is :"<<dist[a][1]<<G4endl;
	
	G4cout << "\nDir   : x=" << Dir[a].x() 
	     << " y=" << Dir[a].y() 
	     << " z=" << Dir[a].z()<<G4endl;
	dist[a][2] = sol->DistanceToOut(Pt[a], Dir[a]);
	G4cout<<"Distance to out is :"<<dist[a][2]<<G4endl;
      }
      else
      {
	G4cout <<" is outside"<<G4endl;

	dist[a][1] = sol->DistanceToIn(Pt[a]);
	G4cout<<"\nDistance to in is :"<<dist[a][1];
	
	G4cout << "\nDir   : x=" << Dir[a].x() 
	     << " y=" << Dir[a].y() 
	     << " z=" << Dir[a].z()<<G4endl;
	dist[a][2] = sol->DistanceToIn(Pt[a], Dir[a]);
	G4cout<<"Distance to in is :"<<dist[a][2];
      }
      
      G4cout<<G4endl;
    }  
  }

  return 0;
}
