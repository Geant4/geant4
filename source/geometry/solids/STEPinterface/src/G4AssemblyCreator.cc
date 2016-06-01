#include "G4AssemblyCreator.hh"
#include "G4Assembly.hh"
extern void HeaderSchemaInit (Registry & reg);

G4AssemblyCreator G4AssemblyCreator::ci;

G4AssemblyCreator::G4AssemblyCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4AssemblyCreator::G4AssemblyCreator(G4String fileName, G4String readerName)
{
  // Name of STEP file to read in, existance of file should be checked
  //here also...
  STEPfileName = fileName;

  // define type of STEP file reader to use
  if(readerName == "NIST")
    StepReader = new G4NISTStepReader();
  //else...
  
}

G4AssemblyCreator::~G4AssemblyCreator()
{
  delete StepReader;
}


void G4AssemblyCreator::ReadStepFile()
{
  StepReader->ReadSTEPFile(STEPfileName);
  instanceManager = StepReader->GetInstanceManager();
}

void G4AssemblyCreator::CreateG4Geometry(STEPentity& Ent)
{  
  // Advanced_Brep_Shape_Representation are created into
  // Context_Dependent_Shape_Representation and
  // Shape_Definition_Representation
  G4int AdvancedBrepShapes = instanceManager.EntityKeywordCount
    ("Advanced_Brep_Shape_Representation");
  
  
  G4int ConDepShapes = instanceManager.EntityKeywordCount
    ("Context_Dependent_Shape_Representation");
  
  G4int ShapeDefReps = instanceManager.EntityKeywordCount
    ("Shape_Definition_Representation"); 
  
  G4int instanceCount = instanceManager.InstanceCount();
  
  STEPentity* ent=0;
  index = 0;
  int tmpindex;
  
  void *tmp = 0;
  G4PlacedSolidVector* psv = new G4PlacedSolidVector();
  G4PlacedSolid*       ps  = new G4PlacedSolid(); 
  G4int a;

  // L. Broglia
  // Total assumption ! 
  
  if(ConDepShapes>0)
  {
    G4cout<<"\n\n Creating the Context_Dependent_Shape_Representation"<<endl;
    index = 0;
  
    for( a=0; a< ConDepShapes; a++)
    {
      G4cout<<"loop "<<a+1<<" of "<<ConDepShapes<<endl;
      
      // Be careful, tmpindex not correspond to STEPfile_id !
      tmpindex = 
	instanceManager.GetIndex("Context_Dependent_Shape_Representation",
				 index                                    );
      
      ent = instanceManager.GetSTEPentity(tmpindex);
      
      if(ent!= ENTITY_NULL)
      {
	tmp =G4GeometryTable::CreateObject(*ent);
	
	G4PlacedSolidVector* tmpV = (G4PlacedSolidVector*)tmp; 
	G4int entr = tmpV->entries();

	for(G4int b=0; b<entr; b++)
	{
	  ps =  tmpV->at(b);
	  psv->append(ps);
	}

	index = ent->STEPfile_id ;
	G4cout<<" Context_Dependent_Shape_Representation find in index "
	      <<index<<endl;
      }
      
      // Set index to the true value
      index = tmpindex + 1;
      ent=0;
    }

  }
  else
  {       
    G4cout<<"\n Creating the Shape_Definition_Representation"<<endl;
    for(a=0; a<  ShapeDefReps ; a++)
    {
      G4cout<<"loop "<<a+1<<" of "<<ShapeDefReps<<endl;
      
      // Be careful, tmpindex not correspond to STEPfile_id !
      tmpindex = instanceManager.GetIndex("Shape_Definition_Representation",
					  index                             );
      
      ent = instanceManager.GetSTEPentity(tmpindex);
    
      if(ent!= ENTITY_NULL)
      {
	tmp = G4GeometryTable::CreateObject(*ent);

	G4PlacedSolidVector* tmpV = (G4PlacedSolidVector*)tmp; 
	G4int entr = tmpV->entries();

	for(G4int b=0; b<entr; b++)
	{
	  ps =  tmpV->at(b);
	  psv->append(ps);
	}

	index = ent->STEPfile_id ;
	G4cout<<"  Shape_Definition_Representation find in index "<<index<<endl;
      }
      
      // Set index to the true value
      index = tmpindex + 1;    
      ent=0;
    }
  } 
  
  createdObject = psv;
}


void G4AssemblyCreator::CreateSTEPGeometry(void* G4obj)
{
  Registry reg(&::HeaderSchemaInit);
  SdaiCONFIG_CONTROL_DESIGNInit (reg);

  //G4Placement* plc = (G4Placement*)G4obj;
  //G4String name("Axis2Placement3d");
  //  G4FCylindricalSurface *fCyl = (G4FCylindricalSurface *)G4obj;
  //  G4String name("Cylindrical_Surface");
  G4ToroidalSurface *tor = (G4ToroidalSurface *)G4obj;
  G4String name("Toroidal_Surface");
  {
    void *tmp =G4GeometryTable::CreateSTEPObject(G4obj, name);
  }
}




