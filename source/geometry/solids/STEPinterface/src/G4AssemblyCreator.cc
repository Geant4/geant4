// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AssemblyCreator.cc,v 1.6 2000-11-09 16:35:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4AssemblyCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4AssemblyCreator.hh"
#include "G4Assembly.hh"

#include "G4GeometryTable.hh"
#include "G4StepFileReader.hh"
#include "G4NISTStepReader.hh"
#include "G4ToroidalSurface.hh"

extern void HeaderSchemaInit (Registry & reg);

G4AssemblyCreator G4AssemblyCreator::ci;

G4AssemblyCreator::G4AssemblyCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4AssemblyCreator::G4AssemblyCreator(const G4AssemblyCreator& c)
{
  index = c.index;
  StepReader = c.StepReader;
  STEPfileName = c.STEPfileName;
}

G4AssemblyCreator& G4AssemblyCreator::operator=(const G4AssemblyCreator& c)
{
  if (&c == this) return *this;
  
  index = c.index;
  StepReader = c.StepReader;
  STEPfileName = c.STEPfileName;
  
  return *this;
}

G4AssemblyCreator::G4AssemblyCreator(const G4String& fileName,
                                     const G4String& readerName)
{
  // Name of STEP file to read in, existance of file should be checked...
  //
  STEPfileName = fileName;

  // define type of STEP file reader to use (use else...)
  //
  if (readerName == "NIST")  StepReader = new G4NISTStepReader();
  
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

void G4AssemblyCreator::CreateG4Geometry(STEPentity& sEntity)
{
  G4cout << G4endl
         << "Creating assembly of G4 solids ..." << G4endl;
  
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
  int tmpindex = 0;
  
  void *tmp = 0;
  G4PlacedSolidVector* psv = new G4PlacedSolidVector();
  G4PlacedSolid*       ps  = new G4PlacedSolid(); 
  G4int a;

  // L. Broglia
  if(ConDepShapes>0)
  {
    index = 0;
    const char* keyw = "Context_Dependent_Shape_Representation";

#ifdef G4_STEPINTERFACE_DEBUG
    G4cout << "\n\n Creating the " << keyw << G4endl;
#endif
  
    for( a=0; a< ConDepShapes; a++)
    {
#ifdef G4_STEPINTERFACE_DEBUG
      G4cout << "loop " << a+1 << " of " << ConDepShapes << G4endl;
#endif     

      ent = instanceManager.GetApplication_instance(keyw, index);

      if(ent!= ENTITY_NULL)
      {
        // Be careful, tmpindex not correspond to STEPfile_id !
        tmpindex = instanceManager.GetIndex(ent);

      	tmp = G4GeometryTable::CreateObject(*ent);	
	if (!tmp)
	{
	   G4cerr << "--- GeometryTable::CreateG4Geometry()" << G4endl;
	   G4cerr << "    No correspondent G4 geometry created for entity: "
	          << ent->EntityName() << G4endl;
	}
	else
	{
	  G4PlacedSolidVector* tmpV = (G4PlacedSolidVector*)tmp; 
	  G4int entr = tmpV->entries();

	  for(G4int b=0; b<entr; b++)
	  {
	    ps =  tmpV->at(b);
	    psv->append(ps);
	  }
	}

	index = ent->STEPfile_id ;

#ifdef G4_STEPINTERFACE_DEBUG
	G4cout << keyw << " found in index " << index << G4endl;
#endif

      }
      // Set index to the true value
      index = tmpindex + 1;
      ent=0;
    }

  }
  else if (ShapeDefReps>0)
  {       

    const char* key = "Shape_Definition_Representation";
    
#ifdef G4_STEPINTERFACE_DEBUG
    G4cout << "\n Creating the " << key << G4endl;
#endif

    for(a=0; a<  ShapeDefReps ; a++)
    {
#ifdef G4_STEPINTERFACE_DEBUG
      G4cout << "loop " << a+1 << " of " << ShapeDefReps << G4endl;
#endif
      
      ent = instanceManager.GetApplication_instance(key, index);

      if(ent!= ENTITY_NULL)
      {
        // Be careful, tmpindex not correspond to STEPfile_id !
        tmpindex = instanceManager.GetIndex(ent);
    
	tmp = G4GeometryTable::CreateObject(*ent);
	if (!tmp)
	{
	   G4cerr << "--- GeometryTable::CreateG4Geometry()" << G4endl;
	   G4cerr << "    No correspondent G4 geometry created for entity: "
	          << ent->EntityName() << G4endl;
	}
	else
	{
	  G4PlacedSolidVector* tmpV = (G4PlacedSolidVector*)tmp; 
	  G4int entr = tmpV->entries();

	  for(G4int b=0; b<entr; b++)
	  {
	    ps =  tmpV->at(b);
	    psv->append(ps);
	  }
	}
	index = ent->STEPfile_id ;

#ifdef G4_STEPINTERFACE_DEBUG
	G4cout << key << " found in index "<< index << G4endl;
#endif

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
