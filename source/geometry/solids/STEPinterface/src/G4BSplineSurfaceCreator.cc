#include "G4BSplineSurfaceCreator.hh"
#include "G4ControlPoints.hh"
G4BSplineSurfaceCreator G4BSplineSurfaceCreator::csc;

G4BSplineSurfaceCreator::G4BSplineSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4BSplineSurfaceCreator::~G4BSplineSurfaceCreator(){}

void G4BSplineSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  SdaiB_spline_surface* bSpline = new SdaiB_spline_surface(&Ent);
  G4int u,v;
  STEPentity* ent=0;
  G4String attrName("u_degree");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  u = *Attr->ptr.i;
  bSpline->U_degree(*Attr->ptr.i);
  attrName = "v_degree";
  Attr = GetNamedAttribute(attrName, Ent);
  v=*Attr->ptr.i;
  bSpline->V_degree(*Attr->ptr.i);

  attrName = "control_points_list";
  Attr = GetNamedAttribute(attrName, Ent);
  STEPaggregate *Aggr = Attr->ptr.a;
  GenericAggregate* gAggr  =  (GenericAggregate*)Attr->ptr.a;
  bSpline->Control_points_list(gAggr);

  // Get control points
  
  G4int cols,rows;
  cols = v+1;
  rows = u+1;

  STEPentity* entity;
  G4int Index;
  STEPentity *Entity;
  char tmp[16];
  SCLstring s;
  const char *Str = Aggr->asStr(s);

  G4int stringlength = strlen(Str);  
  G4ControlPoints controlPoints(4,rows, cols);
  RealAggregate rationalAggr;
  RealNode* rNode =0;
  for(G4int a=0;a<rows;a++)
    for(G4int b=0;b<cols;b++)    
      {
	// get points
	
	// temp version until the NIST toolkit can handle two dimensional aggregates
	// The string Str contains the STEP file id:s of the underlying point
	// entities so well have to parse the string to get them out...arghhh!
	char c = ' ';
	int Count=0;
	// Loop to find the entities


	// Fill points
	//Temporary solution until the STEP toolkit has been updated:

	while(c != '#')
	  {
	    c = Str[Count];
	    Count++;
	    if(Count>stringlength)
	      {
		G4cout << "\nString index overflow in G4ControlPoints:116";
		exit(0);
	      }
	  }

	c = Str[Count];
	int Index=0;

	while(c != ',' && c != ')')
	  {
	    tmp[Index]=c;
	    Index++;
	    Count++;
	    c = Str[Count];
	  }
	tmp[Index]='\0';
	Index = atoi(tmp);
	//delete [] tmp;
	//Entity = InstanceList.GetSTEPentity(Index);
	MgrNode* MgrTmp = instanceManager.FindFileId(Index);
	Index = instanceManager.GetIndex(MgrTmp);
	Entity = instanceManager.GetSTEPentity(Index);
	void *tmp =G4GeometryTable::CreateObject(*Entity);
	controlPoints.put(a,b,*(G4PointRat*)tmp);
      }  
  
  createdObject = bSpline;
}




void G4BSplineSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{

}

