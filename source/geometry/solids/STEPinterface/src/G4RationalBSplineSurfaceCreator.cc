#include "G4RationalBSplineSurfaceCreator.hh"

G4RationalBSplineSurfaceCreator G4RationalBSplineSurfaceCreator::csc;

G4RationalBSplineSurfaceCreator::G4RationalBSplineSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4RationalBSplineSurfaceCreator::~G4RationalBSplineSurfaceCreator(){}

void G4RationalBSplineSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("weights_data");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPaggregate *weightAggr = Attr->ptr.a;
  SdaiRational_b_spline_surface* bSpline = new SdaiRational_b_spline_surface(&Ent);
  bSpline->Weights_data((GenericAggregate*)Attr->ptr.a);
 
  // Temp solution until NIST supports
  // two dimensional instances.grrh
  SCLstring s;
  const char* Str = weightAggr->asStr(s);
  G4int rows = weightAggr->EntryCount();   
  G4int cols = 1;
  G4int counter=0;
  while (Str[counter] != ')')
    {
      if(Str[counter]==',')
	cols++;
      counter++;
    }
  
  char c;
  G4int Count=0;
  G4double* ratVector = new G4double[cols*rows];
  for(int a=0;a<rows*cols;a++)
    {
      c = '(';
      while(c == '(' || c==',')
	{
	  Count++;
	  c = Str[Count];
	}
      int Index=0;
      char *tmp = new char[16];
      while(c != ',' && c != ')')
	{
	  tmp[Index]=c;
	  Index++;
	  Count++;
	  c = Str[Count];
	}
      tmp[Index]='\0';
      G4double Value = atof(tmp);
      delete [] tmp;
      ratVector[a] = Value;

    }
  /*
  // base data
  
  // U dir
  attrName = "u_multiplicities";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      bSpline->U_multiplicities((IntAggregate*)Attr->ptr.a);
      STEPaggregate *multAggr = Attr->ptr.a;
      G4int uMultCount = multAggr->EntryCount();
    }

  attrName = "u_knots";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {

      STEPaggregate *knotAggr = Attr->ptr.a;
      G4int uKnotCount = knotAggr->EntryCount();

      G4int totalUKnotCount = 0;
      IntNode* multiNode = (IntNode*)multAggr->GetHead();
      
      for(G4int a=0;a<uMultCount;a++)
	{
	  totalUKnotCount += multiNode->value;
	  multiNode = (IntNode*)multiNode->NextNode();
	}
      
      G4KnotVector *uKnots  = new G4KnotVector(totalUKnotCount);

      RealNode* knotNode = (RealNode*)knotAggr->GetHead();
      multiNode = (IntNode*)multAggr->GetHead();

      bSpline->U_knots((RealAggregate*)knotAggr);
      


      G4int multValue=0;
      G4double knotValue=0;
      G4int index=0;
      for(a=0;a<uKnotCount;a++)
	{
	  multValue = multiNode->value;
	  knotValue = knotNode->value;

	  for(G4int b=0;b<multValue;b++)
	    {
	      uKnots->PutKnot(index, knotValue);
	      index++;
	    }
	  knotNode = (RealNode*)knotNode->NextNode();
	}

    }
  // V dir
  attrName = "v_multiplicities";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      multAggr = Attr->ptr.a;
      bSpline->U_multiplicities((IntAggregate*)Attr->ptr.a);
      G4int vMultCount = multAggr->EntryCount();
    }
  attrName = "v_knots";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      knotAggr = Attr->ptr.a;
      bSpline->V_knots((RealAggregate*)knotAggr);
      G4int vKnotCount = knotAggr->EntryCount();

      G4int totalVKnotCount = 0;
      multiNode = (IntNode*)multAggr->GetHead();
      
      for(a=0;a<uMultCount;a++)
	{
	  totalVKnotCount += multiNode->value;
	  multiNode = (IntNode*)multiNode->NextNode();
	}
      G4KnotVector *vKnots  = new G4KnotVector(totalVKnotCount);

      knotNode = (RealNode*)knotAggr->GetHead();
      multiNode = (IntNode*)multAggr->GetHead();

      multValue=0;
      knotValue=0;
      index=0;

      for(a=0;a<uKnotCount;a++)
	{
	  multValue = multiNode->value;
	  knotValue = knotNode->value;

	  for(G4int b=0;b<multValue;b++)
	    {
	      vKnots->PutKnot(index, knotValue);
	      index++;
	    }
	  knotNode = (RealNode*)knotNode->NextNode();
	}
    }
  // b_spline base parts
  
  G4int u,v;
  STEPentity* ent=0;
  G4String attrName("u_degree");
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      u = *Attr->ptr.i;
      bSpline->U_degree(*Attr->ptr.i);
    }


  attrName = "v_degree";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      v=*Attr->ptr.i;
      bSpline->V_degree(*Attr->ptr.i);
    }

  attrName = "control_points_list";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
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
    }  
    */


  
  createdObject = bSpline;
}

void G4RationalBSplineSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{

}


