// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BSplineSurfaceWithKnotsCreator.cc,v 1.4 2000-11-20 18:17:28 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BSplineSurfaceWithKnotsCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include <instmgr.h>

#include "G4BSplineSurfaceWithKnotsCreator.hh"
#include "G4GeometryTable.hh"
#include "G4ControlPoints.hh"
#include "G4KnotVector.hh"
#include "G4BSplineSurface.hh"

G4BSplineSurfaceWithKnotsCreator G4BSplineSurfaceWithKnotsCreator::csc;

G4BSplineSurfaceWithKnotsCreator::G4BSplineSurfaceWithKnotsCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4BSplineSurfaceWithKnotsCreator::~G4BSplineSurfaceWithKnotsCreator() {}

void G4BSplineSurfaceWithKnotsCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4int rows=0, cols=0;
  G4KnotVector* uKnots=0;
  G4KnotVector* vKnots=0;
  G4ControlPoints* controlPoints=0;

  SdaiB_spline_surface_with_knots bSpline(&Ent);
  // U dir
  G4String attrName("u_multiplicities");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  bSpline.u_multiplicities_((IntAggregate*)Attr->ptr.a);
  STEPaggregate *multAggr  = Attr->ptr.a;
  G4int uMultCount = multAggr->EntryCount();
  
  attrName = "u_knots";
  Attr = GetNamedAttribute(attrName, Ent);

  STEPaggregate *knotAggr = Attr->ptr.a;
  G4int uKnotCount = knotAggr->EntryCount();

  G4int totalUKnotCount = 0;
  IntNode* multiNode = (IntNode*)multAggr->GetHead();
  
  G4int a;
  for(a=0;a<uMultCount;a++)
    {
      totalUKnotCount += multiNode->value;
      multiNode = (IntNode*)multiNode->NextNode();
    }
  
  uKnots  = new G4KnotVector(totalUKnotCount);

  RealNode* knotNode = (RealNode*)knotAggr->GetHead();
  multiNode = (IntNode*)multAggr->GetHead();

  bSpline.u_knots_((RealAggregate*)knotAggr);

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

  
  // V dir
  attrName = "v_multiplicities";
  Attr = GetNamedAttribute(attrName, Ent);
  multAggr = Attr->ptr.a;
  bSpline.v_multiplicities_((IntAggregate*)Attr->ptr.a);
  G4int vMultCount = multAggr->EntryCount();
  
  attrName = "v_knots";
  Attr = GetNamedAttribute(attrName, Ent);
  knotAggr = Attr->ptr.a;
  bSpline.v_knots_((RealAggregate*)knotAggr);
  // G4int vKnotCount = knotAggr->EntryCount();

  G4int totalVKnotCount = 0;
  multiNode = (IntNode*)multAggr->GetHead();
  
  for(a=0;a<vMultCount;a++)
    {
      totalVKnotCount += multiNode->value;
      multiNode = (IntNode*)multiNode->NextNode();
    }
  vKnots  = new G4KnotVector(totalVKnotCount);

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

  // b_spline base parts
  
  G4int u=0,v=0;
  attrName = "u_degree";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      u = *Attr->ptr.i;
      bSpline.u_degree_(*Attr->ptr.i);
    }
  else
    G4cerr << "WARNING - G4BSplineSurfaceWithKnotsCreator::CreateG4Geometry" << G4endl
           << "\tSurface attribute u_degree not valid." << G4endl;

  attrName = "v_degree";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      v=*Attr->ptr.i;
      bSpline.v_degree_(*Attr->ptr.i);
    }
  else
    G4cerr << "WARNING - G4BSplineSurfaceWithKnotsCreator::CreateG4Geometry" << G4endl
           << "\tSurface attribute v_degree not valid." << G4endl;

  attrName = "control_points_list";
  Attr = GetNamedAttribute(attrName, Ent);
  if(Attr)
    {
      STEPaggregate *Aggr = Attr->ptr.a;
      GenericAggregate* gAggr  =  (GenericAggregate*)Attr->ptr.a;
      bSpline.control_points_list_(gAggr);
      // Get control points
      
      G4int cols,rows;
      cols = v+1;
      rows = u+1;

      STEPentity *Entity;
      char tmp[16];
      SCLstring s;
      const char *Str = Aggr->asStr(s);

      G4int stringlength = strlen(Str);  
      controlPoints = new G4ControlPoints(4,rows, cols);
      RealAggregate rationalAggr;
      for(G4int a=0;a<rows;a++)
	for(G4int b=0;b<cols;b++)    
	  {
	    // get points
	    
	    // temp version until the NIST toolkit can handle two dimensional aggregates
	    // The string Str contains the STEP file id:s of the underlying point
	    // entities so well have to parse the string to get them out...arghhh!
	    char c = ' ';
	    G4int Count=0;
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
	    G4int Index=0;

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
//          Entity = instanceManager.GetSTEPentity(Index);
            Entity = instanceManager.GetApplication_instance(Index);
	    void *tmp =G4GeometryTable::CreateObject(*Entity);
	    if (tmp)
	      controlPoints->put(a,b,*(G4PointRat*)tmp);
	    else
              G4cerr << "WARNING - G4BSplineSurfaceWithKnotsCreator::CreateG4Geometry" << G4endl
                     << "\tNULL control point (G4PointRat) detected." << G4endl;
	  }  
    }
    else
      G4cerr << "WARNING - G4BSplineSurfaceWithKnotsCreator::CreateG4Geometry" << G4endl
             << "\tSurface attribute control_points_list not valid." << G4endl;

  if (uKnots && vKnots && controlPoints)
  {
    createdObject = new G4BSplineSurface(rows, cols, *uKnots, *vKnots, *controlPoints);
    delete uKnots;
    delete vKnots;
    delete controlPoints;
  }
  else
  {
    createdObject = 0;
    G4cerr << "\tG4BSplineSurface not created !" << G4endl;
  }
   
}

void G4BSplineSurfaceWithKnotsCreator::CreateSTEPGeometry(void* G4obj)
{
}
