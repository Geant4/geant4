// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BoundedSurfaceCreator.cc,v 1.3 2000-01-21 13:45:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BoundedSurfaceCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include <instmgr.h>
#include <STEPcomplex.h>

#include "G4BoundedSurfaceCreator.hh"
#include "G4GeometryTable.hh"
#include "G4ControlPoints.hh"
#include "G4BSplineSurface.hh"

G4BoundedSurfaceCreator G4BoundedSurfaceCreator::csc;

G4BoundedSurfaceCreator::G4BoundedSurfaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4BoundedSurfaceCreator::~G4BoundedSurfaceCreator() {}

void G4BoundedSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  STEPcomplex* complexEnt = (STEPcomplex*)&Ent;
  void* tmp=0;
  STEPentity* subEnt=0;
  SdaiB_spline_surface *bSpline=0;
  SdaiB_spline_surface_with_knots *bSplineWithKnots=0;
  SdaiRational_b_spline_surface *rationalBSpline =0;

  if(complexEnt->EntityExists("B_Spline_Surface"))
    {
      subEnt = complexEnt->EntityPart("B_Spline_Surface");
      bSpline =(SdaiB_spline_surface*)G4GeometryTable::CreateObject(*subEnt);      
    }
  
  if(complexEnt->EntityExists("B_Spline_Surface_With_Knots"))
    {
      subEnt = complexEnt->EntityPart("B_Spline_Surface_With_Knots");
      bSplineWithKnots =(SdaiB_spline_surface_with_knots*)G4GeometryTable::CreateObject(*subEnt);
    }

  if(complexEnt->EntityExists("Rational_B_Spline_Surface"))
    {
      subEnt = complexEnt->EntityPart("Rational_B_Spline_Surface");
      rationalBSpline =(SdaiRational_b_spline_surface*)G4GeometryTable::CreateObject(*subEnt);
    }


  G4int u,v;
  u=bSpline->u_degree_();
  v=bSpline->v_degree_();


  // Get control points
  
  G4int cols,rows;
  cols = v+1;
  rows = u+1;
  char tmpstr[16];
  STEPentity* entity;
  G4int Index;
  STEPentity *Entity;
  SCLstring s;
  STEPaggregate *Aggr=bSpline->control_points_list_();
  const char *Str = Aggr->asStr(s);

  G4int stringlength = strlen(Str);  
  G4ControlPoints controlPoints(4,rows, cols);
  RealAggregate rationalAggr;
  RealNode* rNode =0;
  G4int a;
  for(a=0;a<rows;a++)
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
	    tmpstr[Index]=c;
	    Index++;
	    Count++;
	    c = Str[Count];
	  }
	tmpstr[Index]='\0';
	Index = atoi(tmpstr);
	//delete [] tmp;
	//Entity = InstanceList.GetSTEPentity(Index);
	MgrNode* MgrTmp = instanceManager.FindFileId(Index);
	Index = instanceManager.GetIndex(MgrTmp);
//      Entity = instanceManager.GetSTEPentity(Index);
        Entity = instanceManager.GetApplication_instance(Index);
	void *tmp =G4GeometryTable::CreateObject(*Entity);
	controlPoints.put(a,b,*(G4PointRat*)tmp);
      }  
  
  
  // Get knot vectors
  STEPaggregate *multAggr = bSplineWithKnots->u_multiplicities_();
  G4int uMultCount = multAggr->EntryCount();
  STEPaggregate *knotAggr = bSplineWithKnots->u_knots_();
  G4int uKnotCount = knotAggr->EntryCount();

  G4int totalUKnotCount = 0;
  IntNode* multiNode = (IntNode*)multAggr->GetHead();
  
  for(a=0;a<uMultCount;a++)
    {
      totalUKnotCount += multiNode->value;
      multiNode = (IntNode*)multiNode->NextNode();
    }
  
  G4KnotVector *uKnots  = new G4KnotVector(totalUKnotCount);

  RealNode* knotNode = (RealNode*)knotAggr->GetHead();
  multiNode = (IntNode*)multAggr->GetHead();
  
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
  multAggr = bSplineWithKnots->v_multiplicities_();
  G4int vMultCount = multAggr->EntryCount();
  
  knotAggr = bSplineWithKnots->v_knots_();
  G4int vKnotCount = knotAggr->EntryCount();

  G4int totalVKnotCount = 0;
  multiNode = (IntNode*)multAggr->GetHead();
  
  for(a=0;a<vMultCount;a++)
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



  // copy weights data
  STEPaggregate *weightAggr =  rationalBSpline->weights_data_();

  // Temp solution until NIST supports
  // two dimensional instances.grrh
  const char* Str2 = weightAggr->asStr(s);
  rows = weightAggr->EntryCount();   
  cols = 1;
  G4int counter=0;
  while (Str2[counter] != ')')
    {
      if(Str2[counter]==',')
	cols++;
      counter++;
    }
  
  char c;
  G4int Count=0;
  G4double* ratVector = new G4double[cols*rows];
  for(a=0;a<rows*cols;a++)
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

  controlPoints.SetWeights(ratVector);

  
  // create BSpline
  G4BSplineSurface* bSplineSrf = new G4BSplineSurface(
						   u,
						   v,
						   *uKnots,
						   *vKnots,
						   controlPoints
						     );  
}

void G4BoundedSurfaceCreator::CreateSTEPGeometry(void * G4obj)
{
}
