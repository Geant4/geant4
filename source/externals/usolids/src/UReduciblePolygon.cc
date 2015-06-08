//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UReduciblePolygon
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include <iostream>

#include "UReduciblePolygon.hh"

//
// Constructor: with simple arrays
//
UReduciblePolygon::UReduciblePolygon(const double a[],
                                     const double b[],
                                     int n)
  : aMin(0.), aMax(0.), bMin(0.), bMax(0.),
    vertexHead(0)
{
  //
  // Do all of the real work in Create
  //
  Create(a, b, n);
}


//
// Constructor: special PGON/PCON case
//
UReduciblePolygon::UReduciblePolygon(const double rmin[],
                                     const double rmax[],
                                     const double z[], int n)
  : aMin(0.), aMax(0.), bMin(0.), bMax(0.),
    vertexHead(0)
{
  //
  // Translate
  //
  double* a = new double[n * 2];
  double* b = new double[n * 2];

  double* rOut = a + n,
          *zOut = b + n,
           *rIn = rOut - 1,
            *zIn = zOut - 1;

  int i;
  for (i = 0; i < n; i++, rOut++, zOut++, rIn--, zIn--)
  {
    *rOut = rmax[i];
    *rIn  = rmin[i];
    *zOut = *zIn = z[i];
  }

  Create(a, b, n * 2);

  delete [] a;
  delete [] b;
}


//
// Create
//
// To be called by constructors, fill in the list and statistics for a new
// polygon
//
void UReduciblePolygon::Create(const double a[],
                               const double b[], int n)
{
  if (n < 3)
    UUtils::Exception("UReduciblePolygon::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, "Less than 3 vertices specified.");

  const double* anext = a, *bnext = b;

  ABVertex* prev = 0;
  do
  {
    ABVertex* newVertex = new ABVertex;
    newVertex->a = *anext;
    newVertex->b = *bnext;
    newVertex->next = 0;
    if (prev == 0)
    {
      vertexHead = newVertex;
    }
    else
    {
      prev->next = newVertex;
    }

    prev = newVertex;
  }
  while (++anext, ++bnext < b + n);

  numVertices = n;

  CalculateMaxMin();
}



//
// Destructor
//
UReduciblePolygon::~UReduciblePolygon()
{
  ABVertex* curr = vertexHead;
  while (curr)
  {
    ABVertex* toDelete = curr;
    curr = curr->next;
    delete toDelete;
  }
}


//
// CopyVertices
//
// Copy contents into simple linear arrays.
// ***** CAUTION ***** Be care to declare the arrays to a large
// enough size!
//
void UReduciblePolygon::CopyVertices(double a[], double b[]) const
{
  double* anext = a, *bnext = b;
  ABVertex* curr = vertexHead;
  while (curr)
  {
    *anext++ = curr->a;
    *bnext++ = curr->b;
    curr = curr->next;
  }
}


//
// ScaleA
//
// Multiply all a values by a common scale
//
void UReduciblePolygon::ScaleA(double scale)
{
  ABVertex* curr = vertexHead;
  while (curr)
  {
    curr->a *= scale;
    curr = curr->next;
  }
}


//
// ScaleB
//
// Multiply all b values by a common scale
//
void UReduciblePolygon::ScaleB(double scale)
{
  ABVertex* curr = vertexHead;
  while (curr)
  {
    curr->b *= scale;
    curr = curr->next;
  }
}


//
// RemoveDuplicateVertices
//
// Remove adjacent vertices that are equal. Returns "false" if there
// is a problem (too few vertices remaining).
//
bool UReduciblePolygon::RemoveDuplicateVertices(double tolerance)
{
  ABVertex* curr = vertexHead,
            *prev = 0, *next = 0;
  while (curr)
  {
    next = curr->next;
    if (next == 0) next = vertexHead;

    if (std::fabs(curr->a - next->a) < tolerance &&
        std::fabs(curr->b - next->b) < tolerance)
    {
      //
      // Duplicate found: do we have > 3 vertices?
      //
      if (numVertices <= 3)
      {
        CalculateMaxMin();
        return false;
      }

      //
      // Delete
      //
      ABVertex* toDelete = curr;
      curr = curr->next;
      delete toDelete;

      numVertices--;

      if (prev) prev->next = curr;
      else vertexHead = curr;
    }
    else
    {
      prev = curr;
      curr = curr->next;
    }
  }

  //
  // In principle, this is not needed, but why not just play it safe?
  //
  CalculateMaxMin();

  return true;
}


//
// RemoveRedundantVertices
//
// Remove any unneeded vertices, i.e. those vertices which
// are on the line connecting the previous and next vertices.
//
bool UReduciblePolygon::RemoveRedundantVertices(double tolerance)
{
  //
  // Under these circumstances, we can quit now!
  //
  if (numVertices <= 2) return false;

  double tolerance2 = tolerance * tolerance;

  //
  // Loop over all vertices
  //
  ABVertex* curr = vertexHead, *next = 0;
  while (curr)
  {
    next = curr->next;
    if (next == 0) next = vertexHead;

    double da = next->a - curr->a,
           db = next->b - curr->b;

    //
    // Loop over all subsequent vertices, up to curr
    //
    for (;;)
    {
      //
      // Get vertex after next
      //
      ABVertex* test = next->next;
      if (test == 0) test = vertexHead;

      //
      // If we are back to the original vertex, stop
      //
      if (test == curr) break;

      //
      // Test for parallel line segments
      //
      double dat = test->a - curr->a,
             dbt = test->b - curr->b;

      if (std::fabs(dat * db - dbt * da) > tolerance2) break;

      //
      // Redundant vertex found: do we have > 3 vertices?
      //
      if (numVertices <= 3)
      {
        CalculateMaxMin();
        return false;
      }

      //
      // Delete vertex pointed to by next. Carefully!
      //
      if (curr->next)
      {
        // next is not head
        if (next->next)
          curr->next = test;  // next is not tail
        else
          curr->next = 0;   // New tail
      }
      else
        vertexHead = test;  // New head

      if ((curr != next) && (next != test)) delete next;

      numVertices--;

      //
      // Replace next by the vertex we just tested,
      // and keep on going...
      //
      next = test;
      da = dat;
      db = dbt;
    }
    curr = curr->next;
  }

  //
  // In principle, this is not needed, but why not just play it safe?
  //
  CalculateMaxMin();

  return true;
}


//
// ReverseOrder
//
// Reverse the order of the vertices
//
void UReduciblePolygon::ReverseOrder()
{
  //
  // Loop over all vertices
  //
  ABVertex* prev = vertexHead;
  if (prev == 0) return;  // No vertices

  ABVertex* curr = prev->next;
  if (curr == 0) return;  // Just one vertex

  //
  // Our new tail
  //
  vertexHead->next = 0;

  for (;;)
  {
    //
    // Save pointer to next vertex (in original order)
    //
    ABVertex* save = curr->next;

    //
    // Replace it with a pointer to the previous one
    // (in original order)
    //
    curr->next = prev;

    //
    // Last vertex?
    //
    if (save == 0) break;

    //
    // Next vertex
    //
    prev = curr;
    curr = save;
  }

  //
  // Our new head
  //
  vertexHead = curr;
}

// StartWithZMin
//
// Starting alway with Zmin=bMin
// This method is used for GenericPolycone
//
void UReduciblePolygon::StartWithZMin()
{
  ABVertex* curr = vertexHead;
  double bcurr = curr->b;
  ABVertex* prev = curr;
  while (curr)
  {
    if (curr->b < bcurr)
    {
      bcurr = curr->b;
      ABVertex* curr1 = curr;
      while (curr1)
      {
        if (curr1->next == 0)
        {
          curr1->next = vertexHead;
          break;
        }
        curr1 = curr1->next;
      }
      vertexHead = curr;
      prev->next = 0;
    }
    prev = curr;
    curr = curr->next;
  }
}



//
// CrossesItself
//
// Return "true" if the polygon crosses itself
//
// Warning: this routine is not very fast (runs as N**2)
//
bool UReduciblePolygon::CrossesItself(double tolerance)
{
  double tolerance2 = tolerance * tolerance;
  double one  = 1.0 - tolerance,
         zero = tolerance;
  //
  // Top loop over line segments. By the time we finish
  // with the second to last segment, we're done.
  //
  ABVertex* curr1 = vertexHead, *next1 = 0;
  while (curr1->next)
  {
    next1 = curr1->next;
    double da1 = next1->a - curr1->a,
           db1 = next1->b - curr1->b;

    //
    // Inner loop over subsequent line segments
    //
    ABVertex* curr2 = next1->next;
    while (curr2)
    {
      ABVertex* next2 = curr2->next;
      if (next2 == 0) next2 = vertexHead;
      double da2 = next2->a - curr2->a,
             db2 = next2->b - curr2->b;
      double a12 = curr2->a - curr1->a,
             b12 = curr2->b - curr1->b;

      //
      // Calculate intersection of the two lines
      //
      double deter = da1 * db2 - db1 * da2;
      if (std::fabs(deter) > tolerance2)
      {
        double s1, s2;
        s1 = (a12 * db2 - b12 * da2) / deter;

        if (s1 >= zero && s1 < one)
        {
          s2 = -(da1 * b12 - db1 * a12) / deter;
          if (s2 >= zero && s2 < one) return true;
        }
      }
      curr2 = curr2->next;
    }
    curr1 = next1;
  }
  return false;
}



//
// BisectedBy
//
// Decide if a line through two points crosses the polygon, within tolerance
//
bool UReduciblePolygon::BisectedBy(double a1, double b1,
                                   double a2, double b2,
                                   double tolerance)
{
  int nNeg = 0, nPos = 0;

  double a12 = a2 - a1, b12 = b2 - b1;
  double len12 = std::sqrt(a12 * a12 + b12 * b12);
  a12 /= len12;
  b12 /= len12;

  ABVertex* curr = vertexHead;
  do
  {
    double av = curr->a - a1,
           bv = curr->b - b1;

    double Cross = av * b12 - bv * a12;

    if (Cross < -tolerance)
    {
      if (nPos) return true;
      nNeg++;
    }
    else if (Cross > tolerance)
    {
      if (nNeg) return true;
      nPos++;
    }
    curr = curr->next;
  }
  while (curr);

  return false;
}



//
// Area
//
// Calculated signed polygon area, where polygons specified in a
// clockwise manner (where x==a, y==b) have negative area
//
//    References: [O' Rourke (C)] pp. 18-27; [Gems II] pp. 5-6:
//    "The Area of a Simple Polygon", Jon Rokne.
//
double UReduciblePolygon::Area()
{
  double answer = 0;

  ABVertex* curr = vertexHead, *next;
  do
  {
    next = curr->next;
    if (next == 0) next = vertexHead;

    answer += curr->a * next->b - curr->b * next->a;
    curr = curr->next;
  }
  while (curr);

  return 0.5 * answer;
}


//
// Print
//
void UReduciblePolygon::Print()
{
  ABVertex* curr = vertexHead;
  do
  {
    std::cerr << curr->a << " " << curr->b << std::endl;
    curr = curr->next;
  }
  while (curr);
}


//
// CalculateMaxMin
//
// To be called when the vertices are changed, this
// routine re-calculates global values
//
void UReduciblePolygon::CalculateMaxMin()
{
  ABVertex* curr = vertexHead;
  aMin = aMax = curr->a;
  bMin = bMax = curr->b;
  curr = curr->next;
  while (curr)
  {
    if (curr->a < aMin)
      aMin = curr->a;
    else if (curr->a > aMax)
      aMax = curr->a;

    if (curr->b < bMin)
      bMin = curr->b;
    else if (curr->b > bMax)
      bMax = curr->b;

    curr = curr->next;
  }
}
