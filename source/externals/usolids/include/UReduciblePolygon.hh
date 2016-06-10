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
// Class description:
//
//   Utility class used to specify, test, reduce, and/or otherwise
//   manipulate a 2D polygon.
//
//   For this class, a polygon consists of n > 2 points in 2D
//   space (a,b). The polygon is always closed by connecting the
//   last point to the first. A UReduciblePolygon is guaranteed
//   to fulfill this definition in all instances.
//
//   Illegal manipulations (such that a valid polygon would be
//   produced) result in an error return if possible and
//   otherwise a // UException.
//
//   The Set of manipulations is limited currently to what
//   is needed for UPolycone and UPolyhedra.
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UReduciblePolygon_hh
#define UReduciblePolygon_hh

#include "UTypes.hh"

class UReduciblePolygon
{
    friend class UReduciblePolygonIterator;

  public:
    //
    // Creator: via simple a/b arrays
    //
    UReduciblePolygon(const double a[], const double b[], int n);

    //
    // Creator: a special version for UPolygon and UPolycone
    // that takes two a points at planes of b
    // (where a==r and b==z for the GEANT3 classic PCON and PGON)
    //
    UReduciblePolygon(const double rmin[], const double rmax[],
                      const double z[], int n);

    virtual ~UReduciblePolygon();

    //
    // Queries
    //
    inline int NumVertices() const
    {
      return numVertices;
    }

    inline double Amin() const
    {
      return aMin;
    }
    inline double Amax() const
    {
      return aMax;
    }
    inline double Bmin() const
    {
      return bMin;
    }
    inline double Bmax() const
    {
      return bMax;
    }

    void CopyVertices(double a[], double b[]) const;

    //
    // Manipulations
    //
    void ScaleA(double scale);
    void ScaleB(double scale);

    bool RemoveDuplicateVertices(double tolerance);
    bool RemoveRedundantVertices(double tolerance);

    void ReverseOrder();
    void StartWithZMin();
    //
    // Tests
    //
    double Area();
    bool CrossesItself(double tolerance);
    bool BisectedBy(double a1, double b1,
                    double a2, double b2, double tolerance);

    void Print(); // Debugging only

  public: // without description

  protected:

    void Create(const double a[], const double b[], int n);

    void CalculateMaxMin();

    //
    // Below are member values that are *always* kept up to date (please!)
    //
    double aMin, aMax, bMin, bMax;
    int  numVertices;

    //
    // A subclass which holds the vertices in a single-linked list
    //
    // Yeah, call me an old-fashioned c hacker, but I cannot make
    // myself use the rogue tools for this trivial list.
    //
    struct ABVertex;              // Secret recipe for allowing
    friend struct ABVertex;      // protected nested structures
    struct ABVertex
    {
      ABVertex() : a(0.), b(0.), next(0) {}
      double a, b;
      ABVertex* next;
    };

    ABVertex* vertexHead;

  private:

    UReduciblePolygon(const UReduciblePolygon&);
    UReduciblePolygon& operator=(const UReduciblePolygon&);
    // Private copy constructor and assignment operator.
};


//
// A companion class for iterating over the vertices of our polygon.
// It is simple enough that all routines are declared inline here.
//
class UReduciblePolygonIterator
{
  public:

    UReduciblePolygonIterator(const UReduciblePolygon* theSubject)
    {
      subject = theSubject;
      current = 0;
    }

    void  Begin()
    {
      current = subject->vertexHead;
    }
    bool  Next()
    {
      if (current) current = current->next;
      return Valid();
    }

    bool  Valid() const
    {
      return current != 0;
    }

    double GetA() const
    {
      return current->a;
    }
    double GetB() const
    {
      return current->b;
    }

  protected:

    const UReduciblePolygon* subject;     // Who are we iterating over
    UReduciblePolygon::ABVertex* current; // Current vertex
};

#endif
