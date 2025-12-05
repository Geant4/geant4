//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4ReduciblePolygon
//
// Class description:
//
// Utility class used to specify, test, reduce, and/or otherwise
// manipulate a 2D polygon.
//
// For this class, a polygon consists of n > 2 points in 2D
// space (a,b). The polygon is always closed by connecting the
// last point to the first. A G4ReduciblePolygon is guaranteed
// to fulfill this definition in all instances. 
//
// Illegal manipulations (such that a valid polygon would be
// produced) result in an error return if possible and 
// otherwise a G4Exception.
//
// The set of manipulations is limited currently to what
// is needed for G4Polycone and G4Polyhedra.

// Author: David C. Williams (UCSC), 1998
// --------------------------------------------------------------------
#ifndef G4REDUCIBLEPOLYGON_HH
#define G4REDUCIBLEPOLYGON_HH

#include "G4Types.hh"

/**
 * @brief G4ReduciblePolygon is a utility class used to specify, test, reduce,
 * and/or otherwise manipulate a 2D polygon.
 */

class G4ReduciblePolygon
{
  friend class G4ReduciblePolygonIterator;

  public:

    /**
     * Constructor of G4ReduciblePolygon via simple a/b arrays.
     *  @param[in] a First array of points.
     *  @param[in] b Second array of points.
     *  @param[in] n The number of vertices of the polygon (has to be >=3).
     */
    G4ReduciblePolygon( const G4double a[], const G4double b[], G4int n );
  
    /**
     * Special constructor version for G4Polyhedra and G4Polycone, that takes
     * two a points at planes of b (where a==r and b==z).
     *  @param[in] rmin Array of r-min coordinates of corners.
     *  @param[in] rmax Array of r-max coordinates of corners.
     *  @param[in] z Array of Z coordinates of corners.
     *  @param[in] n The number of vertices of the polygon.
     */
    G4ReduciblePolygon( const G4double rmin[], const G4double rmax[],
                        const G4double z[], G4int n );

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4ReduciblePolygon(const G4ReduciblePolygon&) = delete;
    G4ReduciblePolygon& operator=(const G4ReduciblePolygon&) = delete;

    /**
     * Destructor, taking care to clear allocated lists.
     */
    ~G4ReduciblePolygon();
  
    /**
     * Accessors.
     */
    inline G4int NumVertices() const { return numVertices; }
    inline G4double Amin() const { return aMin; }
    inline G4double Amax() const { return aMax; }
    inline G4double Bmin() const { return bMin; }
    inline G4double Bmax() const { return bMax; }
  
    /**
     * Copies contents of provided arrays into simple linear arrays.
     */
    void CopyVertices( G4double a[], G4double b[] ) const;

    /**
     * Methods to multiply all a or b values by a common scale.
     */
    void ScaleA( G4double scale );
    void ScaleB( G4double scale );
  
    /**
     * Removes adjacent vertices that are equal.
     *  @param[in] tolerance Provided tolerance for adjacent vertices.
     *  @returns false, if there is a problem (too few vertices remaining).
     */
    G4bool RemoveDuplicateVertices( G4double tolerance );

    /**
     * Removes any unneeded vertices, i.e. those vertices which are on the
     * line connecting the previous and next vertices.
     *  @param[in] tolerance Provided tolerance for parallel line segments.
     *  @returns false, if there is a problem (too few vertices remaining).
     */
    G4bool RemoveRedundantVertices( G4double tolerance );
  
    /**
     * Reverses the order of the vertices.
     */
    void ReverseOrder();

    /**
     * Method is used for G4GenericPolycone; starting always with Zmin=bMin.
     */
    void StartWithZMin();

    // Methods for tests

    /**
     * Calculates signed polygon area, where polygons specified in a
     * clockwise manner have negative area.
     */
    G4double Area();

    /**
     * Returns "true" if the polygon crosses itself.
     */
    G4bool CrossesItself( G4double tolerance );

    /**
     * Decides if a line through two points crosses the polygon,
     * within tolerance.
     */
    G4bool BisectedBy( G4double a1, G4double b1,
                       G4double a2, G4double b2, G4double tolerance );
   
    /**
     * Print function for debugging.
     */
    void Print();
  
    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4ReduciblePolygon(__void__&);

  private:
  
    /**
     * Create the polygon; used in constructors.
     */
    void Create( const G4double a[], const G4double b[], G4int n );
  
    /**
     * Re-calculates global values. To be called when the vertices are changed.
     */
    void CalculateMaxMin();
  
  private:
  
    // Below are member values that are *always* kept up to date
    //
    G4double aMin, aMax, bMin, bMax;
    G4int numVertices = 0;
  

    /**
     * A subclass which holds the vertices in a single-linked list.
     */
    struct ABVertex;              // Secret recipe for allowing
    friend struct ABVertex;       // protected nested structures
    struct ABVertex
    {
      ABVertex()  = default;
      G4double a{0.}, b{0.};
      ABVertex *next{nullptr};
    };
  
    ABVertex* vertexHead = nullptr;
};

// A companion class for iterating over the vertices of our polygon.
// It is simple enough that all routines are declared inline here.
//

/**
 * @brief G4ReduciblePolygonIterator is companion class for iterating over
 * the vertices of a polygon.
 */

class G4ReduciblePolygonIterator
{
  public:

    inline G4ReduciblePolygonIterator( const G4ReduciblePolygon* theSubject )
    {
      subject = theSubject; current = nullptr;
    }
  
    inline void  Begin() { current = subject->vertexHead; }  

    inline G4bool  Next()
    {
      if (current != nullptr) { current=current->next; }
      return Valid();
    }
  
    inline G4bool  Valid() const { return current != nullptr; }  
  
    inline G4double GetA() const { return current->a; }
    inline G4double GetB() const { return current->b; }
  
  private:

    const G4ReduciblePolygon* subject = nullptr;  // Who are we iterating over
    G4ReduciblePolygon::ABVertex* current = nullptr;  // Current vertex
};

#endif
