// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReduciblePolygon.hh,v 1.2 2000-09-12 07:34:17 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4ReduciblePolygon.hh
//
// Class description:
//
//   Utility class used to specify, test, reduce, and/or otherwise
//   manipulate a 2D polygon.
//
//   For this class, a polygon consists of n > 2 points in 2D
//   space (a,b). The polygon is always closed by connecting the
//   last point to the first. A G4ReduciblePolygon is guaranteed
//   to fulfill this definition in all instances. 
//
//   Illegal manipulations (such that a valid polygon would be
//   produced) result in an error return if possible and 
//   otherwise a G4Exception.
//
//   The set of manipulations is limited currently to what
//   is needed for G4Polycone and G4Polyhedra.

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4ReduciblePolygon_hh
#define G4ReduciblePolygon_hh

#include "globals.hh"

class G4ReduciblePolygon {
	friend class G4ReduciblePolygonIterator;

	public:
	//
	// Creator: via simple a/b arrays
	//
	G4ReduciblePolygon( const G4double a[], const G4double b[], const G4int n );
	
	//
	// Creator: a special version for G4Polygon and G4Polycone
	// that takes two a points at planes of b
	// (where a==r and b==z for the GEANT3 classic PCON and PGON)
	//
	G4ReduciblePolygon( const G4double rmin[], const G4double rmax[], const G4double z[], const G4int n );
	
	virtual ~G4ReduciblePolygon();
	
	//
	// Queries
	//
	inline G4int NumVertices() const { return numVertices; }
	
	inline G4double Amin() const { return aMin; }
	inline G4double Amax() const { return aMax; }
	inline G4double Bmin() const { return bMin; }
	inline G4double Bmax() const { return bMax; }
	
	void CopyVertices( G4double a[], G4double b[] ) const;

	//
	// Manipulations
	//
	void ScaleA( const G4double scale );
	void ScaleB( const G4double scale );
	
	G4bool RemoveDuplicateVertices( const G4double tolerance );
	G4bool RemoveRedundantVertices( const G4double tolerance );
	
	void ReverseOrder();

	//
	// Tests
	//
	G4double Area();
	G4bool CrossesItself( const G4double tolerance );
	G4bool BisectedBy( const G4double a1, const G4double b1,
			   const G4double a2, const G4double b2, const G4double tolerance );
	
	void Print();	// Debugging only
	
	protected:
	
	void Create( const G4double a[], const G4double b[], const G4int n );
	
	void CalculateMaxMin();
	
	//
	// Below are member values that are *always* kept up to date (please!)
	//
	G4double aMin, aMax, bMin, bMax;
	G4int	 numVertices;
	
	//
	// A subclass which holds the vertices in a single-linked list
	//
	// Yeah, call me an old-fashioned c hacker, but I cannot make
	// myself use the rogue tools for this trivial list.
	//
	struct ABVertex;
	friend struct ABVertex;
	struct ABVertex {
		G4double a, b;
		ABVertex *next;
	};
	
	ABVertex *vertexHead;
};


//
// A companion class for iterating over the vertices of our polygon.
// It is simple enough that all routines are declared inline here.
//
class G4ReduciblePolygonIterator {
	public:
	G4ReduciblePolygonIterator( const G4ReduciblePolygon *theSubject ) { subject = theSubject; current=0; }
	
	inline void	Begin() { current = subject->vertexHead; }	
	inline G4bool	Next()  { if (current) current = current->next; return Valid(); }		
	
	inline G4bool	Valid() const { return current!=0; }	
	
	inline G4double GetA() const { return current->a; }
	inline G4double GetB() const { return current->b; }
	
	protected:
	const G4ReduciblePolygon	*subject;	// Who are we iterating over
	G4ReduciblePolygon::ABVertex	*current;	// Current vertex
};

#endif
