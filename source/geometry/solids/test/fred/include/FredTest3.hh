//
// FredTest3.hh
//
// Definition of fred's 3rd test suite
//

#ifndef FredTest3_hh
#define FredTest3_hh

#include "g4std/iostream"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"


//
// This is a list of points we keep track of. The definition below
// is a little sloppy.
//
class FredTest3PointList {
	public:
	FredTest3PointList( G4int size );
	~FredTest3PointList();
	
	void AddPoint( G4ThreeVector newPoint );
	inline G4ThreeVector operator[] (G4int i ) const { return pointList[i]; }
	
	inline G4int NumPoints() const { return numPoints; }
	
	protected:
	G4ThreeVector	*pointList;
	G4int		numPoints;
	G4int		maxPoints;
};



class FredTest3 {

	public:
	FredTest3();
	~FredTest3();
	void SetDefaults();
	
	void RunTest( const G4VSolid *testVolume, G4std::ostream &logger );

	G4int DebugError( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;
	G4int DebugInside( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToInP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToInPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToOutP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToOutPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const;

	inline void SetTarget( const G4ThreeVector &newTarget ) { target = newTarget; }
	inline G4ThreeVector GetTarget() const { return target; }

	inline void SetWidths( const G4ThreeVector &newWidths ) { widths = newWidths; }
	inline G4ThreeVector GetWidths() const { return widths; }
	
	inline void SetGrids( const G4ThreeVector &newGrids ) { grids = newGrids; }
	inline G4ThreeVector GetGrids() const { return grids; }
	
	inline void SetMaxPoints( const G4int newMaxPoints ) { maxPoints = newMaxPoints; }
	inline G4int GetMaxPoints() const { return maxPoints; }
	
	inline void SetMaxErrors( const G4int newMaxErrors ) { maxErrors = newMaxErrors; }
	inline G4int GetMaxErrors() const { return maxErrors; }

	protected:
	G4ThreeVector	GetRandomPoint() const;
	G4double	GaussianRandom(const G4double cutoff) const;
	
	void	TestOutsidePoint( const G4VSolid *testVolume, G4int *nError,
				  const FredTest3PointList *inside, const G4ThreeVector point, G4std::ostream &logger );
	void	TestInsidePoint(  const G4VSolid *testVolume, G4int *nError,
				  const FredTest3PointList *inside, const G4ThreeVector point, G4std::ostream &logger );

	void	ReportError( G4int *nError, const G4ThreeVector p, 
			     const G4ThreeVector v, const G4String comment, G4std::ostream &logger );
	void 	ClearErrors();		
	
	G4int	GetLoggedPV( G4std::istream &logger, const G4int errorIndex,
			     G4ThreeVector &p, G4ThreeVector &v        ) const;
	
	protected:
	G4ThreeVector	target,
			widths,
			grids;
	G4int		maxPoints,
			maxErrors;
			
	
	typedef struct sFredTest3ErrorList {
		G4String	message;
		G4int		nUsed;
		struct sFredTest3ErrorList *next;
	} FredTest3ErrorList;
	
	FredTest3ErrorList *errorList;

};

#endif
