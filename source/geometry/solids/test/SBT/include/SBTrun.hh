//
// SBTrun.hh
//
// Definition of the batch solid test
//

#ifndef SBTrun_hh
#define SBTrun_hh

#include "g4std/iostream"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"


class SBTVisManager;

//
// This is a list of points we keep track of. The definition below
// is a little sloppy.
//
class SBTrunPointList {
	public:
	SBTrunPointList( G4int size );
	~SBTrunPointList();
	
	void AddPoint( G4ThreeVector newPoint );
	inline G4ThreeVector operator[] (G4int i ) const { return pointList[i]; }
	
	inline G4int NumPoints() const { return numPoints; }
	
	protected:
	G4ThreeVector	*pointList;
	G4int		numPoints;
	G4int		maxPoints;
};



class SBTrun {

	public:
	SBTrun();
	~SBTrun();
	void SetDefaults();
	
	void RunTest( const G4VSolid *testVolume, G4std::ostream &logger );

	G4int DrawError( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex,
			 SBTVisManager *visManager ) const;
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
				  const SBTrunPointList *inside, const SBTrunPointList *outside, 
				  const G4ThreeVector point, G4std::ostream &logger );
	void	TestInsidePoint(  const G4VSolid *testVolume, G4int *nError,
				  const SBTrunPointList *inside, const G4ThreeVector point, G4std::ostream &logger );

	void	ReportError( G4int *nError, const G4ThreeVector p, 
			     const G4ThreeVector v, const G4String comment, G4std::ostream &logger );
	void 	ClearErrors();		
	G4int 	CountErrors() const;		
	
	G4int	GetLoggedPV( G4std::istream &logger, const G4int errorIndex,
			     G4ThreeVector &p, G4ThreeVector &v        ) const;
	
	protected:
	G4ThreeVector	target,
			widths,
			grids;
	G4int		maxPoints,
			maxErrors;
			
	
	typedef struct sSBTrunErrorList {
		G4String	message;
		G4int		nUsed;
		struct sSBTrunErrorList *next;
	} SBTrunErrorList;
	
	SBTrunErrorList *errorList;

};

#endif
