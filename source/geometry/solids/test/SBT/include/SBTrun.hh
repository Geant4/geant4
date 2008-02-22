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
//
// SBTrun.hh
//
// Definition of the batch solid test
//

#ifndef SBTrun_hh
#define SBTrun_hh

#include <iostream>
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
	
	void RunTest( const G4VSolid *testVolume, std::ostream &logger );

	G4int DrawError( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex,
			 SBTVisManager *visManager ) const;
	G4int DebugInside( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToInP( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToInPV( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToOutP( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;
	G4int DebugToOutPV( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;
	G4int DebugSurfNorm( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const;

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

  static G4String GetCurrentSolid(void) ; 
  static void SetCurrentSolid(G4String SolidName) ;

	protected:
	G4ThreeVector	GetRandomPoint() const;
	G4double	GaussianRandom(const G4double cutoff) const;
	
	void	TestOutsidePoint( const G4VSolid *testVolume, G4int *nError,
				  const SBTrunPointList *inside, const SBTrunPointList *outside, 
				  const G4ThreeVector point, std::ostream &logger );
	void	TestInsidePoint(  const G4VSolid *testVolume, G4int *nError,
				  const SBTrunPointList *inside, const G4ThreeVector point, std::ostream &logger );

	void	ReportError( G4int *nError, const G4ThreeVector p, 
			     const G4ThreeVector v, G4double distance,
			     const G4String comment, std::ostream &logger );
	void 	ClearErrors();		
	G4int 	CountErrors() const;		
	
	G4int	GetLoggedPV( std::istream &logger, const G4int errorIndex,
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

        private:
  static G4String CurrentSolid; /* Solid to be tested, used to visualize */
    

};

#endif
