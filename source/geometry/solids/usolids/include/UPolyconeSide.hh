//
// ********************************************************************
// * License and Disclaimer																					 *
// *																																	*
// * The	Geant4 software	is	copyright of the Copyright Holders	of *
// * the Geant4 Collaboration.	It is provided	under	the terms	and *
// * conditions of the Geant4 Software License,	included in the file *
// * LICENSE and available at	http://cern.ch/geant4/license .	These *
// * include a list of copyright holders.														 *
// *																																	*
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work	make	any representation or	warranty, express or implied, *
// * regarding	this	software system or assume any liability for its *
// * use.	Please see the license in the file	LICENSE	and URL above *
// * for the full disclaimer and the limitation of liability.				 *
// *																																	*
// * This	code	implementation is the result of	the	scientific and *
// * technical work of the GEANT4 collaboration.											*
// * By using,	copying,	modifying or	distributing the software (or *
// * any work based	on the software)	you	agree	to acknowledge its *
// * use	in	resulting	scientific	publications,	and indicate your *
// * acceptance of all terms of the Geant4 Software license.					*
// ********************************************************************
//
//
// $Id: UPolyconeSide.hh 66885 2013-01-16 17:37:13Z gunter $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// UPolyconeSide
//
// Class description:
//
//	 Class implmenting a face that represents one conical side
//	 of a polycone:
//
//	 UPolyconeSide( const UPolyconeSideRZ *prevRZ,
//									 const UPolyconeSideRZ *tail,
//									 const UPolyconeSideRZ *head,
//									 const UPolyconeSideRZ *nextRZ,
//												 double phiStart, double deltaPhi, 
//												 bool phiIsOpen, bool isAllBehind=false )
//
//	 Values for r1,z1 and r2,z2 should be specified in clockwise
//	 order in (r,z).

// Author: 
//	 David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef UPolyconeSide_hh
#define UPolyconeSide_hh

#include "UVCSGface.hh"

class UIntersectingCone;

struct UPolyconeSideRZ
{
	double r, z;	// start of vector
};

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The class PolyconeSidePrivateSubclass is introduced to
//encapsulate the fields of the class UPolyconeSide that may not
//be read-only.
#ifndef POLYCONESIDEPRIVATESUBCLASS_HH
#define POLYCONESIDEPRIVATESUBCLASS_HH

class PolyconeSidePrivateSubclass
{
public:
	std::pair<UVector3, double> fPhi;	// Cached value for phi

	void initialize() {
		fPhi.first = UVector3(0,0,0);
		fPhi.second= 0.0;
	};
};
#endif

//01.25.2009 Xin Dong: Phase II change for Geant4 multithreading.
//The class UPolyconeSideSubInstanceManager is introduced to 
//encapsulate the methods used by both the master thread and 
//worker threads to allocate memory space for the fields encapsulated
//by the class PolyconeSidePrivateSubclass. When each thread
//initializes the value for these fields, it refers to them using a macro
//definition defined below. For every UPolyconeSide instance, there is
//a corresponding PolyconeSidePrivateSubclass instance. All
//PolyconeSidePrivateSubclass instances are organized by the
//class UPolyconeSideSubInstanceManager as an array. The field "	
//int g4polyconeSideSubInstanceID" is added to the class UPolyconeSide.
//The value of this field in each UPolyconeSide instance is the subscript
//of the corresponding PolyconeSidePrivateSubclass instance. In order
//to use the class UPolyconeSideSubInstanceManager, we add a static member in
//the class UPolyconeSide as follows: "	
//static UPolyconeSideSubInstanceManager g4polyconeSideSubInstanceManager".
//For the master thread, the array for PolyconeSidePrivateSubclass 
//instances grows dynamically along with UPolyconeSide instances are
//created. For each worker thread, it copies the array of 
//PolyconeSidePrivateSubclass instances from the master thread.
//In addition, it invokes a method similiar to the constructor explicitly
//to achieve the partial effect for each instance in the array.

#ifndef UPolyconeSideSUBINSTANCEMANAGER_HH
#define UPolyconeSideSUBINSTANCEMANAGER_HH

//#include "UMTTransitory.hh"
//typedef UMTPrivateSubInstanceManager<PolyconeSidePrivateSubclass> UPolyconeSideSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//These macros changes the references to fields that are now encapsulated
//in the class PolyconeSidePrivateSubclass.
#define fPhiPCSUMTThreadPrivate ((g4polyconeSideSubInstanceManager.offset[g4polyconeSideInstanceID]).fPhi)

#endif

class UPolyconeSide : public UVCSGface
{
	public:

		//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
		//This new field is used as instance ID.
		int g4polyconeSideInstanceID;

		//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
		//This new field helps to use the class UPolyconeSideSubInstanceManager
		//introduced above.
		
    //static UPolyconeSideSubInstanceManager g4polyconeSideSubInstanceManager;


		UPolyconeSide( const UPolyconeSideRZ *prevRZ,
										const UPolyconeSideRZ *tail,
										const UPolyconeSideRZ *head,
										const UPolyconeSideRZ *nextRZ,
													double phiStart, double deltaPhi, 
													bool phiIsOpen, bool isAllBehind=false );
		virtual ~UPolyconeSide();
	
		UPolyconeSide( const UPolyconeSide &source );
		UPolyconeSide& operator=( const UPolyconeSide &source );
	
		bool Distance( const UVector3 &p, const UVector3 &v,	
														bool outgoing, double surfTolerance,
														double &distance, double &distFromSurface,
														UVector3 &normal, bool &isAllBehind );

		double Safety( const UVector3 &p, bool outgoing );
	
		VUSolid::EnumInside Inside( const UVector3 &p, double tolerance, 
													double *bestDistance );
	
		UVector3 Normal( const UVector3 &p,	double *bestDistance );

		double Extent( const UVector3 axis );

    /*
		void CalculateExtent( const EAxisType axis, 
													const UVoxelLimits &voxelLimit,
													const UAffineTransform &tranform,
																USolidExtentList &extentList			 );
                                */

		UVCSGface *Clone() { return new UPolyconeSide( *this ); }

		double SurfaceArea();
		UVector3 GetPointOnFace();
	
	public:	// without description

		UPolyconeSide(__void__&);
			// Fake default constructor for usage restricted to direct object
			// persistency for clients requiring preallocation of memory for
			// persistifiable objects.

	protected:

		double DistanceAway( const UVector3 &p, bool opposite,
																 double &distOutside2, double *rzNorm=0 );
			
		bool PointOnCone( const UVector3 &hit, double normSign,
												const UVector3 &p,
												const UVector3 &v, UVector3 &normal );

		void CopyStuff( const UPolyconeSide &source );
	
		static void FindLineIntersect( double x1, double y1,
																	 double tx1, double ty1,
																	 double x2, double y2,
																 double tx2, double ty2,
																 double &x, double &y );

		double GetPhi( const UVector3& p );

	protected:

		double r[2], z[2]; // r, z parameters, in specified order
		double startPhi,	 // Start phi (0 to 2pi), if phiIsOpen
						 deltaPhi;	 // Delta phi (0 to 2pi), if phiIsOpen
		bool	 phiIsOpen;	// True if there is a phi slice
		bool	 allBehind;	// True if the entire solid is "behind" this face
	
		UIntersectingCone *cone;	// Our intersecting utility class
	
		double rNorm, zNorm;	// Normal to surface in r,z space
		double rS, zS;				// Unit vector along surface in r,z space
		double length;				// Length of face in r,z space
		double prevRS,
						 prevZS;				// Unit vector along previous polyconeSide
		double nextRS,
						 nextZS;				// Unit vector along next polyconeSide
	
		double rNormEdge[2],
						 zNormEdge[2];	// Normal to edges

		int ncorners;
		UVector3 *corners; // The coordinates of the corners (if phiIsOpen)

	private:
		//Change to thread private
		//		std::pair<UVector3, double> fPhiPCSUMTThreadPrivate;	// Cached value for phi
		double tolerance; // Geometrical surface thickness
		double fSurfaceArea;	// Used for surface calculation 
};

#endif
