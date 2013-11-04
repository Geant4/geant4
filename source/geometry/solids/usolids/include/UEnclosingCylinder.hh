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
// $Id: UEnclosingCylinder.hh 66241 2012-12-13 18:34:42Z gunter $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// UEnclosingCylinder
//
// Class description:
//
//	 Definition of a utility class for quickly deciding if a point
//	 is clearly outside a polyhedra or polycone or deciding if
//	 a trajectory is clearly going to miss those shapes.

// Author: 
//	 David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef UEnclosingCylinder_hh
#define UEnclosingCylinder_hh

#include "UTypes.hh"
#include "UTubs.hh"

class UReduciblePolygon;

class UEnclosingCylinder
{
	public:	// with description

		UEnclosingCylinder( /*const UReduciblePolygon *rz*/  double r, double lo, double hi,
															 bool phiIsOpen, 
															 double startPhi, double totalPhi );
		~UEnclosingCylinder();
	
		bool MustBeOutside( const UVector3 &p ) const;
			// Decide very rapidly if the point is outside the cylinder.
			// If one is not certain, return false.

		bool ShouldMiss( const UVector3 &p, const UVector3 &v ) const;
			// Decide very rapidly if the trajectory is going to miss the cylinder.
			// If one is not sure, return false.

    double DistanceTo( const UVector3 &p, const UVector3 &v ) const;

    double SafetyFromOutside( const UVector3 &p) const;

	public:	// without description

    void Extent( UVector3 &aMin, UVector3 &aMax ) const;

    double radius;		// radius of our cylinder

	protected:

		double zLo, zHi;	// z extent
	
		bool		phiIsOpen; // true if there is a phi segment
		double	startPhi,	// for isPhiOpen==true, starting of phi segment
							totalPhi;	// for isPhiOpen==true, size of phi segment

		double rx1, ry1,
						 dx1, dy1;
		double rx2, ry2,
						 dx2, dy2;
		 
		bool	 concave;	// true, if x/y Cross section is concave
		
    UTubs *tube;
};

#endif
