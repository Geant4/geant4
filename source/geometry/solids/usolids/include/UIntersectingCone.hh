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
// $Id: UIntersectingCone.hh,v 1.13 2010-07-12 15:33:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// UIntersectingCone
//
// Class description:
//
//	 Utility class which calculates the intersection
//	 of an arbitrary line with a fixed cone

// Author: 
//	 David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef UIntersectingCone_hh
#define UIntersectingCone_hh

#include "UTypes.hh"
#include "geomdefs.hh"


class UIntersectingCone
{
	public:

		UIntersectingCone( const double r[2], const double z[2] );
		virtual ~UIntersectingCone();
	
		int LineHitsCone( const UVector3 &p, const UVector3 &v, double &s1, double &s2);
	
		bool HitOn( const double r, const double z );
	
		inline double RLo() const { return rLo; }
		inline double RHi() const { return rHi; }
		inline double ZLo() const { return zLo; }
		inline double ZHi() const { return zHi; }
	
	public:	// without description

		/*
		UIntersectingCone(__void__&);
			// Fake default constructor for usage restricted to direct object
			// persistency for clients requiring preallocation of memory for
			// persistifiable objects.
		*/


	protected:

		double zLo, zHi,	// Z bounds of side
						 rLo, rHi;	// R bounds of side

		bool	 type1;		// True if cone is type 1
											 //	(std::fabs(z1-z2)>std::fabs(r1-r2))
		double A, B;		 // Cone radius parameter:
											 //	type 1: r = A + B*z
											 //	type 2: z = A + B*r

		// int Solution (const UVector3 &p, const UVector3 &v, double a, double b, double c, double &s1, double &s2);

		int LineHitsCone1( const UVector3 &p, const UVector3 &v,
                           double &s1, double &s2 );

		int LineHitsCone1Optimized( const UVector3 &p, const UVector3 &v,
			double &s1, double &s2 );

		int LineHitsCone2( const UVector3 &p, const UVector3 &v,
                           double &s1, double &s2 );

//		const double kInfinity;
		const static double EpsilonQuad;
};

#endif
