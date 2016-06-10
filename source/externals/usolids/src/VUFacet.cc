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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: VUFacet.cc,v 1.11 2010-09-23 10:30:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
// Author: Marek Gayer, started from original implementation by P R Truscott, 2004
//

#include "VUFacet.hh"
#include "VUSolid.hh"

using namespace std;

bool VUFacet::operator== (const VUFacet &right) const
{
	double tolerance = kCarTolerance*kCarTolerance/4.0;

	if (GetNumberOfVertices() != right.GetNumberOfVertices())
		return false;
	else if ((GetCircumcentre()-right.GetCircumcentre()).Mag2() > tolerance)
		return false;
	else if (std::fabs((right.GetSurfaceNormal()).Dot(GetSurfaceNormal())) < 0.9999999999)
		return false;

	bool coincident  = true;
	int i = 0;
	do
	{
		coincident = false;
		int j   = 0; 
		do
		{
			coincident = (GetVertex(i)-right.GetVertex(j)).Mag2() < tolerance;
		} while (!coincident && ++j < GetNumberOfVertices());
	} while (coincident && ++i < GetNumberOfVertices());

	return coincident;
}

///////////////////////////////////////////////////////////////////////////////
//
void VUFacet::ApplyTranslation(const UVector3 v)
{
	int n = GetNumberOfVertices();
	for (int i = 0; i < n; ++i)
		SetVertex(i, GetVertex(i) + v);
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream &VUFacet::StreamInfo(std::ostream &os) const
{
	os << endl;
	os << "*********************************************************************" << endl;
	os << "FACET TYPE       = " << GetEntityType() << endl;
	os << "ABSOLUTE VECTORS = " << endl;
	int n = GetNumberOfVertices();
	for (int i = 0; i < n; ++i)
		os << "P[" << i << "]      = " << GetVertex(i) << endl;

	/*
	os << "RELATIVE VECTORS = " << endl;
	for (vector<UVector3>::const_iterator it=E.begin(); it!=E.end(); it++)
	{ os << "E[" << it-E.begin()+1 << "]      = " << *it << endl; }
	*/

	os << "*********************************************************************" << endl;

	return os;
}

bool VUFacet::IsInside (const UVector3 &p) const
{
	UVector3 d =  p - GetVertex(0);
	double displacement = d.Dot(GetSurfaceNormal());
	return displacement <= 0.0;
}

const double VUFacet::dirTolerance = 1.0E-14;
const double VUFacet::kCarTolerance = VUSolid::Tolerance();
