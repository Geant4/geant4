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
// $Id: UFacet.hh,v 1.8 2010-09-23 10:27:25 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Marek Gayer, started from original implementation by P R Truscott, 2004
//
//
// Class description:
//
//   Base class defining the facets which are components of a
//   UTessellatedSolid shape.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef UFacet_hh
#define UFacet_hh

#include <iostream>
#include <vector>

#include "UVector3.hh"
#include "UTypes.hh"

enum UFacetVertexType {UABSOLUTE, URELATIVE};

class UTessellatedSolid;

class VUFacet
{
public:

	virtual ~VUFacet () {};

	virtual int GetNumberOfVertices () const = 0;
	virtual UVector3 GetVertex (int i) const = 0;
	virtual void SetVertex (int i, const UVector3 &val) = 0;
	virtual UGeometryType GetEntityType () const = 0;
	virtual UVector3 GetSurfaceNormal () const = 0;
	virtual bool IsDefined () const = 0;
	virtual UVector3 GetCircumcentre () const = 0;
	virtual double GetRadius () const = 0;
	virtual VUFacet *GetClone () = 0;
	virtual double Distance (const UVector3&, const double) = 0;
	virtual double Distance (const UVector3&, const double, const bool) = 0;
	virtual double Extent (const UVector3) = 0;
	virtual bool Intersect (const UVector3&, const UVector3 &, const bool , double &, double &, UVector3 &) = 0;
	virtual double GetArea() = 0;
	virtual UVector3 GetPointOnFace() const = 0;

	bool operator== (const VUFacet &right) const;
	void ApplyTranslation (const UVector3 v);
	std::ostream &StreamInfo(std::ostream &os) const;
	bool IsInside(const UVector3 &p) const;

	virtual int AllocatedMemory() = 0;
	virtual void SetVertexIndex (const int i, const int j) = 0;
	virtual int GetVertexIndex (const int i) const = 0;

	virtual void SetVertices(std::vector<UVector3> *vertices) = 0;

protected:

	static const double dirTolerance;
	static const double kCarTolerance;
};

#endif
