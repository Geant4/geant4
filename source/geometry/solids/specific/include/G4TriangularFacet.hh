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
// $Id: G4TriangularFacet.hh,v 1.8 2007-12-10 16:30:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TriangularFacet.hh
//
// Date:                15/06/2005
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:            C/MAT/N03517
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class description:
//
//   The G4TriangularFacet class is used for the contruction of
//   G4TessellatedSolid.
//   It is defined by three vertices, which shall be supplied in anti-clockwise
//   order looking from the outsider of the solid where it belongs.
//   Its constructor:
//   
//      G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
//                         const G4ThreeVector vt2, G4FacetVertexType);
//
//   takes 4 parameters to define the three vertices:
//      1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1 and vt2 are 
//         the 3 vertices in anti-clockwise order looking from the outsider.
//      2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//         the second vertex is Pt0+vt1 and the third vertex is Pt0+vt2, all  
//         in anti-clockwise order when looking from the outsider.

///////////////////////////////////////////////////////////////////////////////
#ifndef G4TriangularFacet_hh
#define G4TriangularFacet_hh 1

#include "globals.hh"
#include "G4VFacet.hh"
#include "G4ThreeVector.hh"
#include "G4TessellatedGeometryAlgorithms.hh"

class G4TriangularFacet : public G4VFacet
{
  public:  // with description

    G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
                       const G4ThreeVector vt2, G4FacetVertexType);
    ~G4TriangularFacet ();
    
    G4TriangularFacet (const G4TriangularFacet &right);
    const G4TriangularFacet &operator=(G4TriangularFacet &right);    

    G4VFacet *GetClone ();
    G4TriangularFacet *GetFlippedFacet ();
    
    G4ThreeVector Distance (const G4ThreeVector &p);
    G4double Distance (const G4ThreeVector &p, const G4double minDist);
    G4double Distance (const G4ThreeVector &p, const G4double minDist,
                       const G4bool outgoing);
    G4double Extent   (const G4ThreeVector axis);
    G4bool Intersect  (const G4ThreeVector &p, const G4ThreeVector &v,
                       const G4bool outgoing, G4double &distance,
                             G4double &distFromSurface, G4ThreeVector &normal);
    G4double GetArea ();
    G4ThreeVector GetPointOnFace () const;

  private:

    G4double a;
    G4double b;
    G4double c;
    G4double det;
    
    G4double sMin, sMax;
    G4double tMin;

    G4double sqrDist;

    G4TessellatedGeometryAlgorithms *tGeomAlg;
};

#endif
