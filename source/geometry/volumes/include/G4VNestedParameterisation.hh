//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VNestedParameterisation.hh,v 1.1 2005-02-23 10:53:53 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VNestedParamterisation
//
// Class description:
//
// Parameterisation class, that can use parent volume information 
// eg to compute the transformation and the dimensions of parameterised 
// volumes, given a replication number and parent volume information.

// History:
// 04 Feb 05 J.Apostolakis Re-enabling the parameterisation using parent info
// --------------------------------------------------------------------
#ifndef G4VNESTEDPARAMETERISATION_HH
#define G4VNESTEDPARAMETERISATION_HH

#include "G4Types.hh"
#include "G4VPVParameterisation.hh" 

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

// CSG Entities which may be parameterised/replicated
//
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

class G4VNestedParameterisation: public G4VPVParameterisation
{
  public:
    // Key method, required
    virtual void ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume *currentPV, 
				       const G4VTouchable* parentTouch) const = 0;

    // Optional methods, possibly to be overridden by concrete parameterisation
    virtual G4VSolid*   ComputeSolid(const G4int, 
				     G4VPhysicalVolume  *pvol,  // currentVol
				     const G4VTouchable *);     // parentTouch
    virtual G4Material* ComputeMaterial(const G4int, 
					G4VPhysicalVolume *,     // currentVol
					const G4VTouchable *);   // parentTouch

    // Implements the following methods, in terms of the above 
    //  providing 're-interpretation' to add information of parent
    virtual void ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume *currPhysTouch ) const;

    virtual G4VSolid*   ComputeSolid(const G4int no, G4VPhysicalVolume *thisVol);

    virtual G4Material* ComputeMaterial(const G4int, G4VPhysicalVolume *);

    virtual void ComputeDimensions(G4Box &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const;

    virtual void ComputeDimensions(G4Box &,
                                   const G4int,
                                   const G4VPhysicalVolume *,
				   const G4VTouchable *     ) const {}

    // Similar 're-interpretation' could be done 
    //     of the ComputeDimensions methods below
				       
    virtual void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

  private: 
    void ReportErrorInTouchable(const G4String method) const; 

    G4VPhysicalVolume*  ObtainParts( G4VPhysicalVolume* thisVolPlus,  // PhysicalTouchable
				     const G4String method, 
				     const G4VTouchable* pTouchableParent) const;
    const G4VPhysicalVolume*  ObtainParts( const G4VPhysicalVolume* thisVolPlus,  // PhysicalTouchable
				     const G4String method, 
				     const G4VTouchable* pTouchableParent) const;
       //  Demux the expected PhysicalTouchable into 
       //    - pointer to this physical volume
       //    - pointer to parent touchable
};

#endif
