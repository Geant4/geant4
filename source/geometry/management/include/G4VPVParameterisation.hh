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
// $Id: G4VPVParameterisation.hh,v 1.4 2001-07-11 09:59:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VPVParamterisation
//
// Class description:
//
// Parameterisation class, able to compute the transformation and
// (indirectly) the dimensions of parameterised volumes, given a
// replication number.

// History:
// 25.07.96 P.Kent        Initial stub version
// 20.09.96 V.Grichine    Modifications for G4Trap/Cons/Sphere
// 31.10.96 V.Grichine    Modifications for G4Torus/Para
// 17.02.98 J.Apostolakis Allowing the parameterisation of Solid type

#ifndef G4VPVPARAMETERISATION_HH
#define G4VPVPARAMETERISATION_HH

#include "globals.hh"

class G4VPhysicalVolume;

// CSG Entities which may be parameterised/replicated

class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4VSolid;
class G4Material;

class G4VPVParameterisation
{
  public:

    virtual void ComputeTransformation(const G4int,
                                       G4VPhysicalVolume *) const = 0;

    virtual G4VSolid*   ComputeSolid(const G4int, G4VPhysicalVolume *);
				       
    virtual G4Material* ComputeMaterial(const G4int, G4VPhysicalVolume *);
				       
    virtual void ComputeDimensions(G4Box &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

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

    virtual void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
};

#endif
