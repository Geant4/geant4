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
// $Id: ApproxEqual.hh,v 1.4 2006-06-29 18:49:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ApproxEqual Functions for geometry test programs
//
// History:
// 20.07.95 P.Kent Translated from old code

#ifndef APPROXEQUAL_HH
#define APPROXEQUAL_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"

const G4double kApproxEqualTolerance = 1E-6;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is check if less than kApproxEqualTolerance from target

G4bool ApproxEqual(const G4double check,const G4double target)
{
    return (std::fabs(check-target)<kApproxEqualTolerance) ?true:false;
}

// Return true if the 3vector check is approximately equal to target
G4bool ApproxEqual(const G4ThreeVector& check, const G4ThreeVector& target)
{
    return (ApproxEqual(check.x(),target.x())&&
	   ApproxEqual(check.y(),target.y())&&
	    ApproxEqual(check.z(),target.z())) ? true : false;
}


G4bool ApproxEqual(const G4AffineTransform &tf1,
		   const G4AffineTransform &tf2)
{
	for (G4int i=0;i<15;i++)
		{
		if (!ApproxEqual(tf1[i],tf2[i])) return false;
		}
	return true;
}

#endif









