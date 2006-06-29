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
// $Id: G4ThreeMat.hh,v 1.10 2006-06-29 18:40:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ThreeMat
//
// Class description:
// 
// Defines the class G4ThreeMat for three by three matrices.

// Author: Alan Breakstone
// Adaptation by: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __THREEMAT_H
#define __THREEMAT_H

#include "G4Vector3D.hh"

class G4ThreeMat
{

public:  // with description

  G4ThreeMat();
    // Default constructor.

  G4ThreeMat( G4double a[3][3] );
    // Normal constructors with a 3 x 3 array argument.

  virtual ~G4ThreeMat();
    // Destructor.

  G4ThreeMat( const G4ThreeMat& m );
  G4ThreeMat& operator=( const G4ThreeMat& m );
    // Copy constructor and assignment operator.

  G4int operator==( const G4ThreeMat& m ) const;
    // Equality operator.

  friend std::ostream& operator<<( std::ostream& os, const G4ThreeMat& m );

  virtual const char* NameOf() const;
    // Returns the class name.

  virtual void PrintOn( std::ostream& os = G4cout ) const;
    // Printing functions (derived classes do not need to overwrite operator <<).

  G4double Determinant() const;
    // Determinant of matrix.

private:

  G4double element[3][3];
  G4Vector3D row[3], column[3];
    // The elements exist individually and are also aggregated into rows and
    // columns to use operations already written for the G4Vector3Dc class.
};

#endif
