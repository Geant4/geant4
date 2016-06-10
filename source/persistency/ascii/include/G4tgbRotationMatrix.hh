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
// $Id: G4tgbRotationMatrix.hh 68052 2013-03-13 14:38:53Z gcosmo $
//
//
// class G4tgbRotationMatrix
//
// Class description:
//
// Transient class of a rotation matrix; builds a G4RotationMatrix,
// of each rotation matrix.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbRotationMatrix_h
#define G4tgbRotationMatrix_h

#include "globals.hh"

#include <vector>
#include <string>

#include "G4tgrRotationMatrix.hh"
#include "G4RotationMatrix.hh"

class G4tgbRotationMatrix
{
  public:  // with description

    G4tgbRotationMatrix();
   ~G4tgbRotationMatrix();

    G4tgbRotationMatrix( G4tgrRotationMatrix* tgr );
      // Construct the G4tgbRotationMatrix (fill its data members)
      // interpreting the data in the list of words 'wl'

    G4RotationMatrix* BuildG4RotMatrix( );
    G4RotationMatrix* BuildG4RotMatrixFrom3( std::vector<G4double>& values );
    G4RotationMatrix* BuildG4RotMatrixFrom6( std::vector<G4double>& values );
    G4RotationMatrix* BuildG4RotMatrixFrom9( std::vector<G4double>& values );
      // Build a G4RotationMatrix transforming theValues

    G4String GetName() { return theTgrRM->GetName(); }

  private:

    G4tgrRotationMatrix* theTgrRM;
};

#endif
 
