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
// $Id: G4tgrPlaceDivRep.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrPlaceDivRep
//
// Class description:
//
// Class to descripe the position of a G4tgrVolume inside another
// G4tgrVolume as a replica, i.e. filling the whole parent volume.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrPlaceDivRep_h
#define G4tgrPlaceDivRep_h

#include "globals.hh"
#include "geomdefs.hh"

#include <vector>

#include "G4tgrPlace.hh"

enum G4DivType { DivByNdiv, DivByWidth, DivByNdivAndWidth };

class G4tgrPlaceDivRep : public G4tgrPlace
{

  public:  // with description

    G4tgrPlaceDivRep();
   ~G4tgrPlaceDivRep();

    G4tgrPlaceDivRep( const std::vector<G4String>& wl );
      // Creates an object passing the only data that is fixed
      // (ndiv, width, offset may be have to be recalculated)

    EAxis BuildAxis( const G4String& axisName );

    // Access functions

    EAxis GetAxis() const { return theAxis; }
    G4int GetNDiv() const { return theNDiv; }
    G4double GetWidth() const { return theWidth; }
    G4double GetOffset() const { return theOffset; }
    G4DivType GetDivType() const { return theDivType; }

    void SetParentName(const G4String& parentName) { theParentName=parentName; }
    void SetNDiv( G4int ndiv ) { theNDiv = ndiv; }
    void SetWidth( G4double width ) { theWidth = width; }
    void SetAxis( EAxis axis ) { theAxis = axis; }
    void SetOffset( G4double offset ) { theOffset = offset; }
    void SetDivType( G4DivType typ ) { theDivType = typ; }
    
    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrPlaceDivRep& obj);
  private:

    G4int theNDiv;
    G4double theWidth;
    EAxis theAxis;
    G4double theOffset;
    G4DivType theDivType;  
};

#endif
