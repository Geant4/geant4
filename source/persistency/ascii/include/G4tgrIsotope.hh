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
// $Id: G4tgrIsotope.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrIsotope
//
// Class description:
//
// Transient class of a chemical element.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrIsotope_h
#define G4tgrIsotope_h

#include "globals.hh"

#include <vector>

class G4tgrIsotope
{ 
  public:  // with description

    G4tgrIsotope();
   ~G4tgrIsotope();

    G4tgrIsotope( const std::vector<G4String>& wl );
      // Construct the G4tgrIsotope (fill its data members)
      // interpreting the data in the list of words 'wl' 
  
    // Retrieval methods

    const G4String& GetName() const { return theName; }
    G4int    GetZ() const { return theZ; }
    G4int    GetN() const { return theN; }
    G4double GetA() const { return theA; }

    friend std::ostream& operator<<(std::ostream& os, const G4tgrIsotope& obj);
  
 private:

    G4String theName;              // name of the Isotope
    G4int    theZ;                 // atomic number
    G4int    theN;                 // number of nucleons
    G4double theA;                 // mass
};

#endif
