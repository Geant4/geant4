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
// $Id: G4tgrSolid.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrSolid
//
// Class description:
//
// Base class for management of transient solids.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrSolid_h
#define G4tgrSolid_h

#include "globals.hh"

#include <vector>

#include "G4ThreeVector.hh"

class G4tgrSolid
{
  public:  // with description

    G4tgrSolid();
    G4tgrSolid( const std::vector<G4String>& wl);
    virtual ~G4tgrSolid();

    friend std::ostream& operator<<(std::ostream&, const G4tgrSolid&);

    // Accessors  

    const G4String& GetName() const { return theName; }
    const G4String& GetType() const { return theType; }
    const std::vector< std::vector<G4double>* > GetSolidParams() const;
    virtual const G4String& GetRelativeRotMatName() const;
    virtual G4ThreeVector GetRelativePlace() const;

  private:

    void FillSolidParams( const std::vector<G4String>&  wl );

  protected:

    G4String theName;   
      // Name of the solid
    G4String theType;   
      // Type of the solid (Simple, Boolean, Box, Tube, ...)
    std::vector< std::vector<G4double>* > theSolidParams; 
      // Vectors of parameters
};

#endif
