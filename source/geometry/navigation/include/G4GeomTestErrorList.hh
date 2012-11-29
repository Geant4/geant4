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
// $Id$
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestErrorList
//
// Class description:
//
// A list of line segments that are found inside two
// separate daughter volumes, indicating a geometry error.
//
// This class relies on the compiler generated copy constructor
// and assignment operator.

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestErrorList_hh
#define G4GeomTestErrorList_hh

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;

class G4GeomTestErrorList
{
  public:  // with description
  
    G4GeomTestErrorList( const G4VPhysicalVolume *theMother );
    virtual ~G4GeomTestErrorList();
      // Constructor and virtual destructor

    void AddError( const G4ThreeVector &s1, const G4ThreeVector &s2  );
      // Declare a new instance of an error, by specifying
      // two points in the coordinate system of the mother

    const G4VPhysicalVolume *GetMother() const;
      // Return pointers to mother volume

    G4int NumError() const;
      // Return number of errors

    void GetMotherPoints( G4int i, G4ThreeVector &s1, G4ThreeVector &s2 ) const;
    void GetGlobalPoints( G4int i, G4ThreeVector &s1, G4ThreeVector &s2 ) const;
      // Return start and end points in various
      // coordinate systems

    void GetOneDaughtPoints( const G4VPhysicalVolume *daught,
                          G4int i, G4ThreeVector &s1, G4ThreeVector &s2 ) const;
      // Return start and end points in the coordinate system of a 
      // daughter volume

  private:
  
    void FindGlobalCoordinateSystem();
      // Calculate the global coordinate system

    class Segment
    {
      public:
        Segment( const G4ThreeVector &as1, const G4ThreeVector &as2 )
          : s1(as1), s2(as2) {;}
    
        const G4ThreeVector &GetS1() const { return s1; }
        const G4ThreeVector &GetS2() const { return s2; }
    
      private:
        G4ThreeVector s1, s2;
    };
    std::vector<Segment> segments;
      // List of error segments

    const G4VPhysicalVolume *mother;
      // Mother volume

    G4ThreeVector globalTranslation;
    G4RotationMatrix globalRotation;
      // Global coordinate system with respect to mother
};

#endif
