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
// $Id: G4GeomTestErrorList.hh,v 1.1 2001/10/17 12:59:53 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-00 $
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

#ifndef G4GeomTestErrorList_hh
#define G4GeomTestErrorList_hh

#include "g4std/vector"

#include "globals.hh"
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
    G4std::vector<Segment> segments;
      // List of error segments

    const G4VPhysicalVolume *mother;
      // Mother volume

    G4ThreeVector globalTranslation;
    G4RotationMatrix globalRotation;
      // Global coordinate system with respect to mother
};

#endif
