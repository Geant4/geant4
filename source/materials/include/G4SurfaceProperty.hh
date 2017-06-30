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
// $Id: G4SurfaceProperty.hh 103256 2017-03-23 08:53:38Z gcosmo $
//
// 
////////////////////////////////////////////////////////////////////////
// G4SurfaceProperty Definition
////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// A base class describing a surface property.
// Derived classes are G4Opticalsurface, G4Firovsurface, etc.      
// Contains the enumeration G4SurfaceType.

// File:        G4SurfaceProperty.hh
// Description: A base class for for descriping surface property such
//              as G4OpticalSurface, G4FirsovSurface, G4X-raySurface
// Version:     1.0
// Created:     13-10-2003
// Author:      Fan Lei
// Updated:     Mariele Stockhoff 2017-02-24 add DAVIS model 
////////////////////////////////////////////////////////////////////////

#ifndef G4SurfaceProperty_h
#define G4SurfaceProperty_h 1

/////////////
// Includes
/////////////

#include <vector>

#include "G4Types.hh"
#include "G4String.hh"

class G4SurfaceProperty;

typedef std::vector<G4SurfaceProperty*> G4SurfacePropertyTable;

enum G4SurfaceType
{
   dielectric_metal,            // dielectric-metal interface
   dielectric_dielectric,       // dielectric-dielectric interface
   dielectric_LUT,              // dielectric-Look-Up-Table interface
   dielectric_LUTDAVIS,         // dielectric-Look-Up-Table DAVIS interface
   dielectric_dichroic,         // dichroic filter interface
   firsov,                      // for Firsov Process
   x_ray                        // for x-ray mirror process
};

/////////////////////
// Class Definition
/////////////////////

class G4SurfaceProperty
{
  public: // Without description

     //////////////
     // Operators
     //////////////

    // G4SurfaceProperty(const G4SurfaceProperty &right);
    // const G4SurfaceProperty & operator=(const G4SurfaceProperty &right);

    // G4int operator==(const G4SurfaceProperty &right) const;
    // G4int operator!=(const G4SurfaceProperty &right) const;

  public: // With description

    ////////////////////////////////
    // Constructors and Destructor
    ////////////////////////////////

    G4SurfaceProperty(const G4String& name, G4SurfaceType type = x_ray);
    // Constructor of a X-ray optical surface object.

  public: // Without description

    G4SurfaceProperty();
    virtual ~G4SurfaceProperty();

    ////////////
    // Methods
    ////////////

  public: // With description

    const G4String& GetName() const { return theName; }
    // Returns the surface name.
    void     SetName(const G4String& name) { theName = name; }
    // Sets the surface name.

    const G4SurfaceType& GetType() const { return theType; }
    // Returns the surface type.
    void     SetType(const G4SurfaceType& type) { theType = type; }
    // Sets the surface type.        

    static void CleanSurfacePropertyTable();
    static const G4SurfacePropertyTable* GetSurfacePropertyTable();
    static size_t GetNumberOfSurfaceProperties();
    static void DumpTableInfo();
    // To handle the table of surface properties.

  protected:

    // ------------------
    // Basic data members ( To define surface property)
    // ------------------

    G4String theName;                // Surface name

    G4SurfaceType theType;           // Surface type

    static G4SurfacePropertyTable theSurfacePropertyTable;
    // The static Table of SurfaceProperties.
};

#endif /* G4SurfaceProperty_h */
