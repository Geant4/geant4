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
// G4LogicalSurface
//
// Class description:
//
// An abstraction of a geometrical surface, it is an abstract 
// base class for different implementations of surfaces.
// Its primary function is to hold pointers to objects that describe the
// surface's physical properties. For example it holds a pointer to a
// surface's optical properties, and because of this it is used in processes
// like G4OpBoundaryProcess. 
//
// Methods:
//   G4SurfaceProperty*  GetSurfaceProperty() const
//   void     SetSurfaceProperty(G4SurfaceProperty*)
//
//   G4String GetName() const
//   void     SetName(const G4String&)
//
//   G4TransitionRadiationSurface*  GetTransitionRadiationSurface() const
//   void SetTransitionRadiationSurface(G4TransitionRadiationSurface*)
//
// Data members:
//   G4String                       theName
//   G4SurfaceProperty*             theSurfaceProperty
//   G4TransitionRadiationSurface*  theTransRadSurface

// Authors: John Apostolakis (CERN) & Peter Gumplinger (TRIUMF), 26.06.1997
// ------------------------------------------------------------------------
#ifndef G4LOGICALSURFACE_HH
#define G4LOGICALSURFACE_HH

#include "G4Types.hh"
#include "G4String.hh"

class G4SurfaceProperty;
class G4TransitionRadiationSurface;

/**
 * @brief G4LogicalSurface is an abstraction of a geometrical surface, it is
 * an abstract base class for different implementations of surfaces.
 * Its primary function is to hold pointers to objects that describe the
 * surface's physical properties. For example it holds a pointer to a
 * surface's optical properties, and because of this it is used in processes
 * like G4OpBoundaryProcess.
 */

class G4LogicalSurface
{
  public:

    /**
     * Default Destructor.
     */
    virtual ~G4LogicalSurface() = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4LogicalSurface(const G4LogicalSurface&) = delete;
    G4LogicalSurface& operator=(const G4LogicalSurface&) = delete;

    /**
     * Getter and setter for the surface property object.
     */
    inline G4SurfaceProperty* GetSurfaceProperty() const;
    inline void SetSurfaceProperty(G4SurfaceProperty* ptrSurfaceProperty);

    /**
     * Getter and setter for the surface name.
     */
    inline const G4String& GetName() const;
    inline void SetName(const G4String& name);

    /**
     * Getter and setter for the transition radiation surface attribute.
     */
    inline G4TransitionRadiationSurface* GetTransitionRadiationSurface() const;
    inline void SetTransitionRadiationSurface(G4TransitionRadiationSurface* trs);

    /**
     * Equality and disequality operators.
     */
    inline G4bool operator==(const G4LogicalSurface &right) const;
    inline G4bool operator!=(const G4LogicalSurface &right) const;

  protected:

    /**
     * Protected constructor. There should be no free instances of this class.
     */
    G4LogicalSurface(const G4String& name, G4SurfaceProperty* prop); 

  private:

    /** Surface name. */
    G4String theName;

    G4SurfaceProperty* theSurfaceProperty = nullptr;
    G4TransitionRadiationSurface* theTransRadSurface = nullptr;
};

#include "G4LogicalSurface.icc"

#endif
