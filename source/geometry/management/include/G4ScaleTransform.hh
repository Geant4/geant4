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
// G4ScaleTransform
//
// Class description:
//
// A class for geometric scaling transformations.
// Supports efficient arbitrary transformation of points, vectors and
// normals and the computation of compound & inverse transformations.
//
// Interfaces to the CLHEP class G4ThreeVector
//
// For member function descriptions, see comments by declarations. For
// additional clarification, also check the `const' declarations for
// functions & their parameters.
//
// Member data:
//
//    G4ThreeVector fScale;  // scale transformation 
//    G4ThreeVector fIScale; // inverse scale (avoid divisions)
//    G4double flFactor;     // factor for conversion to local frame
//    G4double fgFactor;     // factor for conversion to global frame

// Author: Gabriele Cosmo (CERN), 18.02.2016 - Initial version
//         Evgueni Tcherniaev (CERN), 11.03.2016 - Added normals transforms
// --------------------------------------------------------------------
#ifndef G4SCALETRANSFORM_HH
#define G4SCALETRANSFORM_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

/**
 * @brief G4ScaleTransform is a class for geometric scaling transformations.
 * It supports efficient arbitrary transformation of points, vectors and
 * normals and the computation of compound and inverse transformations.
 */

class G4ScaleTransform
{
  public:

    /**
     * Default Constructor.
     */
    inline G4ScaleTransform();
    
    /**
     * Constructor with scale parameters on each axis.
     *  @param[in] sx Scaling in X.
     *  @param[in] sy Scaling in Y.
     *  @param[in] sz Scaling in Z.
     */
    inline G4ScaleTransform(G4double sx, G4double sy, G4double sz);

    /**
     * Constructor taking a 3-vector.
     *  @param[in] scale Scaling transformation.
     */
    inline G4ScaleTransform(const G4ThreeVector& scale);

    /**
     * Constructor taking a Scale3D.
     *  @param[in] scale Scaling transformation.
     */
    inline G4ScaleTransform(const G4Scale3D& scale);
        
    /**
     * Copy constructor and assignment operator.
     */
    inline G4ScaleTransform(const G4ScaleTransform& right);
    inline G4ScaleTransform& operator=(const G4ScaleTransform& right);

    /**
     * Updates the backed-up inverse scale and special conversion factors
     * based on the values of the scale. Needed at initialisation and
     * whenever the scale has changed value.
     */
    inline void Init();

    /**
     * Accessors returning a reference to the scale and inverse scale
     * transformations.
     */
    inline const G4ThreeVector& GetScale() const;
    inline const G4ThreeVector& GetInvScale() const;
 
    /**
     * Modifiers for the scale transformation.
     */
    inline void SetScale(const G4ThreeVector& scale);
    inline void SetScale(const G4Scale3D& scale);
    inline void SetScale(G4double sx, G4double sy, G4double sz);

    /**
     * Methods to transform a point from global to local frame.
     */
    inline void Transform(const G4ThreeVector& global,
                                G4ThreeVector& local) const;  
    inline G4ThreeVector Transform(const G4ThreeVector& global) const;

    /**
     * Methods to transform a point from local to global frame.
     */
    inline void InverseTransform(const G4ThreeVector& local, 
                                       G4ThreeVector& global) const;
    inline G4ThreeVector InverseTransform(const G4ThreeVector& local) const;

    /**
     * Methods to transform a normal from global to local frame.
     */
    inline void TransformNormal(const G4ThreeVector& global,
                                      G4ThreeVector& local) const;  
    inline G4ThreeVector TransformNormal(const G4ThreeVector& global) const;

    /**
     * Methods to transform a normal from local to global frame.
     */
    inline void InverseTransformNormal(const G4ThreeVector& local, 
                                             G4ThreeVector& global) const;
    inline G4ThreeVector InverseTransformNormal(const G4ThreeVector& local) const;

    /**
     * Transforms a distance 'dist' along a given direction 'dir'
     * from global to local frame.
     */
    inline G4double TransformDistance(G4double dist,
                                      const G4ThreeVector& dir) const;

    /**
     * Transforms a 'safety' distance from global to local frame (conservative).
     */
    inline G4double TransformDistance(G4double safety) const;

    /**
     * Transforms a distance 'dist' along a given direction 'dir'
     * from local to global frame.
     */
    inline G4double InverseTransformDistance(G4double dist,
                                             const G4ThreeVector& dir) const;

    /**
     * Transforms a 'safety' distance from local to global frame (conservative).
     */
    inline G4double InverseTransformDistance(G4double safety) const;

  private:

    G4ThreeVector fScale;  // scale transformation 
    G4ThreeVector fIScale; // inverse scale (avoid divisions)

    /** Conversion factors to local/global frames. */
    G4double flFactor = 1.0, fgFactor = 1.0;
};

std::ostream& operator<<(std::ostream& os, const G4ScaleTransform& scale);

#include "G4ScaleTransform.icc"

#endif
