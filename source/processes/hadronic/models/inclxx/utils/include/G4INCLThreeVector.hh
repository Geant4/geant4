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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * ThreeVector.hh
 *
 *  \date 4 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLThreeVector_hh
#define G4INCLThreeVector_hh 1

#include <string>
#include <sstream>
#include <cmath>

namespace G4INCL {

  class ThreeVector {
    public:
      ThreeVector()
        :x(0.0), y(0.0), z(0.0)
      {}

      ThreeVector(G4double ax, G4double ay, G4double az)
        :x(ax), y(ay), z(az)
      {}

      inline G4double getX() const { return x; }
      inline G4double getY() const { return y; }
      inline G4double getZ() const { return z; }

      inline G4double perp() const { return std::sqrt(x*x + y*y); }
      inline G4double perp2() const { return x*x + y*y; }
      /**
       * Get the length of the vector.
       */
      inline G4double mag() const { return std::sqrt(x*x + y*y + z*z); }

      /**
       * Get the square of the length.
       */
      inline G4double mag2() const { return (x*x + y*y + z*z); }

      /**
       * Theta angle
       */
      inline G4double theta() const {
        return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : std::atan2(perp(),z);
      }

      /**
       * Phi angle
       */
      inline G4double phi() const {
        return x == 0.0 && y == 0.0 ? 0.0 : std::atan2(y,x);
      }

      /**
       * Dot product.
       */
      inline G4double dot(const ThreeVector &v) const {
        return (x*v.x + y*v.y + z*v.z);
      }

      /**
       * Vector product.
       */
      ThreeVector vector(const ThreeVector &v) const {
        return ThreeVector(
            y*v.z - z*v.y,
            z*v.x - x*v.z,
            x*v.y - y*v.x
            );
      }

      /// \brief Set the x coordinate
      inline void setX(G4double ax) { x =  ax; }

      /// \brief Set the y coordinate
      inline void setY(G4double ay) { y =  ay; }

      /// \brief Set the z coordinate
      inline void setZ(G4double az) { z =  az; }

      /// \brief Set all the coordinates
      inline void set(const G4double ax, const G4double ay, const G4double az) { x=ax; y=ay; z=az; }

      inline void operator+= (const ThreeVector &v) {
        x += v.x;
        y += v.y;
        z += v.z;
      }

      /// \brief Unary minus operator
      inline ThreeVector operator- () const {
        return ThreeVector(-x,-y,-z);
      }

      inline void operator-= (const ThreeVector &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
      }

      template<typename T>
        inline void operator*= (const T &c) {
          x *= c;
          y *= c;
          z *= c;
        }

      template<typename T>
        inline void operator/= (const T &c) {
          const G4double oneOverC = 1./c;
          this->operator*=(oneOverC);
        }

      inline ThreeVector operator- (const ThreeVector &v) const {
        return ThreeVector(x-v.x, y-v.y, z-v.z);
      }

      inline ThreeVector operator+ (const ThreeVector &v) const {
        return ThreeVector(x+v.x, y+v.y, z+v.z);
      }

      /**
       * Divides all components of the vector with a constant number.
       */
      inline ThreeVector operator/ (const G4double C) const {
        const G4double oneOverC = 1./C;
        return ThreeVector(x*oneOverC, y*oneOverC, z*oneOverC);
      }

      inline ThreeVector operator* (const G4double C) const {
        return ThreeVector(x*C, y*C, z*C);
      }

      /** \brief Rotate the vector by a given angle around a given axis
       *
       * \param angle the rotation angle
       * \param axis the rotation axis, which must be a unit vector
       */
      inline void rotate(const G4double angle, const ThreeVector &axis) {
        // Use Rodrigues' formula
        const G4double cos = std::cos(angle);
        const G4double sin = std::sin(angle);
        (*this) = (*this) * cos + axis.vector(*this) * sin + axis * (axis.dot(*this)*(1.-cos));
      }

      /** \brief Return a vector orthogonal to this
       *
       * Simple algorithm from Hughes and Moeller, J. Graphics Tools 4 (1999)
       * 33.
       */
      ThreeVector anyOrthogonal() const {
        if(x<=y && x<=z)
          return ThreeVector(0., -z, y);
        else if(y<=x && y<=z)
          return ThreeVector(-z, 0., x);
        else
          return ThreeVector(-y, x, 0.);
      }

      std::string print() const {
        std::stringstream ss;
        ss <<"(x = " << x << "   y = " << y << "   z = " << z <<")";
        return ss.str();
      }

      std::string dump() const {
        std::stringstream ss;
        ss <<"(vector3 " << x << " " << y << " " << z << ")";
        return ss.str();
      }

    private:
      G4double x, y, z; //> Vector components
  };

}

#endif
