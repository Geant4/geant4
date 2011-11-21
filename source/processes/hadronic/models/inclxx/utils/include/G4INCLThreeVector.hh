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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * ThreeVector.hh
 *
 *  Created on: 4 June 2009
 *      Author: Pekka Kaitaniemi
 */

#ifndef G4INCLThreeVector_hh
#define G4INCLThreeVector_hh 1

#include <string>
#include <sstream>
#include <cmath>

namespace G4INCL {

#ifdef INCLXX_IN_GEANT4_MODE
#define INCL_SIMPLEVECTOR_USE 1
#endif // INCLXX_IN_GEANT4_MODE

#ifdef INCL_SIMPLEVECTOR_USE // Use simple vector
  class ThreeVector {
  public:
    ThreeVector()
      :x(0.0), y(0.0), z(0.0)
    {};

    ThreeVector(G4double ax, G4double ay, G4double az)
      :x(ax), y(ay), z(az)
    {};

    ThreeVector(const ThreeVector& v)
      :x(v.getX()), y(v.getY()), z(v.getZ()) {};

    virtual ~ThreeVector() {};

    inline G4double getX() const { return x; };
    inline G4double getY() const { return y; };
    inline G4double getZ() const { return z; };

    inline G4double perp() const { return std::sqrt(x*x + y*y); };
    inline G4double perp2() const { return x*x + y*y; };
    /**
     * Get the length of the vector.
     */
    inline G4double mag() const { return std::sqrt(x*x + y*y + z*z); };

    /**
     * Get the square of the length.
     */
    inline G4double mag2() const { return (x*x + y*y + z*z); };

    /**
     * Theta angle
     */
    inline G4double theta() const {
      return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : std::atan2(perp(),z);
    };

    /**
     * Phi angle
     */
    inline G4double phi() const {
      return x == 0.0 && y == 0.0 ? 0.0 : std::atan2(y,x);
    };

    /**
     * Dot product.
     */
    inline G4double dot(const ThreeVector &v) const { return (x*v.getX() + y*v.getY() + z*v.getZ()); };

    /**
     * Vector product.
     */
    ThreeVector vector(const ThreeVector &v) const {
      return ThreeVector(
          y*v.getZ() - z*v.getY(),
          z*v.getX() - x*v.getZ(),
          x*v.getY() - y*v.getX()
          );
    }

    /// \brief Set the x coordinate
    inline void setX(G4double ax) { x =  ax; }

    /// \brief Set the y coordinate
    inline void setY(G4double ay) { y =  ay; }

    /// \brief Set the z coordinate
    inline void setZ(G4double az) { z =  az; }

    inline void operator+= (const ThreeVector &v) {
      x += v.getX();
      y += v.getY();
      z += v.getZ();
    }

    /// \brief Unary minus operator
    inline ThreeVector operator- () const {
      return ThreeVector(-x,-y,-z);
    }

    inline void operator-= (const ThreeVector &v) {
      x -= v.getX();
      y -= v.getY();
      z -= v.getZ();
    }

    template<typename T>
    inline void operator*= (const T c) {
      x *= c;
      y *= c;
      z *= c;
    }

    template<typename T>
    inline void operator/= (const T c) {
      x /= c;
      y /= c;
      z /= c;
    }

    inline ThreeVector operator- (const ThreeVector &v) const {
      ThreeVector w(*this);
      w -= v;
      return w;
    };

    inline ThreeVector operator+ (const ThreeVector &v) const {
      ThreeVector w(*this);
      w += v;
      return w;
    };

    /**
     * Divides all components of the vector with a constant number.
     */
    inline ThreeVector operator/ (const G4double C) const {
      ThreeVector w(*this);
      w /= C;
      return w;
    };

    inline ThreeVector operator* (const G4double C) const {
      ThreeVector w(*this);
      w *= C;
      return w;
    };

    std::string prG4int() const {
      std::stringstream ss;
      ss <<"(x = " << x << "   y = " << y << "   z = " << z <<")";
      return ss.str();
    };

    std::string dump() const {
      std::stringstream ss;
      ss <<"(vector3 " << x << " " << y << " " << z << ")";
      return ss.str();
    }

  private:
    G4double x, y, z; //> Vector components
  };

#else // Use expression templates
  // Vector argument
  template< class ta_a >
  class vecarg {
    const ta_a& argv;
  public:
    inline vecarg(const ta_a &a)
      :argv(a) {};
    inline G4double eval(const G4int i) const {
      return argv.eval(i);
    };
  };
  
  template<>
  class vecarg< const G4double > {
    const G4double &argv;
  public:
    inline vecarg(const G4double &a)
      :argv(a)
    {};
    inline const G4double eval(const G4int i) const {
      return argv;
    };
  };

  // Expressions
  template< class ta_a, class ta_b, class ta_eval >
  class vecexp_bin {
    const vecarg<ta_a> arg0;
    const vecarg<ta_b> arg1;

  public:
    inline vecexp_bin(const ta_a &a0, const ta_b &a1)
      :arg0(a0), arg1(a1)
    {};
    inline const G4double eval(const G4int i) const {
      return ta_eval::eval(i, arg0, arg1);
    };
  };

  template< class ta_a, class ta_eval >
  class vecexp_unary {
    const vecarg<ta_a> arg0;

  public:
    inline vecexp_unary(const ta_a &a0)
      :arg0(a0)
    {};
    inline const G4double eval(const G4int i) const {
      return ta_eval::eval(i, arg0.eval(i));
    };
  };


  // Sum

  struct sum {
    template< class ta_a, class ta_b >
    inline static const G4double eval(const G4int i, const ta_a &A, const ta_b &B) {
      return A.eval(i) + B.eval(i);
    }
  };

  template< class ta_c1, class ta_c2>
  inline const vecexp_bin< const ta_c1, const ta_c2, sum >
  operator + (const ta_c1 &Pa, const ta_c2 &Pb) {
    return vecexp_bin< const ta_c1, const ta_c2, sum>(Pa, Pb);
  }

  // Subtraction
  struct subtraction {
    template< class ta_a, class ta_b >
    inline static const G4double eval(const G4int i, const ta_a &A, const ta_b &B) {
      return A.eval(i) - B.eval(i);
    }
  };

  template< class ta_c1, class ta_c2 >
  inline const vecexp_bin< const ta_c1, const ta_c2, subtraction >
  operator - (const ta_c1 &Pa, const ta_c2 &Pb) {
    return vecexp_bin< const ta_c1, const ta_c2, subtraction >(Pa, Pb);
  }

  // Multiplication
  struct multiplication {
    template< class ta_a, class ta_b>
    inline static const G4double eval(const G4int i, const ta_a &A, const ta_b &B) {
      return A.eval(i) * B.eval(i);
    }
  };

  template< class ta_c1, class ta_c2>
  inline const vecexp_bin< const ta_c1, const ta_c2, multiplication >
  operator * (const ta_c1 &Pa, const ta_c2 &Pb) {
    return vecexp_bin< const ta_c1, const ta_c2, multiplication>(Pa, Pb);
  }

  // Division
  struct division {
    template< class ta_a, class ta_b>
    inline static const G4double eval(const G4int i, const ta_a &A, const ta_b &B) {
      return A.eval(i) / B.eval(i);
    }
  };

  template< class ta_c1, class ta_c2>
  inline const vecexp_bin< const ta_c1, const ta_c2, division >
  operator / (const ta_c1 &Pa, const ta_c2 &Pb) {
    return vecexp_bin< const ta_c1, const ta_c2, division>(Pa, Pb);
  }

  struct ThreeVector {
    G4double x, y, z;

    ThreeVector()
      :x(0.0), y(0.0), z(0.0)
    {};

    ThreeVector(const G4double X, const G4double Y, const G4double Z)
    {
      x = X; y = Y; z = Z;
    };

    template< class ta>
    inline ThreeVector(const ta &expr)
      :x(expr.eval(0)), y(expr.eval(1)), z(expr.eval(2))
    {
      //      x = expr.eval(0); y = expr.eval(1); z = expr.eval(2);
    }

    inline G4double eval(const G4int i) const {
      if(i == 0) return x;
      else if(i == 1) return y;
      else if(i == 2) return z;
      return 0.0; // Undefined element index
    }

    template< class ta >
    inline void operator = (const ta &expr) {
      x = expr.eval(0);
      y = expr.eval(1);
      z = expr.eval(2);
    }

    template< class ta >
    inline void operator+= (const ta &expr) {
      x += expr.eval(0); y += expr.eval(1); z += expr.eval(2);
    }

    template< class ta >
    inline void operator-= (const ta &expr) {
      x -= expr.eval(0); y -= expr.eval(1); z -= expr.eval(2);
    }

    template< class ta >
    inline void operator/= (const ta &expr) {
      x /= expr.eval(0); y /= expr.eval(1); z /= expr.eval(2);
    }

    template< class ta >
    inline const G4double dot(const ta &expr) const {
      return x * expr.eval(0) + y * expr.eval(1) + z * expr.eval(2);
    }

    inline G4double perp() const {
      return std::sqrt(x*x + y*y);
    };

    inline G4double perp2() const {
      return x*x + y*y;
    };

    inline G4double mag2() const {
      return x*x + y*y + z*z;
    }

    inline G4double mag() const {
      return std::sqrt(mag2());
    }

    inline G4double getX() const {
      return x;
    }

    inline G4double getY() const {
      return x;
    }

    inline G4double getZ() const {
      return x;
    }

    inline void setX(G4double ax) {
      x =  ax;
    }

    inline void setY(G4double ay) {
      y =  ay;
    }

    inline void setZ(G4double az) {
      z =  az;
    }

    inline G4double theta() const {
      return std::acos(z/mag2());
    }

    inline G4double phi () const {
      return std::atan2(y, x);
    }

    std::string prG4int() const {
      std::stringstream ss;
      ss <<" (x = " << x << " y = " << y << " z = " << z << ")" << std::endl;
      return ss.str();
    }

    std::string dump() const {
      std::stringstream ss;
      ss <<"(vector3 " << x << " " << y << " " << z << ")";
      return ss.str();
    }
  };
#endif
}

#endif
