// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4SetRotation.hh,v 1.1 1999-01-07 16:06:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _G3TOG4ROTATION_
#define _G3TOG4ROTATION_

#include "G4RotationMatrix.hh"
#include "globals.hh"

class G3toG4SetRotation : public G4RotationMatrix 
{
public:
    G3toG4SetRotation();

    G3toG4SetRotation(const G4ThreeVector& col1,
                      const G4ThreeVector& col2,
                      const G4ThreeVector& col3);
    
    ~G3toG4SetRotation(){;}
    
    inline void Setxx(const G4double axx) {rxx = axx;
    }
    
    inline void Setxy(const G4double axy) {rxy = axy;
    }
    
    inline void Setxz(const G4double axz) {rxz = axz;
    }
    
    inline void Setyx(const G4double ayx) {ryx = ayx;
    }
    
    inline void Setyy(const G4double ayy) {ryy = ayy;
    }
    
    inline void Setyz(const G4double ayz) {ryz = ayz;
    }
    
    inline void Setzx(const G4double azx) {rzx = azx;
    }
    
    inline void Setzy(const G4double azy) {rzy = azy;
    }
    
    inline void Setzz(const G4double azz) {rzz = azz;
    }
};

#endif    
