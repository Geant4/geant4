#include "G3toG4SetRotation.hh"

G3toG4SetRotation::G3toG4SetRotation()
{
    rxx = 1;
    ryx = 0;
    rzx = 0;
    rxy = 0;
    ryy = 1;
    rzy = 0;
    
        rxz = 0;
        ryz = 0;
        rzz = 1;
    
    
}
G3toG4SetRotation::G3toG4SetRotation(const G4ThreeVector& col1,
                                     const G4ThreeVector& col2,
                                     const G4ThreeVector& col3)
{
    rxx = col1.x();
    ryx = col1.y();
    rzx = col1.z();
    
    rxy = col2.x();
    ryy = col2.y();
    rzy = col2.z();
    
    rxz = col3.x();
    ryz = col3.y();
    rzz = col3.z();
    
}
