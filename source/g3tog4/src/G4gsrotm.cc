// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsrotm.cc,v 1.8 1999-12-15 14:49:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3toG4.hh"
#include "G3RotTable.hh"
#include "G4ThreeVector.hh"
#include "G3toG4RotationMatrix.hh"

void PG4gsrotm(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsrotm);

    // interpret the parameters
    G4int irot = Ipar[0];
    
    // the angles in Geant are in degrees
    G4double theta1 = Rpar[0];
    G4double phi1   = Rpar[1];
    G4double theta2 = Rpar[2];
    G4double phi2   = Rpar[3];
    G4double theta3 = Rpar[4];
    G4double phi3   = Rpar[5];

    G4gsrotm(irot, theta1,phi1, theta2,phi2, theta3,phi3);
}

void G4gsrotm(G4int irot, G4double theta1, G4double phi1,
              G4double theta2, G4double phi2, G4double theta3, G4double phi3)
{
    G4double degrad = pi/180;
    
    G4double th1r = theta1*degrad;
    G4double th2r = theta2*degrad;
    G4double th3r = theta3*degrad;
    
    G4double phi1r = phi1*degrad;
    G4double phi2r = phi2*degrad;
    G4double phi3r = phi3*degrad;
    
        // Construct unit vectors 
    
    G4ThreeVector x(sin(th1r)*cos(phi1r), sin(th1r)*sin(phi1r), cos(th1r));
    G4ThreeVector y(sin(th2r)*cos(phi2r), sin(th2r)*sin(phi2r), cos(th2r));
    G4ThreeVector z(sin(th3r)*cos(phi3r), sin(th3r)*sin(phi3r), cos(th3r));

        // check for orthonormality and left-handedness

    G4double check = (x.cross(y))*z;
    G4double tol = 1.0e-3;
        
    if (1-abs(check)>tol) {
        G4cerr << "Coordinate axes forming rotation matrix "
               << irot << " are not orthonormal.(" << 1-abs(check) << ")" 
         << G4endl;
        G4cerr << " theta1=" << theta1;
        G4cerr << " phi1=" << phi1;
        G4cerr << " theta2=" << theta2;
        G4cerr << " phi2=" << phi2;
        G4cerr << " theta3=" << theta3;
        G4cerr << " phi3=" << phi3;
        G4cerr << G4endl;
        G4Exception("G4gsrotm error");
    }
    else if (1+check<=tol) {
        G4cerr << "G4gsrotm warning: coordinate axes forming rotation "
               << "matrix " << irot << " are left-handed" << G4endl;
    }
    
    G3toG4RotationMatrix* rotp = new G3toG4RotationMatrix;

    rotp->SetRotationMatrixByRow(x, y, z);
    
        // add it to the List

    G3Rot.Put(irot, rotp);
}
