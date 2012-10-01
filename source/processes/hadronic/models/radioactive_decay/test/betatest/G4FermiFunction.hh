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
#ifndef G4FermiFunction_h
#define G4FermiFunction_h 1

//#include "globals.hh"

class G4FermiFunction
{
  // class  description
  // It is to calculate the Coulomb correction to beta particles
  // 

public: // with description
  
  G4FermiFunction(int const fA, int const fZ) :
    A(fA), Z(fZ)
  {};
  // constructor: fA the daughter nucleus mass.
  //              fZ the daughter nucleus charge. Negative value for 
  //                 beta+ decays.
  //
  ~G4FermiFunction() ;
  // desctructor
  //
  double GetFF(const double E);
  // Returns the Fermi factor at energy E.
  // E total energy of the beta particle in unit of Me.
  //
  double GetFFN(const double E0);
  // Returns the Fermi factor normalisation, i.e. the maximum
  // value for beta decay with end-point energy E0.
  // E0 is the total energy including beta particle rest mass in unit
  // of Me.
  //
private:

  int A;
  int Z;
  
  static const double   PI = 3.14159;
  /*
  //              COEFFICIENTS FOR MINIMAX 
  //              APPROXIMATION TO GAMMA(X),
  //              2.0 .LE. X .LE. 3.0
  static const G4double P[5] = {-51.49952, 80.05398,-201.4659,-1.889439, 9.895546};
  static const G4double Q[4] = {130.5263, -303.5898, 26.84174, -19.52375};
  //                                  APPROXIMATION TO LN(GAMMA(X)),
  //                               12.0 .LE. X
  static const G4double P4[3]={.9189385, .8333332E-01, -.2770927E-02};
  //
  static const G4int IEND = 4 , IEND1 = 3,  IEND2 = 2;
  static const G4double  XINF = 1.7E+38;
  //  GAMMA(XMIN) .APPROX. XINF
  //      GAMMA(BIG1) .APPROX. XINF
  static const G4double             XMIN = 5.8775E-39;
  static const G4double             BIG1 = 34.844; */

private:

  double Gamma(double X);
  
};
#endif
 






