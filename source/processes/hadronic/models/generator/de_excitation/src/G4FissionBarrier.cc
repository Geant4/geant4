// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FissionBarrier.cc,v 1.1 1999/01/07 16:11:54 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)


#include "G4FissionBarrier.hh"

G4FissionBarrier::G4FissionBarrier(const G4FissionBarrier & right)
{
  G4Exception("G4FissionBarrier::copy_constructor meant to not be accessable.");
}


const G4FissionBarrier & G4FissionBarrier::operator=(const G4FissionBarrier & right)
{
 G4Exception("G4FissionBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4FissionBarrier::operator==(const G4FissionBarrier & right) const 
{
 return false;
}

G4bool G4FissionBarrier::operator!=(const G4FissionBarrier & right) const 
{
 return true;
}



G4double G4FissionBarrier::FissionBarrier(const G4int A, const G4int Z)
  // Compute fission barrier according with Barashenkov's prescription for A >= 65
{
  if (A >= 65) return BarashenkovFissionBarrier(A,Z)*MeV;
  else return 1.0*GeV;
}


G4double G4FissionBarrier::BarashenkovFissionBarrier(const G4int A, const G4int Z)
  // Calculates Fission Barrier heights (in MeV), which are function of the nuclear 
  // fissility parameteter x = Z*Z/A. 
  // Barashenkov V., Iljinov A. and Toneev V. parametrization (1972)
{
  const G4double ShellCorr1[130] = {
    20.80,  15.80,  21.00,  16.80,  19.80,
    16.50,  18.80,  16.50,  18.50,  17.20,
    18.26,  15.05,  16.01,  12.04,  13.27,
    11.09,  12.17,  10.26,  11.04,   8.41,
    9.79,   7.36,   8.15,   5.63,   5.88,
    3.17,   3.32,   0.82,   1.83,   0.97,
    2.33,   1.27,   2.92,   1.61,   2.91,
    1.35,   2.40,   0.89,   1.74,   0.36,
    0.95,  -0.65,  -0.04,  -1.73,  -0.96,
    -2.87,  -2.05,  -4.05,  -3.40,  -5.72,
    -3.75,  -4.13,  -2.42,  -2.85,  -1.01,
    -1.33,   0.54,  -0.02,   1.74,   0.75,
    2.24,   1.00,   1.98,   0.79,   1.54,
    0.39,   1.08,   0.00,   0.78,  -0.35,
    0.58,  -0.55,   0.59,  -0.61,   0.59,
    -0.35,   0.32,  -0.96,  -0.52,  -2.08,
    -2.46,  -3.64,  -1.55,  -0.96,   0.97,
    0.88,   2.37,   1.75,   2.72,   1.90,
    2.55,   1.46,   1.93,   0.86,   1.17,
    0.08,   0.39,  -0.76,  -0.39,  -1.51,
    -1.17,  -2.36,  -1.95,  -3.06,  -2.62,
    -3.55,  -2.95,  -3.75,  -3.07,  -3.79,
    -3.06,  -3.77,  -3.05,  -3.78,  -3.12,
    -3.90,  -3.35,  -4.24,  -3.86,  -4.92,
    -5.06,  -6.77,  -7.41,  -9.18, -10.16,
    -11.12,  -9.76,  -9.23,  -7.96,  -7.65};
  const G4double ShellCorr2[200] = {
    -8.40, -12.90,  -8.00, -11.90,  -9.20,
    -12.50, -10.80, -13.60, -11.20, -12.20,
    -12.81, -15.40, -13.07, -15.80, -13.81,
    -14.98, -12.63, -13.76, -11.37, -12.38,
    -9.23,  -9.65,  -7.64,  -9.17,  -8.05,
    -9.72,  -8.87, -10.76,  -8.64,  -8.89,
    -6.60,  -7.13,  -4.77,  -5.33,  -3.06,
    -3.79,  -1.72,  -2.79,  -0.93,  -2.19,
    -0.52,  -1.90,  -0.45,  -2.20,  -1.22,
    -3.07,  -2.42,  -4.37,  -3.94,  -6.08,
    -4.49,  -4.50,  -3.14,  -2.93,  -1.04,
    -1.36,   0.69,   0.21,   2.11,   1.33,
    3.29,   2.46,   4.30,   3.32,   4.79,
    3.62,   4.97,   3.64,   4.63,   3.07,
    4.06,   2.49,   3.30,   1.46,   2.06,
    0.51,   0.74,  -1.18,  -1.26,  -3.54,
    -3.97,  -5.26,  -4.18,  -3.71,  -2.10,
    -1.70,  -0.08,  -0.18,   0.94,   0.27,
    1.13,   0.08,   0.91,  -0.31,   0.49,
    -0.78,   0.08,  -1.15,  -0.23,  -1.41,
    -0.42,  -1.55,  -0.55,  -1.66,  -0.66,
    -1.73,  -0.75,  -1.74,  -0.78,  -1.69,
    -0.78,  -1.60,  -0.75,  -1.46,  -0.67,
    -1.26,  -0.51,  -1.04,  -0.53,  -1.84,
    -2.42,  -4.52,  -4.76,  -6.33,  -6.76,
    -7.81,  -5.80,  -5.37,  -3.63,  -3.35,
    -1.75,  -1.88,  -0.61,  -0.90,   0.09,
    -0.32,   0.55,  -0.13,   0.70,  -0.06,
    0.49,  -0.20,   0.40,  -0.22,   0.36,
    -0.09,   0.58,   0.12,   0.75,   0.15,
    0.70,   0.17,   1.11,   0.89,   1.85,
    1.62,   2.54,   2.29,   3.20,   2.91,
    3.84,   3.53,   4.48,   4.15,   5.12,
    4.78,   5.75,   5.39,   6.31,   5.91,
    6.87,   6.33,   7.13,   6.61,   7.30,
    6.31,   6.27,   4.83,   4.49,   2.85,
    2.32,   0.58,  -0.11,  -0.98,   0.81,
    1.77,   3.37,   4.13,   5.60,   6.15,
    7.29,   7.35,   7.95,   7.67,   8.16,
    7.83,   8.31,   8.01,   8.53,   8.27};

  const G4int N = A - Z;
  //  const G4double x = (static_cast<G4double>(Z)*static_cast<G4double>(Z))/
    //    static_cast<G4double>(A);
  const G4double x = (G4double(Z)*G4double(Z))/G4double(A);
  G4double BF0 = 0.0;
  if (x <= 33.5) BF0 = 12.5 + 4.7*pow((33.5-x),0.75);
  else BF0 = 12.5 - 2.7*pow((x-33.5),2.0/3.0);


  // Determine which kind of nucleus is: even-even, odd-odd, even-odd, odd-even
  G4double D = 0.0;
  G4int I = 2*(Z/2); 
  G4int J = 2*(N/2);
  if (I < Z) D = 0.0;
  else D = -0.5;
  if (J < N) D += 1.0;

  if (Z > 130 || N > 200) return BF0 + D;
  else return BF0 + D - ShellCorr1[Z-1] - ShellCorr2[N-1];
}

