// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuonMinusCaptureAtRest physics process --------
//                   by Vladimir Ivanchenko
//                     E-mail: Vladimir.Ivantchenko@cern.ch
//                            April 2000
// **************************************************************
//    25.04.00  V.Ivanchenko replace arrays by STL vectors
//-----------------------------------------------------------------------------

#include "G4StopElementSelector.hh"
#include "g4std/vector"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// constructor
G4StopElementSelector::G4StopElementSelector()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// destructor
G4StopElementSelector::~G4StopElementSelector()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Element* G4StopElementSelector::GetElement(const G4Material* aMaterial)
{
  // Fermi-Teller Z-low of mu- capture and exceptions 
  // for halogens and oxigen.
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  G4int i;
  G4double Z;
  const G4int numberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  if(1 == numberOfElements) return (*theElementVector)(0);
    
  const G4double* theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();

  G4double sum = 0.0;
  for ( i=0; i < numberOfElements; i++ ) {

    Z = (*theElementVector)(i)->GetZ();

      // Halogens
    if( (9.0 == Z) || (17.0 == Z) || (35.0 == Z) || (53.0 == Z) || (85.0 == Z) ) {
      sum += 0.66 * Z * theAtomicNumberDensity[i] ; 

      // Oxigen
    } else if( 8.0 == Z ) {
      sum += 0.56 * Z * theAtomicNumberDensity[i] ; 

      // Others
    } else {
      sum +=        Z * theAtomicNumberDensity[i] ; 
    }
  }

  G4double random = G4UniformRand() * sum;
  sum = 0.0 ;
  i   = -1;

  // Selection of element
  do {
    i++;
    Z = (*theElementVector)(i)->GetZ();

      // Galogens
    if( (9.0 == Z) || (17.0 == Z) || (35.0 == Z) || (53.0 == Z) || (85.0 == Z) ) {
      sum += 0.66 * Z * theAtomicNumberDensity[i] ; 

      // Oxigen
    } else if( 8.0 == Z ) {
      sum += 0.56 * Z * theAtomicNumberDensity[i] ; 

      // Others
    } else {
      sum +=        Z * theAtomicNumberDensity[i] ; 
    }
  } while ( (sum < random) && (i < numberOfElements - 1) );

  return (*theElementVector)(i);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double  G4StopElementSelector::GetMuonCaptureRate(G4double Z, G4double A)
{
  // Initialized data

  //  static G4std::vector<G4double> zeff(100);
    static G4double zeff[100] = {
    1.,1.98,2.95,3.89,4.8,5.72,6.61,7.49,8.32,9.12,9.95,10.69,11.48,12.22,
    12.91,13.64,14.24,14.89,15.53,16.15,16.75,17.38,18.04,18.49,
    19.06,19.59,20.1,20.66,21.12,21.61,22.02,22.43,22.84,23.24,
    23.65,24.06,24.47,24.85,25.23,25.61,25.99,26.37,26.69,27.,
    27.32,27.63,27.95,28.2,28.42,28.64,28.79,29.03,29.27,29.51,
    29.75,29.99,30.2,30.36,30.53,30.69,30.85,31.01,31.18,31.34,
    31.48,31.62,31.76,31.9,32.05,32.19,32.33,32.47,32.61,32.76,
    32.94,33.11,33.29,33.46,33.64,33.81,34.21,34.18,34.,34.1,
    34.21,34.31,34.42,34.52,34.63,34.73,34.84,34.94,35.04,35.15,
    35.25,35.36,35.46,35.57,35.67,35.78 };

  // Mu- capture data from B.B.Balashov, G.Ya.Korenman, P.A.Eramgan
  // Atomizdat, 1978. (Experimental capture velocities)

  const size_t ListZE = 65;
  //  static G4std::vector<G4int> ListZExp[ListZE] = {
  static G4int ListZExp[ListZE] = {
      3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
     13, 14, 15, 16, 17, 18, 19, 20, 22, 23,
     24, 25, 26, 27, 28, 31, 32, 33, 34, 37,
     38, 39, 40, 41, 42, 45, 46, 47, 48, 49,
     50, 51, 52, 53, 55, 56, 57, 58, 59, 60,
     62, 64, 65, 67, 72, 73, 74, 80, 81, 82,
     83, 90, 92, 93};
  //  static G4std::vector<G4double> ListCaptureVel[ListZE] = {
  static G4double ListCaptureVel[ListZE] = {
     0.0057, 0.010, 0.0258, 0.0371, 0.0644,
     0.0974, 0.144, 0.250,  0.386,  0.479,
     0.700,  0.849, 1.119,  1.338,  1.40, 
     1.30,   1.98,  2.45,   2.60,   3.19,
     3.29,   3.91,  4.41,   4.96,   5.74,
     5.68,   5.53,  6.06,   5.69,   6.89,
     7.25,   7.89,  8.59,  10.40,   9.22,
    10.01,  10.00, 10.88,  10.62,  11.37,
    10.68,  10.49,  9.06,  11.20,  10.98,
    10.18,  10.71, 11.44,  13.45,  12.32,
    12.22,  12.09, 12.73,  12.95,  13.03,
    12.86,  13.13, 13.39,  12.74,  13.78,
    13.02,  13.26, 13.10,  14.00,  14.70};


  // Local variables
  G4double zeff2, xmu, a2ze, r1, r2;
  G4double lambda;

  // ==  Effective charges from Ford and Wills Nucl Phys 35(1962)295.
  // ==  Untabulated charges are interpolated.
  // ==  Mu capture lifetime (Goulard and Primakoff PRC10(1974)2034.

  G4int i = G4int(Z) - 1 ;
  if(i > 99) i = 99;

  const G4double b0a = -.03;
  const G4double b0b = -.25;
  const G4double b0c = 3.24;
  const G4double t1 = 875.e-10;
  r1 = zeff[i];
  zeff2 = r1 * r1;
  xmu = zeff2 * 2.663e-4;
  a2ze = 0.5 * A / Z;
  r2 = 1.0 - xmu;
  lambda = t1 * zeff2 * zeff2 * (r2 * r2) * (1.0 - (1.0 - xmu) * .75704) *
          (a2ze * b0a + 1.0 - (a2ze - 1.0) * b0b -
          (2.0 * (A - Z) /  Z  + abs(a2ze - 1.) ) * b0c / (A * 4.) );

  // == Mu capture data are taken if exist 
  for (G4int j = 0; j < ListZE; j++) {
    if( ListZExp[j] == i + 1) {
      lambda = ListCaptureVel[j] / microsecond;
      break;
    }
  }

  return lambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double  G4StopElementSelector::GetMuonDecayRate(G4double Z, G4double A)
{
  // Decay time on K-shell 
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.

  G4double lambda = 1.0 - 2.5 * Z * Z / (137.0*137.0);
  if( 0.5 > lambda ) lambda = 0.5;
  return lambda * 0.445 / microsecond; 
}




































































