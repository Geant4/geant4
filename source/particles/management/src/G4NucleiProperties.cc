// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NucleiProperties.cc,v 1.2 1999-04-13 08:00:22 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
// ------------------------------------------------------------
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
// Migrate into particles category by H.Kurashige (17 Nov. 98)


#include "G4NucleiProperties.hh"

G4double G4NucleiProperties::CameronMassExcess(const G4int A, const G4int Z)
{
  const G4double A03 = pow(A,1.0/3.0);
  const G4double A06 = A03*A03;
  const G4double A13 = A06*A06;
  const G4double Z13 = pow(Z,1.0+1.0/3.0);
  const G4double D = (A - 2.0*Z)/A;
  const G4double DR = 1.0 - 0.32052/A06;
  
  // Surface term
  G4double SurfaceEnergy = (25.8357-44.2355*D*D)*DR*DR*A06;
  // Coulomb term
  G4double CoulombEnergy = 0.779*(Z*(Z-1)/A03)*(1.0-1.5849/A06+1.2273/A+1.5772/A13);
  // Exchange term
  G4double ExchangeEnergy = -0.4323*(Z13/A03)*(1.0-0.57811/A03-0.14518/A06+0.49597/A);
  // Volume term
  G4double VolumeEnergy = A*(-17.0354+31.4506*D*D);
  // Compute Mass Excess
  // Neutron mass 8.07169 MeV Proton mass 7.2892 MeV
  return (8.07169*A - 0.7892*Z + SurfaceEnergy + CoulombEnergy + ExchangeEnergy + VolumeEnergy)*MeV;
}


G4double  G4NucleiProperties::AtomicMass(G4double Z, G4double A)
{
  // derived from original FORTRAN code ATOMAS by H. Fesefeldt (2-Dec-1986)
  //
  // Computes atomic mass in MeV
  // units for A example:  A = material->GetA()/(g/mole);
  //
  // Note:  can't just use aEff and zEff since the Nuclear Reaction
  //        function needs to calculate atomic mass for various values of A and Z
  
  G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4ParticleDefinition* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  if ((proton == 0)||(neutron == 0)||(electron == 0)) {
    G4Exception("G4NucleiProperties: G4Proton or G4Neutron is not defined !!"); 
  }
  const G4double proton_mass = proton->GetPDGMass();
  const G4double neutron_mass =  neutron->GetPDGMass();
  const G4double electron_mass = electron->GetPDGMass();
  //
  // Weitzsaecker's Mass formula
  //
  G4int nNeutron = A-Z;
  G4int ipp = G4int(nNeutron)%2;            // pairing
  G4int izz = G4int(Z)%2;
  G4double mass =
      nNeutron*neutron_mass + Z*proton_mass 
      - 15.67*double(A)                                  // nuclear volume
      + 17.23*pow(double(A),2./3.)                       // surface energy
      + 93.15*(double(A/2.-Z)*double(A/2.-Z))/double(A)  // asymmetry
      + 0.6984523*double(Z*Z)*pow(double(A),-1./3.)      // coulomb
      + Z*electron_mass;                                 // electrons mass 
  if( ipp == izz ) mass += (ipp+izz-1) * 12.0 / sqrt(double(A));  // pairing

  return mass*MeV;
}



// S(Z)+P(Z) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
//                   or Delta M(Z) from Tab. 97 of book [1]
const G4double G4NucleiProperties::daTZ[G4NucleiProperties::NTZ] = {
  20.80, 15.80, 21.00, 16.80, 19.80, 16.50, 18.80, 16.50, 18.50, 17.20, // 1
  18.26, 15.05, 16.01, 12.04, 13.27, 11.09, 12.17, 10.26, 11.04,  8.41, // 2
   9.79,  7.36,  8.15,  5.63,  5.88,  3.17,  3.32,   .82,  1.83,   .97, // 3
   2.33,  1.27,  2.92,  1.61,  2.91,  1.35,  2.40,   .89,  1.74,   .36, // 4
   0.95, -0.65, -0.04, -1.73, -0.96, -2.87, -2.05, -4.05, -3.40, -5.72, // 5
  -3.75, -4.13, -2.42, -2.85, -1.01, -1.33,  0.54, -0.02,  1.74,  0.75, // 6
   2.24,  1.00,  1.98,  0.79,  1.54,  0.39,  1.08,  0.00,  0.78, -0.35, // 7
   0.58, -0.55,  0.59, -0.61,  0.59, -0.35,  0.32, -0.96, -0.52, -2.08, // 8
  -2.46, -3.64, -1.55, -0.96,  0.97,  0.88,  2.37,  1.75,  2.72,  1.90, // 9
   2.55,  1.46,  1.93,  0.86,  1.17,  0.08,  0.39, -0.76, -0.39, -1.51, // 0
  -1.17, -2.36, -1.95, -3.06, -2.62, -3.55, -2.95, -3.75, -3.07, -3.79, // 1
  -3.06, -3.77, -3.05, -3.78, -3.12, -3.90, -3.35, -4.24, -3.86, -4.92, // 2
  -5.06, -6.77, -7.41, -9.18,-10.16,-11.12, -9.76, -9.23, -7.96, -7.65  // 3
};

// S(N)+P(N) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
//                   or Delta M(N) from Tab. 97 of book [1]
const G4double G4NucleiProperties::daTAZ[G4NucleiProperties::NTAZ] = {
  -8.40,-12.90, -8.00, 11.90, -9.20,-12.50,-10.80,-13.60,-11.20,-12.20, // 1
 -12.81,-15.40,-13.07,-15.80,-13.81,-14.98,-12.63,-13.76,-11.37,-12.38, // 2
  -9.23, -9.65, -7.64, -9.17, -8.05, -9.72, -8.87,-10.76, -8.64, -8.89, // 3
  -6.60, -7.13, -4.77, -5.33, -3.06, -3.79, -1.72, -2.79, -0.93, -2.19, // 4
  -0.52, -1.90, -0.45, -2.20, -1.22, -3.07, -2.42, -4.37, -3.94, -6.08, // 5
  -4.49, -4.50, -3.14, -2.93, -1.04, -1.36,  0.69,  0.21,  2.11,  1.33, // 6
   3.29,  2.46,  4.30,  3.32,  4.79,  3.62,  4.97,  3.64,  4.63,  3.07, // 7
   4.06,  2.49,  3.30,  1.46,  2.06,  0.51,  0.74, -1.18, -1.26, -3.54, // 8
  -3.97, -5.26, -4.18, -3.71, -2.10, -1.70, -0.08, -0.18,  0.94,  0.27, // 9
   1.13,  0.08,  0.91, -0.31,  0.49, -0.78,  0.08, -1.15, -0.23, -1.41, // 0
  -0.42, -1.55, -0.55, -1.66, -0.66, -1.73, -0.75, -1.74, -0.78, -1.69, // 1
  -0.78, -1.60, -0.75, -1.46, -0.67, -1.26, -0.51, -1.04, -0.53, -1.84, // 2
  -2.42, -4.52, -4.76, -6.33, -6.76, -7.81, -5.80, -5.37, -3.63, -3.35, // 3
  -1.75, -1.88, -0.61, -0.90,  0.09, -0.32,  0.55, -0.13,  0.70, -0.06, // 4
   0.49, -0.20,  0.40, -0.22,  0.36, -0.09,  0.58,  0.12,  0.75,  0.15, // 5
   0.70,  0.17,  1.11,  0.89,  1.85,  1.62,  2.54,  2.29,  3.20,  2.91, // 6
   3.84,  3.53,  4.48,  4.15,  5.12,  4.78,  5.75,  5.39,  6.31,  5.91, // 7
   6.87,  6.33,  7.13,  6.61,  7.30,  6.31,  6.27,  4.83,  4.49,  2.85, // 8
   2.32,  0.58, -0.11, -0.98,  0.81,  1.77,  3.37,  4.13,  5.60,  6.15, // 9
   7.29,  7.35,  7.95,  7.67,  8.16,  7.83,  8.31,  8.01,  8.53,  8.27  // 0
};




G4double G4NucleiProperties::PSCorrectedCameronMassExcess( const G4int A, const G4int Z )
{
  if( Z >= G4NucleiProperties::NTZ || A - Z >= G4NucleiProperties::NTAZ )
    G4Exception( "G4NucleiProperties::CameronMassExcess: Value parameters error!" );


  const G4double A03 = pow(A,1.0/3.0);
  const G4double A06 = A03*A03;
  const G4double A13 = A06*A06;
  const G4double Z13 = pow(Z,1.0+1.0/3.0);
  const G4double D = (A - 2.0*Z)/A;
  const G4double DR = 1.0 - 0.32052/A06;
  
  // Surface term
  G4double SurfaceEnergy = (25.8357 - 44.2355 * D * D) * DR*DR * A06;
  // Coulomb term
  G4double CoulombEnergy = 0.779* (Z*(Z-1)/A03)* (1.0 - 1.5849/A06 + 1.2273/A + 1.5772/A13);
  // Exchange term
  G4double ExchangeEnergy = -0.4323 * (Z13/A03) * (1.0 - 0.57811/A03 - 0.14518/A06 + 0.49597/A);
  // Volume term
  G4double VolumeEnergy = A* 8.367 + (31.4506/A) * (A -  2*Z)*(A -  2*Z) - 0.783* Z - 17.0354 * A;
  // Compute Mass Excess
  return (SurfaceEnergy + CoulombEnergy + ExchangeEnergy + VolumeEnergy + daTZ[Z] + daTAZ[A - Z])*MeV;

}

