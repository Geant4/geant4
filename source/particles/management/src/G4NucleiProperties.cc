// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NucleiProperties.cc,v 1.4 1999-08-20 14:26:54 larazb Exp $
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
// Added Shell-Pairing corrections to the Cameron mass 
// excess formula by V.Lara (9 May 99)


#include "G4NucleiProperties.hh"



G4double  G4NucleiProperties::AtomicMass(G4double Z, G4double A)
{
  const G4double hydrogen_mass_excess = G4NucleiPropertiesTable::GetMassExcess(1,1);
  const G4double neutron_mass_excess =  G4NucleiPropertiesTable::GetMassExcess(0,1);

  G4double mass =
      (A-Z)*neutron_mass_excess + Z*hydrogen_mass_excess - BindingEnergy(A,Z) + A*amu_c2;

  return mass;
}

G4double  G4NucleiProperties::BindingEnergy(G4double A, G4double Z)
{ 
  //
  // Weitzsaecker's Mass formula
  //
  G4int Npairing = G4int(A-Z)%2;                  // pairing
  G4int Zpairing = G4int(Z)%2;
  G4double binding =
      - 15.67*A                           // nuclear volume
      + 17.23*pow(A,2./3.)                // surface energy
      + 93.15*((A/2.-Z)*(A/2.-Z))/A       // asymmetry
      + 0.6984523*Z*Z*pow(A,-1./3.);      // coulomb
  if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / sqrt(A);  // pairing

  return -binding*MeV;
}


G4double G4NucleiProperties::GetNuclearMass(const G4double A, const G4double Z)
{
	if (A < 1 || Z < 0 || Z > A) {
		G4cout << "G4NucleiProperties::GetNuclearMass: Wrong values for A = " << A 
				 << " and Z = " << Z << endl;
		return 0.0;
	} else {
	 	G4ParticleDefinition * nucleus = 0;
		if ( (Z<=2) ) {
			if ( (Z==1)&&(A==1) ) {
				nucleus = G4ParticleTable::GetParticleTable()->FindParticle("proton"); // proton 
			} else if ( (Z==0)&&(A==1) ) {
			  	nucleus = G4ParticleTable::GetParticleTable()->FindParticle("neutron"); // neutron 
    		} else if ( (Z==1)&&(A==2) ) {
	  			nucleus = G4ParticleTable::GetParticleTable()->FindParticle("deuteron"); // deuteron 
    		} else if ( (Z==1)&&(A==3) ) {
	  			nucleus = G4ParticleTable::GetParticleTable()->FindParticle("triton"); // triton 
    		} else if ( (Z==2)&&(A==4) ) {
	  			nucleus = G4ParticleTable::GetParticleTable()->FindParticle("alpha"); // alpha 
    		} else if ( (Z==2)&&(A==3) ) {
	  			nucleus = G4ParticleTable::GetParticleTable()->FindParticle("He3"); // He3 
			}
  		}
		if (nucleus!=0) {
			return nucleus->GetPDGMass();
  		}else {
			return AtomicMass(A,Z) - Z*electron_mass_c2 + 1.433e-5*MeV*pow(Z,2.39);
		}
	}
}


// G4double G4NucleiProperties::CameronMassExcess(const G4int A, const G4int Z)
// {
// 	const G4double alpha = -17.0354*MeV;
// 	const G4double beta = -31.4506*MeV;
// 	const G4double phi = 44.2355*MeV;
// 	const G4double gamma = 25.8357*MeV;
//   
// 	const G4double A13 = pow(G4double(A),1.0/3.0);
// 	const G4double A23 = A13*A13;
// 	const G4double A43 = A23*A23;
// 	const G4double Z43 = pow(G4double(Z),4.0/3.0);
// 	G4double D = (G4double(A) - 2.0*G4double(Z))/G4double(A);
// 	D *= D; // D = pow((A-2Z)/A,2)
//   
//   
// 	// Surface term
// 	G4double SurfaceEnergy = (gamma - phi*D)*(1.0 - 0.62025/A23)*(1.0 - 0.62025/A23)*A23;
// 
// 	// Coulomb term
// 	G4double CoulombEnergy = 0.779*MeV*(G4double(Z*(Z-1))/A13)*(1.0-1.5849/A23+1.2273/A+1.5772/A43);
// 
// 	// Exchange term
// 	G4double ExchangeEnergy = -0.4323*MeV*(Z43/A13)*(1.0-0.57811/A13-0.14518/A23+0.49597/A);
//  
// 	// Volume term
// 	G4double VolumeEnergy = G4double(A)*(alpha-beta*D);
// 
// 	// Shell+Pairing corrections for protons
// 	G4double SPcorrectionZ;
// 	if (Z <= TableSize) SPcorrectionZ = SPZTable[Z-1]*MeV;
// 	else SPcorrectionZ = 0.0;
// 	
// 	// Shell+Pairing corrections for protons
// 	G4int N = A - Z;
// 	G4double SPcorrectionN;
// 	if (N <= TableSize) SPcorrectionN = SPNTable[N-1]*MeV;
// 	else SPcorrectionN = 0.0;
// 
// 
// 	// Mass Excess
// 	// First two terms give the mass excess of the neutrons and protons in the nucleus
// 	// (see Cameron, Canadian Journal of Physics, 35, 1957 page 1022) 
// 	G4double MassExcess = 8.367*MeV*A - 0.783*MeV*Z + 
// 								 SurfaceEnergy + CoulombEnergy + ExchangeEnergy + VolumeEnergy +
// 								 SPcorrectionZ + SPcorrectionN;
// 			 
// 	return MassExcess;
// }


// S(Z)+P(Z) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
//                   or Delta M(Z) from Tab. 97 of book [1]
// const G4double G4NucleiProperties::SPZTable[TableSize] = {
//   20.80, 15.80, 21.00, 16.80, 19.80, 16.50, 18.80, 16.50, 18.50, 17.20, //   1 - 10
//   18.26, 15.05, 16.01, 12.04, 13.27, 11.09, 12.17, 10.26, 11.04,  8.41, //  11 - 20 
//    9.79,  7.36,  8.15,  5.63,  5.88,  3.17,  3.32,   .82,  1.83,   .97, //  21 - 30
//    2.33,  1.27,  2.92,  1.61,  2.91,  1.35,  2.40,   .89,  1.74,   .36, //  31
//    0.95, -0.65, -0.04, -1.73, -0.96, -2.87, -2.05, -4.05, -3.40, -5.72, //  41
//   -3.75, -4.13, -2.42, -2.85, -1.01, -1.33,  0.54, -0.02,  1.74,  0.75, //  51
//    2.24,  1.00,  1.98,  0.79,  1.54,  0.39,  1.08,  0.00,  0.78, -0.35, //  61
//    0.58, -0.55,  0.59, -0.61,  0.59, -0.35,  0.32, -0.96, -0.52, -2.08, //  71
//   -2.46, -3.64, -1.55, -0.96,  0.97,  0.88,  2.37,  1.75,  2.72,  1.90, //  81
//    2.55,  1.46,  1.93,  0.86,  1.17,  0.08,  0.39, -0.76, -0.39, -1.51, //  91 - 100
//   -1.17, -2.36, -1.95, -3.06, -2.62, -3.55, -2.95, -3.75, -3.07, -3.79, // 101 - 110
//   -3.06, -3.77, -3.05, -3.78, -3.12, -3.90, -3.35, -4.24, -3.86, -4.92, // 111 - 120
//   -5.06, -6.77, -7.41, -9.18,-10.16,-11.12, -9.76, -9.23, -7.96, -7.65, // 121 - 130
//   // --------- from this point there are not tabulated values -----------------------
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 131 - 140
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 141 - 150
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 151
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 161
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 171
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, // 181
//    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00  // 191 - 200
// };

// S(N)+P(N) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
//                   or Delta M(N) from Tab. 97 of book [1]
// const G4double G4NucleiProperties::SPNTable[TableSize] = {
//   -8.40,-12.90, -8.00, 11.90, -9.20,-12.50,-10.80,-13.60,-11.20,-12.20, //   1 - 10
//  -12.81,-15.40,-13.07,-15.80,-13.81,-14.98,-12.63,-13.76,-11.37,-12.38, //  11 - 20
//   -9.23, -9.65, -7.64, -9.17, -8.05, -9.72, -8.87,-10.76, -8.64, -8.89, //  21 - 30
//   -6.60, -7.13, -4.77, -5.33, -3.06, -3.79, -1.72, -2.79, -0.93, -2.19, //  31
//   -0.52, -1.90, -0.45, -2.20, -1.22, -3.07, -2.42, -4.37, -3.94, -6.08, //  41
//   -4.49, -4.50, -3.14, -2.93, -1.04, -1.36,  0.69,  0.21,  2.11,  1.33, //  51
//    3.29,  2.46,  4.30,  3.32,  4.79,  3.62,  4.97,  3.64,  4.63,  3.07, //  61
//    4.06,  2.49,  3.30,  1.46,  2.06,  0.51,  0.74, -1.18, -1.26, -3.54, //  71
//   -3.97, -5.26, -4.18, -3.71, -2.10, -1.70, -0.08, -0.18,  0.94,  0.27, //  81
//    1.13,  0.08,  0.91, -0.31,  0.49, -0.78,  0.08, -1.15, -0.23, -1.41, //  91 - 100
//   -0.42, -1.55, -0.55, -1.66, -0.66, -1.73, -0.75, -1.74, -0.78, -1.69, // 101 - 110
//   -0.78, -1.60, -0.75, -1.46, -0.67, -1.26, -0.51, -1.04, -0.53, -1.84, // 111 - 120
//   -2.42, -4.52, -4.76, -6.33, -6.76, -7.81, -5.80, -5.37, -3.63, -3.35, // 121 - 130
//   -1.75, -1.88, -0.61, -0.90,  0.09, -0.32,  0.55, -0.13,  0.70, -0.06, // 131 - 140
//    0.49, -0.20,  0.40, -0.22,  0.36, -0.09,  0.58,  0.12,  0.75,  0.15, // 141 - 150
//    0.70,  0.17,  1.11,  0.89,  1.85,  1.62,  2.54,  2.29,  3.20,  2.91, // 151
//    3.84,  3.53,  4.48,  4.15,  5.12,  4.78,  5.75,  5.39,  6.31,  5.91, // 161
//    6.87,  6.33,  7.13,  6.61,  7.30,  6.31,  6.27,  4.83,  4.49,  2.85, // 171
//    2.32,  0.58, -0.11, -0.98,  0.81,  1.77,  3.37,  4.13,  5.60,  6.15, // 181
//    7.29,  7.35,  7.95,  7.67,  8.16,  7.83,  8.31,  8.01,  8.53,  8.27  // 191 - 200
// };


