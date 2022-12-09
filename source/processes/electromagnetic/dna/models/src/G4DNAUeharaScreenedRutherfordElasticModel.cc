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

#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

// #define UEHARA_VERBOSE

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAUeharaScreenedRutherfordElasticModel::
G4DNAUeharaScreenedRutherfordElasticModel(const G4ParticleDefinition*,
                                          const G4String& nam) : G4VEmModel(nam)
{
  // Energy limits of the models
  intermediateEnergyLimit = 200. * eV; 
  iLowEnergyLimit = 9.*eV;
  iHighEnergyLimit = 1.*MeV;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  fasterCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAUeharaScreenedRutherfordElasticModel::
Initialise(const G4ParticleDefinition* particle,
           const G4DataVector& /*cuts*/)
{
  if(isInitialised) { return; }
  if (verboseLevel > 3)
  {
  }
  
  if(particle->GetParticleName() != "e-")
  {
    G4Exception("*** WARNING: the G4DNAUeharaScreenedRutherfordElasticModel is "
                "not intented to be used with another particle than the electron",
                "",FatalException,"") ;
  }
  
  
  if( verboseLevel>1 )
  {
    G4cout << "G4DNAUeharaScreenedRutherfordElasticModel::Initialise()"
	   << G4endl;
    G4cout << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / MeV << " MeV"
	   << G4endl;
  }
  // Constants for final state by Brenner & Zaider
  // Note: the instantiation must be placed after if (isInitialised)
  
  betaCoeff=
  {
    7.51525,
    -0.41912,
    7.2017E-3,
    -4.646E-5,
    1.02897E-7};

  deltaCoeff=
  {
    2.9612,
    -0.26376,
    4.307E-3,
    -2.6895E-5,
    5.83505E-8};

  gamma035_10Coeff=
  {
    -1.7013,
    -1.48284,
    0.6331,
    -0.10911,
    8.358E-3,
    -2.388E-4};

  gamma10_100Coeff =
  {
    -3.32517,
    0.10996,
    -4.5255E-3,
    5.8372E-5,
    -2.4659E-7};

  gamma100_200Coeff=
  {
    2.4775E-2,
    -2.96264E-5,
    -1.20655E-7};
  
  // Initialize water density pointer
  fpWaterDensity =
   G4DNAMolecularMaterial::Instance()->
    GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
  
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAUeharaScreenedRutherfordElasticModel::
CrossSectionPerVolume(const G4Material* material,
                      const G4ParticleDefinition* /*particleDefinition*/,
                      G4double ekin,
                      G4double,
                      G4double)
{
#ifdef UEHARA_VERBOSE
  if (verboseLevel > 3)
  {
    G4cout
    << "Calling CrossSectionPerVolume() of G4DNAUeharaScreenedRutherfordElasticModel"
    << G4endl;
  }
#endif
  
  // Calculate total cross section for model
  
  G4double sigma = 0.;
  if(ekin < iLowEnergyLimit || ekin > iHighEnergyLimit) { return sigma; }

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];
  
  G4double z = 7.42; // FROM PMB 37 (1992) 1841-1858 p1842
  G4double n = ScreeningFactor(ekin,z);
  G4double crossSection = RutherfordCrossSection(ekin, z);
  sigma = pi * crossSection / (n * (n + 1.));
    
#ifdef UEHARA_VERBOSE
  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "=== G4DNAUeharaScreenedRutherfordElasticModel - XS INFO START"
           << G4endl;
    G4cout << "=== Kinetic energy(eV)=" << ekin/eV
           << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "=== Cross section per water molecule (cm^2)=" << sigma/cm/cm
           << G4endl;
    G4cout << "=== Cross section per water molecule (cm^-1)="
           << sigma*waterDensity/(1./cm) << G4endl;
    G4cout << "=== G4DNAUeharaScreenedRutherfordElasticModel - XS INFO END"
           << G4endl;
  }
#endif

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNAUeharaScreenedRutherfordElasticModel::RutherfordCrossSection(G4double k,
                                                                  G4double z)
{
  //
  //                               e^4         /      K + m_e c^2      \^2
  // sigma_Ruth(K) = Z (Z+1) -------------------- | --------------------- |
  //                          (4 pi epsilon_0)^2  \  K * (K + 2 m_e c^2)  /
  //
  // Where K is the electron non-relativistic kinetic energy
  //
  // NIM 155, pp. 145-156, 1978

  G4double length = (e_squared * (k + electron_mass_c2))
      / (4 * pi * epsilon0 * k * (k + 2 * electron_mass_c2));
  G4double cross = z * (z + 1) * length * length;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAUeharaScreenedRutherfordElasticModel::ScreeningFactor(G4double k,
                                                                    G4double z)
{
  // From Phys Med Biol 37 (1992) 1841-1858
  // Z=7.42 for water

  const G4double constK(1.7E-5);

  G4double beta2;
  beta2 = 1. - 1. / ((1. + k / electron_mass_c2) * (1. + k / electron_mass_c2));

  G4double etaC;
  if (k < 50 * keV)
    etaC = 1.198;
  else
    etaC = 1.13 + 3.76 * (z * z / (137 * 137 * beta2));

  G4double numerator = etaC * constK * std::pow(z, 2. / 3.);

  k /= electron_mass_c2;

  G4double denominator = k * (2 + k);

  G4double value = 0.;
  if (denominator > 0.)
    value = numerator / denominator; // form 3

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAUeharaScreenedRutherfordElasticModel::
SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                  const G4MaterialCutsCouple* /*couple*/,
                  const G4DynamicParticle* aDynamicElectron,
                  G4double,
                  G4double)
{
#ifdef UEHARA_VERBOSE
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling SampleSecondaries() of G4DNAUeharaScreenedRutherfordElasticModel"
        << G4endl;
  }
#endif

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  if(electronEnergy0 < iLowEnergyLimit || electronEnergy0 > iHighEnergyLimit)
    return;

  G4double cosTheta = 0.;

  if (electronEnergy0<intermediateEnergyLimit)
  {
#ifdef UEHARA_VERBOSE
    if (verboseLevel > 3)
      G4cout << "---> Using Brenner & Zaider model" << G4endl;
#endif
    cosTheta = BrennerZaiderRandomizeCosTheta(electronEnergy0);
  }
  else //if (electronEnergy0>=intermediateEnergyLimit)
  {
#ifdef UEHARA_VERBOSE
    if (verboseLevel > 3)
      G4cout << "---> Using Screened Rutherford model" << G4endl;
#endif
    G4double z = 7.42;  // FROM PMB 37 (1992) 1841-1858 p1842
    cosTheta = ScreenedRutherfordRandomizeCosTheta(electronEnergy0,z);
  }
  
  G4double phi = 2. * pi * G4UniformRand();
  
  G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
  G4ThreeVector xVers = zVers.orthogonal();
  G4ThreeVector yVers = zVers.cross(xVers);
  
  G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
  G4double yDir = xDir;
  xDir *= std::cos(phi);
  yDir *= std::sin(phi);
  
  G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));
  
  fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());
  
  fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAUeharaScreenedRutherfordElasticModel::
BrennerZaiderRandomizeCosTheta(G4double k)
{
  //  d sigma_el                         1                                 beta(K)
  // ------------ (K) ~ --------------------------------- + ---------------------------------
  //   d Omega           (1 + 2 gamma(K) - cos(theta))^2     (1 + 2 delta(K) + cos(theta))^2
  //
  // Maximum is < 1/(4 gamma(K)^2) + beta(K)/((2+2delta(K))^2)
  //
  // Phys. Med. Biol. 29 N.4 (1983) 443-447

  // gamma(K), beta(K) and delta(K) are polynomials with coefficients for energy measured in eV

  k /= eV;

  G4double beta = G4Exp(CalculatePolynomial(k, betaCoeff));
  G4double delta = G4Exp(CalculatePolynomial(k, deltaCoeff));
  G4double gamma;

  if (k > 100.)
  {
    gamma = CalculatePolynomial(k, gamma100_200Coeff);
    // Only in this case it is not the exponent of the polynomial
  } else
  {
    if (k > 10)
    {
      gamma = G4Exp(CalculatePolynomial(k, gamma10_100Coeff));
    }
    else
    {
      gamma = G4Exp(CalculatePolynomial(k, gamma035_10Coeff));
    }
  }

  // ***** Original method
  
  if (fasterCode == false)
  {
    G4double oneOverMax = 1.
    / (1. / (4. * gamma * gamma)
       + beta / ((2. + 2. * delta) * (2. + 2. * delta)));
    
    G4double cosTheta = 0.;
    G4double leftDenominator = 0.;
    G4double rightDenominator = 0.;
    G4double fCosTheta = 0.;
    
    do
    {
      cosTheta = 2. * G4UniformRand()- 1.;
      
      leftDenominator = (1. + 2.*gamma - cosTheta);
      rightDenominator = (1. + 2.*delta + cosTheta);
      if ( (leftDenominator * rightDenominator) != 0. )
      {
        fCosTheta = oneOverMax * (1./(leftDenominator*leftDenominator)
                                  + beta/(rightDenominator*rightDenominator));
      }
    }
    while (fCosTheta < G4UniformRand());
    
    return cosTheta;
  }

  // ***** Alternative method using cumulative probability

  else // if (fasterCode)
  {

   // 
   // modified by Shogo OKADA @ KEK, JP, 2016.2.27(Sat.)
   // 
   // An integral of differential cross-section formula shown above this member function
   // (integral variable: cos(theta), integral interval: [-1, x]) is as follows:
   // 
   //          1.0 + x                beta * (1 + x)
   // I = --------------------- + ----------------------   (1)
   //      (a - x) * (a + 1.0)      (b + x) * (b - 1.0)
   //
   // where a = 1.0 + 2.0 * gamma(K), b = 1.0 + 2.0 * delta(K)
   // 
   // Then, a cumulative probability (cp) is as follows:
   // 
   //  cp       1.0 + x                beta * (1 + x)      
   // ---- = --------------------- + ----------------------  (2)
   //  S      (a - x) * (a + 1.0)      (b + x) * (b - 1.0) 
   //
   // where 1/S is the integral of differnetical cross-section (1) on interval [-1, 1]
   // 
   //   1           2.0                     2.0 * beta
   //  --- = ----------------------- + -----------------------   (3)
   //   S     (a - 1.0) * (a + 1.0)     (b + 1.0) * (b - 1.0)
   //
   // x is calculated from the quadratic equation derived from (2) and (3):
   //
   // A * x^2 + B * x + C = 0
   //
   // where A, B, anc C are coefficients of the equation:
   //  A = S * {(b - 1.0) - beta * (a + 1.0)} + cp * (a + 1.0) * (b - 1.0),
   //  B = S * {(b - 1.0) * (b + 1.0) + beta * (a - 1.0) * (a + 1.0)} - cp * (a + 1.0) * (b - 1.0) * (a - b)
   //  C = S * {b * (b - 1.0) + beta * a * (a + 1.0)} - cp * (a + 1.0) * (b - 1.0) * ab
   //
   
   // sampling cumulative probability
   G4double cp = G4UniformRand();
  
   G4double a = 1.0 + 2.0 * gamma;
   G4double b = 1.0 + 2.0 * delta;
   G4double a1 = a - 1.0;
   G4double a2 = a + 1.0;
   G4double b1 = b - 1.0;
   G4double b2 = b + 1.0;
   G4double c1 = a - b;
   G4double c2 = a * b;

   G4double S = 2.0 / (a1 * a2) + 2.0 * beta / (b1 * b2); S = 1.0 / S;
 
   // coefficients for the quadratic equation
   G4double A = S * (b1 - beta * a2) + cp * a2 * b1;
   G4double B = S * (b1 * b2 + beta * a1 * a2) - cp * a2 * b1 * c1;
   G4double C = S * (b * b1 + beta * a * a2) - cp * a2 * b1 * c2;
  
   // calculate cos(theta)
   return (-1.0 * B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
 
   /*
   G4double cosTheta = -1;
   G4double cumul = 0;
   G4double value = 0;
   G4double leftDenominator = 0.;
   G4double rightDenominator = 0.;

   // Number of integration steps in the -1,1 range
   G4int iMax=200;

   G4double random = G4UniformRand();

   // Cumulate differential cross section
   for (G4int i=0; i<iMax; i++)
   {
   cosTheta = -1 + i*2./(iMax-1);
   leftDenominator = (1. + 2.*gamma - cosTheta);
   rightDenominator = (1. + 2.*delta + cosTheta);
   if ( (leftDenominator * rightDenominator) != 0. )
   {
   cumul = cumul + (1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator));
   }
   }

   // Select cosTheta
   for (G4int i=0; i<iMax; i++)
   {
   cosTheta = -1 + i*2./(iMax-1);
   leftDenominator = (1. + 2.*gamma - cosTheta);
   rightDenominator = (1. + 2.*delta + cosTheta);
   if (cumul !=0 && (leftDenominator * rightDenominator) != 0.)
   value = value + (1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator)) / cumul;
   if (random < value) break;
   }

   return cosTheta;
   */
  }
 
  //return 0.;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNAUeharaScreenedRutherfordElasticModel::
CalculatePolynomial(G4double k,
                    std::vector<G4double>& vec)
{
  // Sum_{i=0}^{size-1} vector_i k^i
  //
  // Phys. Med. Biol. 29 N.4 (1983) 443-447

  G4double result = 0.;
  size_t size = vec.size();

  while (size > 0)
  {
    size--;

    result *= k;
    result += vec[size];
  }

  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNAUeharaScreenedRutherfordElasticModel::
ScreenedRutherfordRandomizeCosTheta(G4double k,
                                    G4double z)
{

  //  d sigma_el                sigma_Ruth(K)
  // ------------ (K) ~ -----------------------------
  //   d Omega           (1 + 2 n(K) - cos(theta))^2
  //
  // We extract cos(theta) distributed as (1 + 2 n(K) - cos(theta))^-2
  //
  // Maximum is for theta=0: 1/(4 n(K)^2) (When n(K) is positive, that is always
  // satisfied within the validity of the process)
  //
  // Phys. Med. Biol. 45 (2000) 3171-3194
  
  // ***** Original method
  
  if (fasterCode == false)
  {
    G4double n = ScreeningFactor(k, z);
    
    G4double oneOverMax = (4. * n * n);
    
    G4double cosTheta = 0.;
    G4double fCosTheta;
    
    do
    {
      cosTheta = 2. * G4UniformRand()- 1.;
      fCosTheta = (1 + 2.*n - cosTheta);
      if (fCosTheta !=0.) fCosTheta = oneOverMax / (fCosTheta*fCosTheta);
    }
    while (fCosTheta < G4UniformRand());
    
    return cosTheta;
  }
  
  // ***** Alternative method using cumulative probability
  else // if (fasterCode)
  {
    
    //
    // modified by Shogo OKADA @ KEK, JP, 2016.2.27(Sat.)
    //
    // The cumulative probability (cp) is calculated by integrating
    // the differential cross-section fomula with cos(theta):
    //
    //         n(K) * (1.0 + cos(theta))
    //  cp = ---------------------------------
    //         1.0 + 2.0 * n(K) - cos(theta)
    //
    // Then, cos(theta) is as follows:
    //
    //               cp * (1.0 + 2.0 * n(K)) - n(K)
    // cos(theta) = --------------------------------
    //                       n(k) + cp
    //
    // where, K is kinetic energy, n(K) is screeing factor, and cp is cumulative probability
    //
    
    G4double n = ScreeningFactor(k, z);
    G4double cp = G4UniformRand();
    G4double numerator = cp * (1.0 + 2.0 * n) - n;
    G4double denominator = n + cp;
    return numerator / denominator;
    
    /*
     G4double cosTheta = -1;
     G4double cumul = 0;
     G4double value = 0;
     G4double n = ScreeningFactor(k, z);
     G4double fCosTheta;
     
     // Number of integration steps in the -1,1 range
     G4int iMax=200;
     
     G4double random = G4UniformRand();
     
     // Cumulate differential cross section
     for (G4int i=0; i<iMax; i++)
     {
     cosTheta = -1 + i*2./(iMax-1);
     fCosTheta = (1 + 2.*n - cosTheta);
     if (fCosTheta !=0.) cumul = cumul + 1./(fCosTheta*fCosTheta);
     }
     
     // Select cosTheta
     for (G4int i=0; i<iMax; i++)
     {
     cosTheta = -1 + i*2./(iMax-1);
     fCosTheta = (1 + 2.*n - cosTheta);
     if (cumul !=0.) value = value + (1./(fCosTheta*fCosTheta)) / cumul;
     if (random < value) break;
     }
     return cosTheta;
     */
  }
 
  //return 0.;
}
