//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4BraggModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// 04-12-02 Fix problem of G4DynamicParticle constructor (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)

// Class Description:
//
// Implementation of energy loss and delta-electron production by
// slow charged heavy particles

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4BraggModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BraggModel::G4BraggModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),
  particle(0),
  highKinEnergy(2.0*MeV),
  lowKinEnergy(0.0*MeV),
  iMolecula(0),
  protonMassAMU(1.007276),
  theZieglerFactor(eV*cm2*1.0e-15)
{
  if(p) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BraggModel::~G4BraggModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BraggModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()/eplus;
  chargeSquare = q*q;
  massRate     = mass/proton_mass_c2;
  highKinEnergy *= massRate;
  ratio = electron_mass_c2/mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::HighEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::LowEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::MinEnergyCut(const G4ParticleDefinition* p,
                                    const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4BraggModel::IsInCharge(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BraggModel::Initialise(const G4ParticleDefinition* p,
                              const G4DataVector&)
{
  if(!particle) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::ComputeDEDX(const G4Material* material,
                                   const G4ParticleDefinition* p,
                                         G4double kineticEnergy,
                                         G4double cutEnergy)
{

  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double dedx  = DEDX(material, kineticEnergy/massRate);

  if (cutEnergy < tmax) {

    G4double x   = cutEnergy/tmax;
    G4double tau = kineticEnergy/mass;
    G4double gam = tau + 1.0;
    G4double beta2 = 1. - 1./(gam*gam);
    //    G4double bg2   = tau * (tau+2.0);
    dedx += (log(x) + (1.0 - x)*beta2) * twopi_mc2_rcl2
          *  (material->GetElectronDensity())/beta2;
  }

  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;

  dedx *= chargeSquare;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::CrossSection(const G4Material* material,
                                    const G4ParticleDefinition* p,
                                          G4double kineticEnergy,
                                          G4double cutEnergy,
                                          G4double maxEnergy)
{

  G4double cross = 0.0;
  G4double tmax = G4std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(cutEnergy < tmax) {
    
    G4double x      = cutEnergy/tmax;
    G4double energy = kineticEnergy + mass;
    G4double gam    = energy/mass;
    G4double beta2  = 1. - 1./(gam*gam);
    cross = (1.0 - x*(1.0 - beta2*log(x)))/cutEnergy;
        
    cross *= twopi_mc2_rcl2*chargeSquare*material->GetElectronDensity()/beta2;
  }
  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax 
  //         << " cross= " << cross << G4endl; 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4BraggModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double xmin = tmin/tmax;
  G4double xmax = G4std::min(tmax, maxEnergy)/tmax;
  if(xmin >= xmax) return 0;

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy = kineticEnergy + mass;
  G4double beta2 = 1. - mass*mass/(energy*energy);
  G4double grej = 1.0 - beta2*xmin;
  G4double z, f;

  G4ThreeVector momentum = dp->GetMomentumDirection();

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    z = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    f = 1.0 - beta2 * z;

    if(f > grej) {
        G4cout << "G4BraggModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for x= " << z
               << G4endl;
    }

  } while( grej*G4UniformRand() >= f );

  G4double deltaKinEnergy = z * tmax;

  G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double totMomentum = sqrt(energy*energy - mass*mass);
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);
  G4double sint = sqrt(1.0 - cost*cost);

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  deltaDirection.rotateUz(momentum);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle();
  delta->SetDefinition(G4Electron::Electron());
  delta->SetKineticEnergy(deltaKinEnergy);
  delta->SetMomentumDirection(deltaDirection);

  return delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4BraggModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  G4DynamicParticle* delta = SampleSecondary(couple, dp, tmin, maxEnergy);
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4BraggModel::HasMaterial(const G4Material* material)
{
  G4String chFormula = material->GetChemicalFormula() ;
  G4String myFormula = G4String(" ") ;

  if (myFormula ==  chFormula ) {
    if(1 == (material->GetNumberOfElements())) return true;
    return false ;
  }

  // ICRU Report N49, 1993. Power's model for He.
  const size_t numberOfMolecula = 11 ;
  static G4String name[numberOfMolecula] = {
    "Al_2O_3",                 "CO_2",                      "CH_4",
    "(C_2H_4)_N-Polyethylene", "(C_2H_4)_N-Polypropylene",  "(C_8H_8)_N",
    "C_3H_8",                  "SiO_2",                     "H_2O",
    "H_2O-Gas",                "Graphite" } ;      
  
  // Special treatment for water in gas state
  const G4State theState = material->GetState() ;

  myFormula = G4String("H_2O");
  if( theState == kStateGas && myFormula == chFormula) {
    chFormula = G4String("H_2O-Gas");
  }
  
  // Search for the material in the table
  for (size_t i=0; i<numberOfMolecula; i++) {
      if (chFormula == name[i]) {
        SetMoleculaNumber(i) ;    
	return true ;
      }
  }
  return false ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::StoppingPower(const G4Material* material,
                                           G4double kineticEnergy) 
{
  G4double ionloss = 0.0 ;

  // pure material (normally not the case for this function)
  if(1 == (material->GetNumberOfElements())) {
    G4double z = material->GetZ() ;
    ionloss = ElectronicStoppingPower( z, kineticEnergy ) ;  

  } else if (iMolecula < 11) {
  
    // The data and the fit from: 
    // ICRU Report N49, 1993. Ziegler's model for protons.
    // Proton kinetic energy for parametrisation (keV/amu)  

    G4double T = kineticEnergy/(keV*protonMassAMU) ; 

     static G4double a[11][5] = {
   {1.187E+1, 1.343E+1, 1.069E+4, 7.723E+2, 2.153E-2}, 
   {7.802E+0, 8.814E+0, 8.303E+3, 7.446E+2, 7.966E-3}, 
   {7.294E+0, 8.284E+0, 5.010E+3, 4.544E+2, 8.153E-3}, 
   {8.646E+0, 9.800E+0, 7.066E+3, 4.581E+2, 9.383E-3}, 
   {1.286E+1, 1.462E+1, 5.625E+3, 2.621E+3, 3.512E-2}, 
   {3.229E+1, 3.696E+1, 8.918E+3, 3.244E+3, 1.273E-1}, 
   {1.604E+1, 1.825E+1, 6.967E+3, 2.307E+3, 3.775E-2}, 
   {8.049E+0, 9.099E+0, 9.257E+3, 3.846E+2, 1.007E-2}, 
   {4.015E+0, 4.542E+0, 3.955E+3, 4.847E+2, 7.904E-3}, 
   {4.571E+0, 5.173E+0, 4.346E+3, 4.779E+2, 8.572E-3}, 
   {2.631E+0, 2.601E+0, 1.701E+3, 1.279E+3, 1.638E-2} };
      

    if ( T < 10.0 ) {
      ionloss = a[iMolecula][0] * sqrt(T) ;
    
    } else if ( T < 10000.0 ) {
      G4double slow  = a[iMolecula][1] * pow(T, 0.45) ;
      G4double shigh = log( 1.0 + a[iMolecula][3]/T  
                     + a[iMolecula][4]*T ) * a[iMolecula][2]/T ;
      ionloss = slow*shigh / (slow + shigh) ;     
    } 

    if ( ionloss < 0.0) ionloss = 0.0 ;
    if ( 10 == iMolecula ) { 
      if (T < 100.0) {    
	ionloss *= (1.0+0.023+0.0066*log10(T));  
      }
      else if (T < 700.0) {   
	ionloss *=(1.0+0.089-0.0248*log10(T-99.));
      } 
      else if (T < 10000.0) {    
	ionloss *=(1.0+0.089-0.0248*log10(700.-99.));
      }
    }
  }
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::ElectronicStoppingPower(G4double z,
                                               G4double kineticEnergy) const
{
  G4double ionloss ;
  G4int i = G4int(z)-1 ;  // index of atom
  if(i < 0)  i = 0 ;
  if(i > 91) i = 91 ;
  
  // The data and the fit from: 
  // ICRU Report 49, 1993. Ziegler's type of parametrisations.
  // Proton kinetic energy for parametrisation (keV/amu)  

  G4double T = kineticEnergy/(keV*protonMassAMU) ; 
  
  static G4double a[92][5] = {
   {1.254E+0, 1.440E+0, 2.426E+2, 1.200E+4, 1.159E-1},
   {1.229E+0, 1.397E+0, 4.845E+2, 5.873E+3, 5.225E-2},
   {1.411E+0, 1.600E+0, 7.256E+2, 3.013E+3, 4.578E-2},
   {2.248E+0, 2.590E+0, 9.660E+2, 1.538E+2, 3.475E-2},
   {2.474E+0, 2.815E+0, 1.206E+3, 1.060E+3, 2.855E-2},
   {2.631E+0, 2.601E+0, 1.701E+3, 1.279E+3, 1.638E-2},
   {2.954E+0, 3.350E+0, 1.683E+3, 1.900E+3, 2.513E-2},
   {2.652E+0, 3.000E+0, 1.920E+3, 2.000E+3, 2.230E-2},
   {2.085E+0, 2.352E+0, 2.157E+3, 2.634E+3, 1.816E-2},
   {1.951E+0, 2.199E+0, 2.393E+3, 2.699E+3, 1.568E-2},
   {2.542E+0, 2.869E+0, 2.628E+3, 1.854E+3, 1.472E-2},
   {3.791E+0, 4.293E+0, 2.862E+3, 1.009E+3, 1.397E-2},
   {4.154E+0, 4.739E+0, 2.766E+3, 1.645E+2, 2.023E-2},
   {4.914E+0, 5.598E+0, 3.193E+3, 2.327E+2, 1.419E-2},
   {3.232E+0, 3.647E+0, 3.561E+3, 1.560E+3, 1.267E-2},
   {3.447E+0, 3.891E+0, 3.792E+3, 1.219E+3, 1.211E-2},
   {5.301E+0, 6.008E+0, 3.969E+3, 6.451E+2, 1.183E-2},
   {5.731E+0, 6.500E+0, 4.253E+3, 5.300E+2, 1.123E-2},
   {5.152E+0, 5.833E+0, 4.482E+3, 5.457E+2, 1.129E-2},
   {5.521E+0, 6.252E+0, 4.710E+3, 5.533E+2, 1.112E-2},
   {5.201E+0, 5.884E+0, 4.938E+3, 5.609E+2, 9.995E-3},
   {4.858E+0, 5.489E+0, 5.260E+3, 6.511E+2, 8.930E-3},
   {4.479E+0, 5.055E+0, 5.391E+3, 9.523E+2, 9.117E-3},
   {3.983E+0, 4.489E+0, 5.616E+3, 1.336E+3, 8.413E-3},
   {3.469E+0, 3.907E+0, 5.725E+3, 1.461E+3, 8.829E-3},
   {3.519E+0, 3.963E+0, 6.065E+3, 1.243E+3, 7.782E-3},
   {3.140E+0, 3.535E+0, 6.288E+3, 1.372E+3, 7.361E-3},
   {3.553E+0, 4.004E+0, 6.205E+3, 5.551E+2, 8.763E-3},
   {3.696E+0, 4.194E+0, 4.649E+3, 8.113E+1, 2.242E-2},
   {4.210E+0, 4.750E+0, 6.953E+3, 2.952E+2, 6.809E-3},
   {5.041E+0, 5.697E+0, 7.173E+3, 2.026E+2, 6.725E-3},
   {5.554E+0, 6.300E+0, 6.496E+3, 1.100E+2, 9.689E-3},
   {5.323E+0, 6.012E+0, 7.611E+3, 2.925E+2, 6.447E-3},
   {5.874E+0, 6.656E+0, 7.395E+3, 1.175E+2, 7.684E-3},
   {6.658E+0, 7.536E+0, 7.694E+3, 2.223E+2, 6.509E-3},
   {6.413E+0, 7.240E+0, 1.185E+4, 1.537E+2, 2.880E-3},
   {5.694E+0, 6.429E+0, 8.478E+3, 2.929E+2, 6.087E-3},
   {6.339E+0, 7.159E+0, 8.693E+3, 3.303E+2, 6.003E-3},
   {6.407E+0, 7.234E+0, 8.907E+3, 3.678E+2, 5.889E-3},
   {6.734E+0, 7.603E+0, 9.120E+3, 4.052E+2, 5.765E-3},
   {6.901E+0, 7.791E+0, 9.333E+3, 4.427E+2, 5.587E-3},
   {6.424E+0, 7.248E+0, 9.545E+3, 4.802E+2, 5.376E-3},
   {6.799E+0, 7.671E+0, 9.756E+3, 5.176E+2, 5.315E-3},
   {6.109E+0, 6.887E+0, 9.966E+3, 5.551E+2, 5.151E-3},
   {5.924E+0, 6.677E+0, 1.018E+4, 5.925E+2, 4.919E-3},
   {5.238E+0, 5.900E+0, 1.038E+4, 6.300E+2, 4.758E-3},
   {5.345E+0, 6.038E+0, 6.790E+3, 3.978E+2, 1.676E-2},
   {5.814E+0, 6.554E+0, 1.080E+4, 3.555E+2, 4.626E-3},
   {6.229E+0, 7.024E+0, 1.101E+4, 3.709E+2, 4.540E-3},
   {6.409E+0, 7.227E+0, 1.121E+4, 3.864E+2, 4.474E-3},
   {7.500E+0, 8.480E+0, 8.608E+3, 3.480E+2, 9.074E-3},
   {6.979E+0, 7.871E+0, 1.162E+4, 3.924E+2, 4.402E-3},
   {7.725E+0, 8.716E+0, 1.183E+4, 3.948E+2, 4.376E-3},
   {8.337E+0, 9.425E+0, 1.051E+4, 2.696E+2, 6.206E-3},
   {7.287E+0, 8.218E+0, 1.223E+4, 3.997E+2, 4.447E-3},
   {7.899E+0, 8.911E+0, 1.243E+4, 4.021E+2, 4.511E-3},
   {8.041E+0, 9.071E+0, 1.263E+4, 4.045E+2, 4.540E-3},
   {7.488E+0, 8.444E+0, 1.283E+4, 4.069E+2, 4.420E-3},
   {7.291E+0, 8.219E+0, 1.303E+4, 4.093E+2, 4.298E-3},
   {7.098E+0, 8.000E+0, 1.323E+4, 4.118E+2, 4.182E-3},
   {6.909E+0, 7.786E+0, 1.343E+4, 4.142E+2, 4.058E-3},
   {6.728E+0, 7.580E+0, 1.362E+4, 4.166E+2, 3.976E-3},
   {6.551E+0, 7.380E+0, 1.382E+4, 4.190E+2, 3.877E-3},
   {6.739E+0, 7.592E+0, 1.402E+4, 4.214E+2, 3.863E-3},
   {6.212E+0, 6.996E+0, 1.421E+4, 4.239E+2, 3.725E-3},
   {5.517E+0, 6.210E+0, 1.440E+4, 4.263E+2, 3.632E-3},
   {5.220E+0, 5.874E+0, 1.460E+4, 4.287E+2, 3.498E-3},
   {5.071E+0, 5.706E+0, 1.479E+4, 4.330E+2, 3.405E-3},
   {4.926E+0, 5.542E+0, 1.498E+4, 4.335E+2, 3.342E-3},
   {4.788E+0, 5.386E+0, 1.517E+4, 4.359E+2, 3.292E-3},
   {4.893E+0, 5.505E+0, 1.536E+4, 4.384E+2, 3.243E-3},
   {5.028E+0, 5.657E+0, 1.555E+4, 4.408E+2, 3.195E-3},
   {4.738E+0, 5.329E+0, 1.574E+4, 4.432E+2, 3.186E-3},
   {4.587E+0, 5.160E+0, 1.541E+4, 4.153E+2, 3.406E-3},
   {5.201E+0, 5.851E+0, 1.612E+4, 4.416E+2, 3.122E-3},
   {5.071E+0, 5.704E+0, 1.630E+4, 4.409E+2, 3.082E-3},
   {4.946E+0, 5.563E+0, 1.649E+4, 4.401E+2, 2.965E-3},
   {4.477E+0, 5.034E+0, 1.667E+4, 4.393E+2, 2.871E-3},
   {4.844E+0, 5.458E+0, 7.852E+3, 9.758E+2, 2.077E-2},
   {4.307E+0, 4.843E+0, 1.704E+4, 4.878E+2, 2.882E-3},
   {4.723E+0, 5.311E+0, 1.722E+4, 5.370E+2, 2.913E-3},
   {5.319E+0, 5.982E+0, 1.740E+4, 5.863E+2, 2.871E-3},
   {5.956E+0, 6.700E+0, 1.780E+4, 6.770E+2, 2.660E-3},
   {6.158E+0, 6.928E+0, 1.777E+4, 5.863E+2, 2.812E-3},
   {6.203E+0, 6.979E+0, 1.795E+4, 5.863E+2, 2.776E-3},
   {6.181E+0, 6.954E+0, 1.812E+4, 5.863E+2, 2.748E-3},
   {6.949E+0, 7.820E+0, 1.830E+4, 5.863E+2, 2.737E-3},
   {7.506E+0, 8.448E+0, 1.848E+4, 5.863E+2, 2.727E-3},
   {7.648E+0, 8.609E+0, 1.866E+4, 5.863E+2, 2.697E-3},
   {7.711E+0, 8.679E+0, 1.883E+4, 5.863E+2, 2.641E-3},
   {7.407E+0, 8.336E+0, 1.901E+4, 5.863E+2, 2.603E-3},
   {7.290E+0, 8.204E+0, 1.918E+4, 5.863E+2, 2.673E-3}
  };

  G4double fac = 1.0 ;

    // Carbon specific case for E < 40 keV
  if ( T < 40.0 && 5 == i) {
    fac = sqrt(T/40.0) ;
    T = 40.0 ;  

    // Free electron gas model
  } else if ( T < 10.0 ) { 
    fac = sqrt(T*0.1) ;
    T =10.0 ;  
  }

  // Main parametrisation
  G4double slow  = a[i][1] * pow(T, 0.45) ;
  G4double shigh = log( 1.0 + a[i][3]/T + a[i][4]*T ) * a[i][2]/T ;
  ionloss = slow*shigh*fac / (slow + shigh) ;     
  
  if ( ionloss < 0.0) ionloss = 0.0 ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::DEDX(const G4Material* material,
                                  G4double kineticEnergy) 
{
  G4double eloss = 0.0;
  const G4int numberOfElements = material->GetNumberOfElements();
  const G4double* theAtomicNumDensityVector =
                                 material->GetAtomicNumDensityVector();
  
  // pure material
  if(1 == numberOfElements) {

    G4double z = material->GetZ();
    eloss = ElectronicStoppingPower(z, kineticEnergy)
                               * (material->GetTotNbOfAtomsPerVolume());

  // compaund material with parametrisation
  } else if( HasMaterial(material) ) {

    eloss = StoppingPower(material, kineticEnergy)
                               * (material->GetTotNbOfAtomsPerVolume());
    G4int nAtoms = 0;
     
    const G4int* theAtomsVector = material->GetAtomsVector();
    for (G4int iel=0; iel<numberOfElements; iel++) {
      nAtoms += theAtomsVector[iel];
    }
    eloss /= nAtoms;

  // Experimental data exist only for kinetic energy 125 keV
  } else if( MolecIsInZiegler1988(material) ) { 

  // Cycle over elements - calculation based on Bragg's rule 
    G4double eloss125 = 0.0 ;
    const G4ElementVector* theElementVector =
                           material->GetElementVector();
  
    //  loop for the elements in the material
    for (G4int i=0; i<numberOfElements; i++) {
      const G4Element* element = (*theElementVector)[i] ;
      G4double z = element->GetZ() ;
      eloss    += ElectronicStoppingPower(z,kineticEnergy)
                                    * theAtomicNumDensityVector[i] ;
      eloss125 += ElectronicStoppingPower(z,125.0*keV)
                                    * theAtomicNumDensityVector[i] ;
    }      

    // Chemical factor is taken into account
    eloss *= ChemicalFactor(kineticEnergy, eloss125) ;
 
  // Brugg's rule calculation
  } else {
    const G4ElementVector* theElementVector =
                           material->GetElementVector() ;
  
    //  loop for the elements in the material
    for (G4int i=0; i<numberOfElements; i++)
    {
      const G4Element* element = (*theElementVector)[i] ;
      eloss   += ElectronicStoppingPower(element->GetZ(), kineticEnergy)
                                   * theAtomicNumDensityVector[i];
    }      
  }
  return eloss*theZieglerFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4BraggModel::MolecIsInZiegler1988(const G4Material* material) 
{
  // The list of molecules from
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
  
  G4String myFormula = G4String(" ") ;
  const G4String chFormula = material->GetChemicalFormula() ;
  if (myFormula == chFormula ) return false ;
  
  //  There are no evidence for difference of stopping power depended on
  //  phase of the compound except for water. The stopping power of the 
  //  water in gas phase can be predicted using Bragg's rule.
  //  
  //  No chemical factor for water-gas 
   
  myFormula = G4String("H_2O") ;
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && myFormula == chFormula) return false ;
    
  const size_t numberOfMolecula = 53 ;

  // The coffecient from Table.4 of Ziegler & Manoyan
  const G4double HeEff = 2.8735 ;    
  
  static G4String name[numberOfMolecula] = {
    "H_2O",      "C_2H_4O",    "C_3H_6O",  "C_2H_2",             "C_H_3OH",
    "C_2H_5OH",  "C_3H_7OH",   "C_3H_4",   "NH_3",               "C_14H_10",
    "C_6H_6",    "C_4H_10",    "C_4H_6",   "C_4H_8O",            "CCl_4",
    "CF_4",      "C_6H_8",     "C_6H_12",  "C_6H_10O",           "C_6H_10",
    "C_8H_16",   "C_5H_10",    "C_5H_8",   "C_3H_6-Cyclopropane","C_2H_4F_2",
    "C_2H_2F_2", "C_4H_8O_2",  "C_2H_6",   "C_2F_6",             "C_2H_6O",
    "C_3H_6O",   "C_4H_10O",   "C_2H_4",   "C_2H_4O",            "C_2H_4S",
    "SH_2",      "CH_4",       "CCLF_3",   "CCl_2F_2",           "CHCl_2F",
    "(CH_3)_2S", "N_2O",       "C_5H_10O", "C_8H_6",             "(CH_2)_N",
    "(C_3H_6)_N","(C_8H_8)_N", "C_3H_8",   "C_3H_6-Propylene",   "C_3H_6O",
    "C_3H_6S",   "C_4H_4S",    "C_7H_8"
  } ;
    
  static G4double expStopping[numberOfMolecula] = {
     66.1,  190.4, 258.7,  42.2, 141.5, 
    210.9,  279.6, 198.8,  31.0, 267.5,
    122.8,  311.4, 260.3, 328.9, 391.3,
    206.6,  374.0, 422.0, 432.0, 398.0,
    554.0,  353.0, 326.0,  74.6, 220.5,
    197.4,  362.0, 170.0, 330.5, 211.3,
    262.3,  349.6,  51.3, 187.0, 236.9,
    121.9,   35.8, 247.0, 292.6, 268.0,
    262.3,   49.0, 398.9, 444.0,  22.91,
     68.0,  155.0,  84.0,  74.2, 254.7,
    306.8,  324.4, 420.0
  } ;

  static G4double expCharge[numberOfMolecula] = {
    HeEff, HeEff, HeEff,   1.0, HeEff, 
    HeEff, HeEff, HeEff,   1.0,   1.0,
      1.0, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff,   1.0, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff,   1.0, HeEff, HeEff,
    HeEff,   1.0, HeEff, HeEff, HeEff,
    HeEff,   1.0, HeEff, HeEff,   1.0,
      1.0,   1.0,   1.0,   1.0, HeEff,
    HeEff, HeEff, HeEff
  } ;

  static G4double numberOfAtomsPerMolecula[numberOfMolecula] = {
    3.0,  7.0, 10.0,  4.0,  6.0,  
    9.0, 12.0,  7.0,  4.0, 24.0,
    12.0, 14.0, 10.0, 13.0,  5.0,
    5.0, 14.0, 18.0, 17.0, 17.0,
    24.0, 15.0, 13.0,  9.0,  8.0,
    6.0, 14.0,  8.0,  8.0,  9.0,
    10.0, 15.0,  6.0,  7.0,  7.0,
    3.0,  5.0,  5.0,  5.0,  5.0,
    9.0,  3.0, 16.0, 14.0,  3.0,
    9.0, 16.0, 11.0,  9.0, 10.0,
    10.0,  9.0, 15.0
  } ;

  // Search for the compaund in the table
  for (size_t i=0; i<numberOfMolecula; i++)
    { 
      if(chFormula == name[i]) { 
        G4double exp125 = expStopping[i] * 
	                  (material->GetTotNbOfAtomsPerVolume()) /
	                  (expCharge[i] * numberOfAtomsPerMolecula[i]) ;
        SetExpStopPower125(exp125);
        return true;
      }
    }
  
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BraggModel::ChemicalFactor(G4double kineticEnergy, 
                                      G4double eloss125) const
{
  // Approximation of Chemical Factor according to
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
  
  G4double gamma    = 1.0 + kineticEnergy/proton_mass_c2 ;    
  G4double gamma25  = 1.0 + 25.0*keV /proton_mass_c2 ;
  G4double gamma125 = 1.0 + 125.0*keV/proton_mass_c2 ;
  G4double beta     = sqrt(1.0 - 1.0/(gamma*gamma)) ;
  G4double beta25   = sqrt(1.0 - 1.0/(gamma25*gamma25)) ;
  G4double beta125  = sqrt(1.0 - 1.0/(gamma125*gamma125)) ;
  
  G4double factor = 1.0 + (expStopPower125/eloss125 - 1.0) *
                   (1.0 + exp( 1.48 * ( beta125/beta25 - 7.0 ) ) ) /
                   (1.0 + exp( 1.48 * ( beta/beta25    - 7.0 ) ) ) ;

  return factor ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

