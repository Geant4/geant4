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
// $Id: G4IonisParamMat.cc 106243 2017-09-26 01:56:43Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 09-07-98, data moved from G4Material, M.Maire
// 18-07-98, bug corrected in ComputeDensityEffect() for gas
// 16-01-01, bug corrected in ComputeDensityEffect() E100eV (L.Urban)
// 08-02-01, fShellCorrectionVector correctly handled (mma)
// 28-10-02, add setMeanExcitationEnergy (V.Ivanchenko)
// 06-09-04, factor 2 to shell correction term (V.Ivanchenko) 
// 10-05-05, add a missing coma in FindMeanExcitationEnergy() - Bug#746 (mma)
// 27-09-07, add computation of parameters for ions (V.Ivanchenko)
// 04-03-08, remove reference to G4NistManager. Add fBirks constant (mma)
// 30-10-09, add G4DensityEffectData class and density effect computation (VI)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4IonisParamMat.hh"
#include "G4Material.hh"
#include "G4DensityEffectData.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4DensityEffectData* G4IonisParamMat::fDensityData = nullptr;

#ifdef G4MULTITHREADED
  G4Mutex G4IonisParamMat::ionisMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamMat::G4IonisParamMat(const G4Material* material)
  : fMaterial(material)
{
  fBirks = 0.;
  fMeanEnergyPerIon = 0.0;
  twoln10 = 2.*G4Pow::GetInstance()->logZ(10);

  // minimal set of default parameters for density effect
  fCdensity = 0.0;
  fD0density = 0.0;
  fAdjustmentFactor = 1.0;
  if(fDensityData == nullptr) { fDensityData = new G4DensityEffectData(); }

  // compute parameters
  ComputeMeanParameters();
  ComputeDensityEffect();
  ComputeFluctModel();
  ComputeIonParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4IonisParamMat::G4IonisParamMat(__void__&)
  : fMaterial(nullptr), fShellCorrectionVector(nullptr)
{
  fMeanExcitationEnergy = 0.0;
  fLogMeanExcEnergy = 0.0;
  fTaul = 0.0;
  fCdensity = 0.0;
  fMdensity = 0.0;
  fAdensity = 0.0;
  fX0density = 0.0;
  fX1density = 0.0;
  fD0density = 0.0;
  fPlasmaEnergy = 0.0;
  fAdjustmentFactor = 0.0;
  fF1fluct = 0.0;          
  fF2fluct = 0.0;                       
  fEnergy1fluct = 0.0;
  fLogEnergy1fluct = 0.0;
  fEnergy2fluct = 0.0;
  fLogEnergy2fluct = 0.0;
  fEnergy0fluct = 0.0;
  fRateionexcfluct = 0.0;
  fZeff = 0.0;
  fFermiEnergy = 0.0;
  fLfactor = 0.0;
  fInvA23 = 0.0;
  fBirks = 0.0;
  fMeanEnergyPerIon = 0.0;
  twoln10 = 2.*G4Pow::GetInstance()->logZ(10);

  if(fDensityData == nullptr) { fDensityData = new G4DensityEffectData(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamMat::~G4IonisParamMat()
{
  delete [] fShellCorrectionVector; 
  delete fDensityData;
  fDensityData = nullptr;
  fShellCorrectionVector = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::ComputeMeanParameters()
{
  // compute mean excitation energy and shell correction vector
  fTaul = (*(fMaterial->GetElementVector()))[0]->GetIonisation()->GetTaul();

  fMeanExcitationEnergy = 0.;
  fLogMeanExcEnergy = 0.;

  size_t nElements = fMaterial->GetNumberOfElements();
  const G4ElementVector* elmVector = fMaterial->GetElementVector();
  const G4double* nAtomsPerVolume = fMaterial->GetVecNbOfAtomsPerVolume();
 
  fMeanExcitationEnergy = FindMeanExcitationEnergy(fMaterial);

  // Chemical formula defines mean excitation energy
  if(fMeanExcitationEnergy > 0.0) {
    fLogMeanExcEnergy = G4Log(fMeanExcitationEnergy);

    // Compute average 
  } else {
    for (size_t i=0; i < nElements; i++) {
      const G4Element* elm = (*elmVector)[i];
      fLogMeanExcEnergy += nAtomsPerVolume[i]*elm->GetZ()
	*G4Log(elm->GetIonisation()->GetMeanExcitationEnergy());
    }
    fLogMeanExcEnergy /= fMaterial->GetTotNbOfElectPerVolume();
    fMeanExcitationEnergy = G4Exp(fLogMeanExcEnergy);
  }

  fShellCorrectionVector = new G4double[3]; 

  for (G4int j=0; j<=2; j++)
  {
    fShellCorrectionVector[j] = 0.;

    for (size_t k=0; k<nElements; k++) {
      fShellCorrectionVector[j] += nAtomsPerVolume[k]
	*(((*elmVector)[k])->GetIonisation()->GetShellCorrectionVector())[j];
    }
    fShellCorrectionVector[j] *= 2.0/fMaterial->GetTotNbOfElectPerVolume();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4DensityEffectData* G4IonisParamMat::GetDensityEffectData()
{
  return fDensityData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
                    
void G4IonisParamMat::ComputeDensityEffect()
{
  G4State State = fMaterial->GetState();

  // Check if density effect data exist in the table
  // R.M. Sternheimer, Atomic Data and Nuclear Data Tables, 30: 261 (1984)
  G4int idx = fDensityData->GetIndex(fMaterial->GetName());
  G4int nelm= fMaterial->GetNumberOfElements();
  G4int Z0  = ((*(fMaterial->GetElementVector()))[0])->GetZasInt();

  // for simple non-NIST materials 
  G4double corr = 0.0;
  if(idx < 0 && 1 == nelm) {
    idx = fDensityData->GetElementIndex(Z0, fMaterial->GetState());

    // Correction for base material or for non-nominal density
    // Except cases of very different density defined in user code
    if(idx >= 0) {
      const G4Material* bmat = fMaterial->GetBaseMaterial();
      if(bmat) {
	corr = G4Log(bmat->GetDensity()/fMaterial->GetDensity());
      } else {
	G4double dens = G4NistManager::Instance()->GetNominalDensity(Z0);
	if(dens <= 0.0) { idx = -1; }
	else {
	  corr = G4Log(dens/fMaterial->GetDensity());
	}
      }
      // 1.0 is an arbitrary empirical limit 
      // parameterisation with very different density is not applicable
      if(std::abs(corr) > 1.0) { idx = -1; }
    }
  }

  //G4cout<<"DensityEffect for "<<fMaterial->GetName()<<"  "<< idx << G4endl; 

  if(idx >= 0) {

    // Take parameters for the density effect correction from
    // R.M. Sternheimer et al. Density Effect For The Ionization Loss 
    // of Charged Particles in Various Substances. 
    // Atom. Data Nucl. Data Tabl. 30 (1984) 261-271. 

    fCdensity  = fDensityData->GetCdensity(idx); 
    fMdensity  = fDensityData->GetMdensity(idx);
    fAdensity  = fDensityData->GetAdensity(idx);
    fX0density = fDensityData->GetX0density(idx);
    fX1density = fDensityData->GetX1density(idx);
    fD0density = fDensityData->GetDelta0density(idx);
    fPlasmaEnergy = fDensityData->GetPlasmaEnergy(idx);
    fAdjustmentFactor = fDensityData->GetAdjustmentFactor(idx);

    // parameter C is computed and not taken from Sternheimer tables
    //fCdensity = 1. + 2*G4Log(fMeanExcitationEnergy/fPlasmaEnergy);
    //G4cout << "IonisParamMat: " << fMaterial->GetName() 
    //	   << "  Cst= " << Cdensity << " C= " << fCdensity << G4endl;

    // correction on nominal density
    fCdensity  += corr;
    fX0density += corr/twoln10;
    fX1density += corr/twoln10;

  } else {

    static const G4double Cd2 = 4*pi*hbarc_squared*classic_electr_radius;
    fPlasmaEnergy = std::sqrt(Cd2*fMaterial->GetTotNbOfElectPerVolume());

    // Compute parameters for the density effect correction in DE/Dx formula.
    // The parametrization is from R.M. Sternheimer, Phys. Rev.B,3:3681 (1971)
    G4int icase;
    
    fCdensity = 1. + 2*G4Log(fMeanExcitationEnergy/fPlasmaEnergy);
    //
    // condensed materials
    //  
    if ((State == kStateSolid)||(State == kStateLiquid)) {

      static const G4double E100eV  = 100.*eV; 
      static const G4double ClimiS[] = {3.681 , 5.215 };
      static const G4double X0valS[] = {1.0   , 1.5   };
      static const G4double X1valS[] = {2.0   , 3.0   };
                                
      if(fMeanExcitationEnergy < E100eV) { icase = 0; }
      else                               { icase = 1; } 

      if(fCdensity < ClimiS[icase])    { fX0density = 0.2; }
      else { fX0density = 0.326*fCdensity - X0valS[icase]; }

      fX1density = X1valS[icase]; fMdensity = 3.0;
      
      //special: Hydrogen
      if (1 == nelm && 1 == Z0) {
         fX0density = 0.425; fX1density = 2.0; fMdensity = 5.949;
      }
    } else {
      //
      // gases
      //
      fMdensity = 3.;
      fX1density = 4.0;
      //static const G4double ClimiG[] = {10.,10.5,11.,11.5,12.25,13.804};
      //static const G4double X0valG[] = {1.6,1.7,1.8,1.9,2.0,2.0};
      //static const G4double X1valG[] = {4.0,4.0,4.0,4.0,4.0,5.0};

      if(fCdensity < 10.) {
	fX0density = 1.6; 
      } else if(fCdensity < 11.5) { 
	fX0density = 1.6 + 0.2*(fCdensity - 10.); 
      } else if(fCdensity < 12.25) { 
	fX0density = 1.9 + (fCdensity - 11.5)/7.5; 
      } else if(fCdensity < 13.804) { 
	fX0density = 2.0; 
	fX1density = 4.0 + (fCdensity - 12.25)/1.554;
      } else {
	fX0density = 0.326*fCdensity-2.5; fX1density = 5.0; 
      }
      
      //special: Hydrogen
      if (1 == nelm && 1 == Z0) {
         fX0density = 1.837; fX1density = 3.0; fMdensity = 4.754;
      }
      
      //special: Helium
      if (1 == nelm && 2 == Z0) {
         fX0density = 2.191; fX1density = 3.0; fMdensity = 3.297;
      }
    }
  }

  // change parameters if the gas is not in STP.
  // For the correction the density(STP) is needed. 
  // Density(STP) is calculated here : 
  
    
  if (State == kStateGas) { 
    G4double Density  = fMaterial->GetDensity();
    G4double Pressure = fMaterial->GetPressure();
    G4double Temp     = fMaterial->GetTemperature();
      
    G4double DensitySTP = Density*STP_Pressure*Temp/(Pressure*NTP_Temperature);

    G4double ParCorr = G4Log(Density/DensitySTP);
  
    fCdensity  -= ParCorr;
    fX0density -= ParCorr/twoln10;
    fX1density -= ParCorr/twoln10;
  }

  // fAdensity parameter can be fixed for not conductive materials 
  if(0.0 == fD0density) {
    G4double Xa = fCdensity/twoln10;
    fAdensity = twoln10*(Xa-fX0density)
      /std::pow((fX1density-fX0density),fMdensity);
  }
  /*  
  G4cout << "G4IonisParamMat: density effect data for <" 
         << fMaterial->GetName() 
	 << "> " << G4endl;
  G4cout << "Eplasma(eV)= " << fPlasmaEnergy/eV
	 << " rho= " << fAdjustmentFactor
	 << " -C= " << fCdensity 
	 << " x0= " << fX0density
	 << " x1= " << fX1density
	 << " a= " << fAdensity
	 << " m= " << fMdensity
	 << G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::ComputeFluctModel()
{
  // compute parameters for the energy loss fluctuation model
  // needs an 'effective Z' 
  G4double Zeff = 0.;
  for (size_t i=0;i<fMaterial->GetNumberOfElements();i++) {
     Zeff += (fMaterial->GetFractionVector())[i]
             *((*(fMaterial->GetElementVector()))[i]->GetZ());
  }
  if (Zeff > 2.) { fF2fluct = 2./Zeff; }
  else           { fF2fluct = 0.; }

  fF1fluct         = 1. - fF2fluct;
  fEnergy2fluct    = 10.*Zeff*Zeff*eV;
  fLogEnergy2fluct = G4Log(fEnergy2fluct);
  fLogEnergy1fluct = (fLogMeanExcEnergy - fF2fluct*fLogEnergy2fluct)
                     /fF1fluct;
  fEnergy1fluct    = G4Exp(fLogEnergy1fluct);
  fEnergy0fluct    = 10.*eV;
  fRateionexcfluct = 0.4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::ComputeIonParameters()
{
  // get elements in the actual material,
  const G4ElementVector* theElementVector = fMaterial->GetElementVector() ;
  const G4double* theAtomicNumDensityVector =
                         fMaterial->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = fMaterial->GetNumberOfElements() ;

  //  loop for the elements in the material
  //  to find out average values Z, vF, lF
  G4double z(0.0), vF(0.0), lF(0.0), norm(0.0), a23(0.0);

  G4Pow* g4pow = G4Pow::GetInstance();
  if( 1 == NumberOfElements ) {
    const G4Element* element = (*theElementVector)[0];
    z = element->GetZ();
    vF= element->GetIonisation()->GetFermiVelocity();
    lF= element->GetIonisation()->GetLFactor();
    a23 = 1.0/g4pow->A23(element->GetN());

  } else {
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        const G4Element* element = (*theElementVector)[iel];
        const G4double weight = theAtomicNumDensityVector[iel];
        norm += weight ;
        z    += element->GetZ() * weight;
        vF   += element->GetIonisation()->GetFermiVelocity() * weight;
        lF   += element->GetIonisation()->GetLFactor() * weight;
	a23  += weight/g4pow->A23(element->GetN());
      }
    z  /= norm;
    vF /= norm;
    lF /= norm;
    a23 /= norm;
  }  
  fZeff        = z;
  fLfactor     = lF;
  fFermiEnergy = 25.*keV*vF*vF;
  fInvA23      = a23;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::SetMeanExcitationEnergy(G4double value)
{
  if(value == fMeanExcitationEnergy || value <= 0.0) { return; }
  if (G4NistManager::Instance()->GetVerbose() > 1) {
    G4cout << "G4Material: Mean excitation energy is changed for "
           << fMaterial->GetName()
           << " Iold= " << fMeanExcitationEnergy/eV
           << "eV; Inew= " << value/eV << " eV;"
           << G4endl;
  }
  
  fMeanExcitationEnergy = value;

  // add corrections to density effect
  G4double newlog = G4Log(value);
  G4double corr = 2*(newlog - fLogMeanExcEnergy);
  fCdensity  += corr;
  fX0density += corr/twoln10;
  fX1density += corr/twoln10;

  // recompute parameters of fluctuation model
  fLogMeanExcEnergy = newlog;
  ComputeFluctModel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::SetDensityEffectParameters(G4double cd, G4double md, 
						 G4double ad, G4double x0, 
						 G4double x1, G4double d0)
{
  // no check on consistence of user parameters  
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&ionisMutex);
#endif
  fCdensity  = cd; 
  fMdensity  = md;
  fAdensity  = ad;
  fX0density = x0;
  fX1density = x1;
  fD0density = d0;
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&ionisMutex);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::SetDensityEffectParameters(const G4Material* bmat)
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&ionisMutex);
#endif
  const G4IonisParamMat* ipm = bmat->GetIonisation();
  fCdensity  = ipm->GetCdensity(); 
  fMdensity  = ipm->GetMdensity();
  fAdensity  = ipm->GetAdensity();
  fX0density = ipm->GetX0density();
  fX1density = ipm->GetX1density();
  fD0density = ipm->GetD0density();

  G4double corr = G4Log(bmat->GetDensity()/fMaterial->GetDensity());
  fCdensity  += corr;
  fX0density += corr/twoln10;
  fX1density += corr/twoln10;
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&ionisMutex);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double G4IonisParamMat::FindMeanExcitationEnergy(const G4Material* mat) const
{
  G4double res = 0.0;
  // data from density effect data
  if(fDensityData) {
    G4int idx = fDensityData->GetIndex(mat->GetName());
    if(idx >= 0) {
      res = fDensityData->GetMeanIonisationPotential(idx);
    }
  }

  // The data on mean excitation energy for compaunds
  // from "Stopping Powers for Electrons and Positrons"
  // ICRU Report N#37, 1984  (energy in eV)
  // this value overwrites Density effect data
  G4String chFormula = mat->GetChemicalFormula();
  if(chFormula != "") { 

    static const size_t numberOfMolecula = 54; 
    static const G4String name[numberOfMolecula] = {
      // gas 0 - 12
      "NH_3",       "C_4H_10",     "CO_2",      "C_2H_6",      "C_7H_16-Gas",
      // "G4_AMMONIA", "G4_BUTANE","G4_CARBON_DIOXIDE","G4_ETHANE", "G4_N-HEPTANE"
      "C_6H_14-Gas",   "CH_4",     "NO",        "N_2O",        "C_8H_18-Gas",      
      // "G4_N-HEXANE" , "G4_METHANE", "x", "G4_NITROUS_OXIDE", "G4_OCTANE"
      "C_5H_12-Gas",   "C_3H_8",   "H_2O-Gas",                        
      // "G4_N-PENTANE", "G4_PROPANE", "G4_WATER_VAPOR"

      // liquid 13 - 39
      "C_3H_6O",    "C_6H_5NH_2",  "C_6H_6",    "C_4H_9OH",    "CCl_4", 
      //"G4_ACETONE","G4_ANILINE","G4_BENZENE","G4_N-BUTYL_ALCOHOL","G4_CARBON_TETRACHLORIDE"  
      "C_6H_5Cl",   "CHCl_3",      "C_6H_12",   "C_6H_4Cl_2",  "C_4Cl_2H_8O",
      //"G4_CHLOROBENZENE","G4_CHLOROFORM","G4_CYCLOHEXANE","G4_1,2-DICHLOROBENZENE",
      //"G4_DICHLORODIETHYL_ETHER"
      "C_2Cl_2H_4", "(C_2H_5)_2O", "C_2H_5OH",  "C_3H_5(OH)_3","C_7H_16",
      //"G4_1,2-DICHLOROETHANE","G4_DIETHYL_ETHER","G4_ETHYL_ALCOHOL","G4_GLYCEROL","G4_N-HEPTANE"
      "C_6H_14",    "CH_3OH",      "C_6H_5NO_2","C_5H_12",     "C_3H_7OH",
      //"G4_N-HEXANE","G4_METHANOL","G4_NITROBENZENE","G4_N-PENTANE","G4_N-PROPYL_ALCOHOL",  
      "C_5H_5N",    "C_8H_8",      "C_2Cl_4",   "C_7H_8",      "C_2Cl_3H",
      //"G4_PYRIDINE","G4_POLYSTYRENE","G4_TETRACHLOROETHYLENE","G4_TOLUENE","G4_TRICHLOROETHYLENE"
      "H_2O",       "C_8H_10",             
      // "G4_WATER", "G4_XYLENE"

      // solid 40 - 53
      "C_5H_5N_5",  "C_5H_5N_5O",  "(C_6H_11NO)-nylon",  "C_25H_52",
      // "G4_ADENINE", "G4_GUANINE", "G4_NYLON-6-6", "G4_PARAFFIN"
      "(C_2H_4)-Polyethylene",     "(C_5H_8O_2)-Polymethil_Methacrylate",
      // "G4_ETHYLENE", "G4_PLEXIGLASS"
      "(C_8H_8)-Polystyrene",      "A-150-tissue",       "Al_2O_3",  "CaF_2",
      // "G4_POLYSTYRENE", "G4_A-150_TISSUE", "G4_ALUMINUM_OXIDE", "G4_CALCIUM_FLUORIDE"
      "LiF",        "Photo_Emulsion",  "(C_2F_4)-Teflon",  "SiO_2"
      // "G4_LITHIUM_FLUORIDE", "G4_PHOTO_EMULSION", "G4_TEFLON", "G4_SILICON_DIOXIDE"
    } ;
    
    static const G4double meanExcitation[numberOfMolecula] = {

      53.7,   48.3,  85.0,  45.4,  49.2,
      49.1,   41.7,  87.8,  84.9,  49.5,
      48.2,   47.1,  71.6,

      64.2,   66.2,  63.4,  59.9,  166.3,
      89.1,  156.0,  56.4, 106.5,  103.3, 
      111.9,   60.0,  62.9,  72.6,   54.4,  
      54.0,  67.6,   75.8,  53.6,   61.1,  
      66.2,  64.0,  159.2,  62.5,  148.1,  
      75.0,  61.8,

      71.4,  75.0,   63.9,  48.3,   57.4,
      74.0,  68.7,   65.1, 145.2,  166.,
      94.0, 331.0,   99.1, 139.2 
    };

    for(size_t i=0; i<numberOfMolecula; i++) {
      if(chFormula == name[i]) {
	res = meanExcitation[i]*eV;
	break;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

