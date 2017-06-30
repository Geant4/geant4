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
// ********************************************************************
// *********************************************************************
// |                                                                   |
// |      G4LowEPPolarizedComptonModel-- Geant4 Monash University      |
// |         polarised low energy Compton scattering model.            |
// |          J. M. C. Brown, Monash University, Australia             |
// |                                                                   |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | The following is a Geant4 class to simulate the process of        |
// | bound electron Compton scattering. General code structure is      |
// | based on G4LowEnergyCompton.cc and                                |
// | G4LivermorePolarizedComptonModel.cc.                              |
// | Algorithms for photon energy, and ejected Compton electron        |
// | direction taken from:                                             |
// |                                                                   |
// | J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin,    |
// | "A low energy bound atomic electron Compton scattering model      |
// |  for Geant4", NIMB, Vol. 338, 77-88, 2014.                        |
// |                                                                   |
// | The author acknowledges the work of the Geant4 collaboration      |
// | in developing the following algorithms that have been employed    |
// | or adapeted for the present software:                             |    
// |                                                                   |
// |  # sampling of photon scattering angle,                           |
// |  # target element selection in composite materials,               |
// |  # target shell selection in element,                             |
// |  # and sampling of bound electron momentum from Compton profiles. |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | History:                                                          |
// | --------                                                          |
// |                                                                   |
// | Jan. 2015 JMCB       - 1st Version based on G4LowEPPComptonModel  |
// | Feb. 2016 JMCB       - Geant4 10.2 FPE fix for bug 1676           |
// | Nov. 2016 JMCB       - Polarisation tracking fix in collaboration |
// |                        of Dr. Merlin Reynaard Kole,               |
// |                        University of Geneva                       |
// |                                                                   |
// *********************************************************************

#include "G4LowEPPolarizedComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//****************************************************************************

using namespace std;

G4int G4LowEPPolarizedComptonModel::maxZ = 99;
G4LPhysicsFreeVector* G4LowEPPolarizedComptonModel::data[] = {0};
G4ShellData*       G4LowEPPolarizedComptonModel::shellData = 0;
G4DopplerProfile*  G4LowEPPolarizedComptonModel::profileData = 0;

static const G4double ln10 = G4Log(10.);

G4LowEPPolarizedComptonModel::G4LowEPPolarizedComptonModel(const G4ParticleDefinition*,
                                                 const G4String& nam)
  : G4VEmModel(nam),isInitialised(false)
{
  verboseLevel=1 ;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(  verboseLevel>1 ) {
    G4cout << "Low energy photon Compton model is constructed " << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  fParticleChange = 0;
  fAtomDeexcitation = 0;
}

//****************************************************************************

G4LowEPPolarizedComptonModel::~G4LowEPPolarizedComptonModel()
{
  if(IsMaster()) {
    delete shellData;
    shellData = 0;
    delete profileData;
    profileData = 0;
  }
}

//****************************************************************************

void G4LowEPPolarizedComptonModel::Initialise(const G4ParticleDefinition* particle,
                                         const G4DataVector& cuts)
{
  if (verboseLevel > 1) {
    G4cout << "Calling G4LowEPPolarizedComptonModel::Initialise()" << G4endl;
  }

  // Initialise element selector

  if(IsMaster()) {

    // Access to elements

    char* path = getenv("G4LEDATA");

    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();

    for(G4int i=0; i<numOfCouples; ++i) {
      const G4Material* material =
        theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();

      for (G4int j=0; j<nelm; ++j) {
        G4int Z = G4lrint((*theElementVector)[j]->GetZ());
        if(Z < 1)        { Z = 1; }
        else if(Z > maxZ){ Z = maxZ; }

        if( (!data[Z]) ) { ReadData(Z, path); }
      }
    }

    // For Doppler broadening
    if(!shellData) {
      shellData = new G4ShellData();
      shellData->SetOccupancyData();
      G4String file = "/doppler/shell-doppler";
      shellData->LoadData(file);
    }
    if(!profileData) { profileData = new G4DopplerProfile(); }

    InitialiseElementSelectors(particle, cuts);
  }

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files" << G4endl;
  }

  if( verboseLevel>1 ) {
    G4cout << "G4LowEPPolarizedComptonModel is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / GeV << " GeV"
           << G4endl;
  }

  if(isInitialised) { return; }
  
  fParticleChange = GetParticleChangeForGamma();
  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();
  isInitialised = true;
}

//****************************************************************************

void G4LowEPPolarizedComptonModel::InitialiseLocal(const G4ParticleDefinition*,
                                              G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//****************************************************************************

void G4LowEPPolarizedComptonModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1)
  {
    G4cout << "G4LowEPPolarizedComptonModel::ReadData()"
           << G4endl;
  }
  if(data[Z]) { return; }
  const char* datadir = path;
  if(!datadir)
  {
    datadir = getenv("G4LEDATA");
    if(!datadir)
    {
      G4Exception("G4LowEPPolarizedComptonModel::ReadData()",
                  "em0006",FatalException,
                  "Environment variable G4LEDATA not defined");
      return;
    }
  }

  data[Z] = new G4LPhysicsFreeVector();

  // Activation of spline interpolation
  data[Z]->SetSpline(false);

  std::ostringstream ost;
  ost << datadir << "/livermore/comp/ce-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());

  if( !fin.is_open())
    {
      G4ExceptionDescription ed;
      ed << "G4LowEPPolarizedComptonModel data file <" << ost.str().c_str()
         << "> is not opened!" << G4endl;
    G4Exception("G4LowEPPolarizedComptonModel::ReadData()",
                "em0003",FatalException,
                ed,"G4LEDATA version should be G4EMLOW6.34 or later");
      return;
    } else {
      if(verboseLevel > 3) {
        G4cout << "File " << ost.str()
               << " is opened by G4LowEPPolarizedComptonModel" << G4endl;
      }
      data[Z]->Retrieve(fin, true);
      data[Z]->ScaleVector(MeV, MeV*barn);
    }
  fin.close();
}

//****************************************************************************


G4double 
G4LowEPPolarizedComptonModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                    G4double GammaEnergy,
                                                    G4double Z, G4double,
                                                    G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "G4LowEPPolarizedComptonModel::ComputeCrossSectionPerAtom()"
           << G4endl;
  }
  G4double cs = 0.0;

  if (GammaEnergy < LowEnergyLimit()) { return 0.0; }

  G4int intZ = G4lrint(Z);
  if(intZ < 1 || intZ > maxZ) { return cs; }

  G4LPhysicsFreeVector* pv = data[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv)
    {
      InitialiseForElement(0, intZ);
      pv = data[intZ];
      if(!pv) { return cs; }
    }

  G4int n = pv->GetVectorLength() - 1;
  G4double e1 = pv->Energy(0);
  G4double e2 = pv->Energy(n);

  if(GammaEnergy <= e1)      { cs = GammaEnergy/(e1*e1)*pv->Value(e1); }
  else if(GammaEnergy <= e2) { cs = pv->Value(GammaEnergy)/GammaEnergy; }
  else if(GammaEnergy > e2)  { cs = pv->Value(e2)/GammaEnergy; }

  return cs;
}

//****************************************************************************

void G4LowEPPolarizedComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* aDynamicGamma,
                                                G4double, G4double)
{

    //Determine number of digits (in decimal base) that G4double can accurately represent
     G4double g4d_order = G4double(numeric_limits<G4double>::digits10);
     G4double g4d_limit = std::pow(10.,-g4d_order);

  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // then accepted or rejected depending on the Scattering Function multiplied
  // by factor from Klein - Nishina formula.
  // Expression of the angular distribution as Klein Nishina
  // angular and energy distribution and Scattering fuctions is taken from
  // D. E. Cullen "A simple model of photon transport" Nucl. Instr. Meth.
  // Phys. Res. B 101 (1995). Method of sampling with form factors is different
  // data are interpolated while in the article they are fitted.
  // Reference to the article is from J. Stepanek New Photon, Positron
  // and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10
  // TeV (draft).
  // The random number techniques of Butcher & Messel are used
  // (Nucl Phys 20(1960),15).


  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy()/MeV;

  if (verboseLevel > 3) {
    G4cout << "G4LowEPPolarizedComptonModel::SampleSecondaries() E(MeV)= "
           << photonEnergy0/MeV << " in " << couple->GetMaterial()->GetName()
           << G4endl;
  }
  // do nothing below the threshold
  // should never get here because the XS is zero below the limit
  if (photonEnergy0 < LowEnergyLimit())
    return ;

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  
  // Polarisation: check orientation of photon propagation direction and polarisation
  // Fix if needed

  G4ThreeVector photonPolarization0 = aDynamicGamma->GetPolarization();

  // Check if polarisation vector is perpendicular and fix if not

  if (!(photonPolarization0.isOrthogonal(photonDirection0, 1e-6))||(photonPolarization0.mag()==0))
    {
      photonPolarization0 = GetRandomPolarization(photonDirection0);
    }

  else
    {
      if ((photonPolarization0.howOrthogonal(photonDirection0) !=0) && (photonPolarization0.howOrthogonal(photonDirection0) > g4d_limit))
        {
          photonPolarization0 = GetPerpendicularPolarization(photonDirection0,photonPolarization0);
        }
     }

  // Select randomly one element in the current material

  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double LowEPPCepsilon0 = 1. / (1. + 2. * e0m);
  G4double LowEPPCepsilon0Sq = LowEPPCepsilon0 * LowEPPCepsilon0;
  G4double alpha1 = -std::log(LowEPPCepsilon0);
  G4double alpha2 = 0.5 * (1. - LowEPPCepsilon0Sq);

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  // Sample the energy of the scattered photon
  G4double LowEPPCepsilon;
  G4double LowEPPCepsilonSq;
  G4double oneCosT;
  G4double sinT2;
  G4double gReject;

  if (verboseLevel > 3) {
    G4cout << "Started loop to sample gamma energy" << G4endl;
  }

  do
    {
      if ( alpha1/(alpha1+alpha2) > G4UniformRand())
        {
          LowEPPCepsilon = G4Exp(-alpha1 * G4UniformRand());
          LowEPPCepsilonSq = LowEPPCepsilon * LowEPPCepsilon;
        }
      else
        {
          LowEPPCepsilonSq = LowEPPCepsilon0Sq + (1. - LowEPPCepsilon0Sq) * G4UniformRand();
          LowEPPCepsilon = std::sqrt(LowEPPCepsilonSq);
        }

      oneCosT = (1. - LowEPPCepsilon) / ( LowEPPCepsilon * e0m);
      sinT2 = oneCosT * (2. - oneCosT);
      G4double x = std::sqrt(oneCosT/2.) / (wlPhoton/cm);
      G4double scatteringFunction = ComputeScatteringFunction(x, Z);
      gReject = (1. - LowEPPCepsilon * sinT2 / (1. + LowEPPCepsilonSq)) * scatteringFunction;

    } while(gReject < G4UniformRand()*Z);

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = std::sqrt(sinT2);
  G4double phi = SetPhi(LowEPPCepsilon,sinT2);
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  // Set outgoing photon polarization

  G4ThreeVector photonPolarization1 = SetNewPolarization(LowEPPCepsilon,
                                                        sinT2,
                                                        phi,
                                                        cosTheta);

  // Scatter photon energy and Compton electron direction - Method based on:
  // J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin'
  // "A low energy bound atomic electron Compton scattering model for Geant4"
  // NIMB, Vol. 338, 77-88, 2014.

  // Set constants and initialize scattering parameters

  const G4double vel_c = c_light / (m/s);
  const G4double momentum_au_to_nat = halfpi* hbar_Planck / Bohr_radius / (kg*m/s);
  const G4double e_mass_kg =  electron_mass_c2 / c_squared / kg ;

  const G4int maxDopplerIterations = 1000;
  G4double bindingE = 0.;
  G4double pEIncident = photonEnergy0 ;
  G4double pERecoil =  -1.;
  G4double eERecoil = -1.;
  G4double e_alpha =0.;
  G4double e_beta = 0.;

  G4double CE_emission_flag = 0.;
  G4double ePAU = -1;
  G4int shellIdx = 0;
  G4double u_temp = 0;
  G4double cosPhiE =0;
  G4double sinThetaE =0;
  G4double cosThetaE =0;
  G4int iteration = 0;

  if (verboseLevel > 3) {
    G4cout << "Started loop to sample photon energy and electron direction" << G4endl;
  }

  do{


      // ******************************************
      // |     Determine scatter photon energy    |
      // ******************************************   

 do
    {
      iteration++;


      // ********************************************
      // |     Sample bound electron information    |
      // ********************************************

      // Select shell based on shell occupancy

      shellIdx = shellData->SelectRandomShell(Z);
      bindingE = shellData->BindingEnergy(Z,shellIdx)/MeV;


      // Randomly sample bound electron momentum (memento: the data set is in Atomic Units)
      ePAU = profileData->RandomSelectMomentum(Z,shellIdx);

      // Convert to SI units     
      G4double ePSI = ePAU * momentum_au_to_nat;

      //Calculate bound electron velocity and normalise to natural units
      u_temp = sqrt( ((ePSI*ePSI)*(vel_c*vel_c)) / ((e_mass_kg*e_mass_kg)*(vel_c*vel_c)+(ePSI*ePSI)) )/vel_c;

      // Sample incident electron direction, amorphous material, to scattering photon scattering plane 

      e_alpha = pi*G4UniformRand();
      e_beta = twopi*G4UniformRand();

      // Total energy of system  

      G4double eEIncident = electron_mass_c2 / sqrt( 1 - (u_temp*u_temp));
      G4double systemE = eEIncident + pEIncident;


      G4double gamma_temp = 1.0 / sqrt( 1 - (u_temp*u_temp));
      G4double numerator = gamma_temp*electron_mass_c2*(1 - u_temp * std::cos(e_alpha));
      G4double subdenom1 =  u_temp*cosTheta*std::cos(e_alpha);
      G4double subdenom2 = u_temp*sinTheta*std::sin(e_alpha)*std::cos(e_beta);
      G4double denominator = (1.0 - cosTheta) +  (gamma_temp*electron_mass_c2*(1 - subdenom1 - subdenom2) / pEIncident);
      pERecoil = (numerator/denominator);
      eERecoil = systemE - pERecoil;
      CE_emission_flag = pEIncident - pERecoil;
    } while ( (iteration <= maxDopplerIterations) && (CE_emission_flag < bindingE));

// End of recalculation of photon energy with Doppler broadening



   // *******************************************************
   // |     Determine ejected Compton electron direction    |
   // *******************************************************      

      // Calculate velocity of ejected Compton electron   

      G4double a_temp = eERecoil / electron_mass_c2;
      G4double u_p_temp = sqrt(1 - (1 / (a_temp*a_temp)));

      // Coefficients and terms from simulatenous equations     

      G4double sinAlpha = std::sin(e_alpha);
      G4double cosAlpha = std::cos(e_alpha);
      G4double sinBeta = std::sin(e_beta);
      G4double cosBeta = std::cos(e_beta);

      G4double gamma = 1.0 / sqrt(1 - (u_temp*u_temp));
      G4double gamma_p = 1.0 / sqrt(1 - (u_p_temp*u_p_temp));

      G4double var_A = pERecoil*u_p_temp*sinTheta;
      G4double var_B = u_p_temp* (pERecoil*cosTheta-pEIncident);
      G4double var_C = (pERecoil-pEIncident) - ( (pERecoil*pEIncident) / (gamma_p*electron_mass_c2))*(1 - cosTheta);

      G4double var_D1 = gamma*electron_mass_c2*pERecoil;
      G4double var_D2 = (1 - (u_temp*cosTheta*cosAlpha) - (u_temp*sinTheta*cosBeta*sinAlpha));
      G4double var_D3 = ((electron_mass_c2*electron_mass_c2)*(gamma*gamma_p - 1)) - (gamma_p*electron_mass_c2*pERecoil);
      G4double var_D = var_D1*var_D2 + var_D3;

      G4double var_E1 = ((gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*cosAlpha);
      G4double var_E2 = gamma_p*electron_mass_c2*pERecoil*u_p_temp*cosTheta;
      G4double var_E = var_E1 - var_E2;

      G4double var_F1 = ((gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*cosBeta*sinAlpha);
      G4double var_F2 = (gamma_p*electron_mass_c2*pERecoil*u_p_temp*sinTheta);
      G4double var_F = var_F1 - var_F2;

      G4double var_G = (gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*sinBeta*sinAlpha;

      // Two equations form a quadratic form of Wx^2 + Yx + Z = 0
      // Coefficents and solution to quadratic

      G4double var_W1 = (var_F*var_B - var_E*var_A)*(var_F*var_B - var_E*var_A);
      G4double var_W2 = (var_G*var_G)*(var_A*var_A) + (var_G*var_G)*(var_B*var_B);
      G4double var_W = var_W1 + var_W2;

      G4double var_Y = 2.0*(((var_A*var_D-var_F*var_C)*(var_F*var_B-var_E*var_A)) - ((var_G*var_G)*var_B*var_C));

      G4double var_Z1 = (var_A*var_D - var_F*var_C)*(var_A*var_D - var_F*var_C);
      G4double var_Z2 = (var_G*var_G)*(var_C*var_C) - (var_G*var_G)*(var_A*var_A);
      G4double var_Z = var_Z1 + var_Z2;
      G4double diff1 = var_Y*var_Y;
      G4double diff2 = 4*var_W*var_Z;
      G4double diff = diff1 - diff2;


     // Check if diff is less than zero, if so ensure it is due to FPE

     //Confirm that diff less than zero is due FPE, i.e if abs of diff / diff1 and diff/ diff2 is less 
     //than 10^(-g4d_order), then set diff to zero

     if ((diff < 0.0) && (abs(diff / diff1) < g4d_limit) && (abs(diff / diff2) < g4d_limit) )
     {
           diff = 0.0;
     }

      // Plus and minus of quadratic
      G4double X_p = (-var_Y + sqrt (diff))/(2*var_W);
      G4double X_m = (-var_Y - sqrt (diff))/(2*var_W);


      // Floating point precision protection
      // Check if X_p and X_m are greater than or less than 1 or -1, if so clean up FPE 
      // Issue due to propagation of FPE and only impacts 8th sig fig onwards

      if(X_p >1){X_p=1;} if(X_p<-1){X_p=-1;}
      if(X_m >1){X_m=1;} if(X_m<-1){X_m=-1;}

      // End of FP protection

      G4double ThetaE = 0.;


      // Randomly sample one of the two possible solutions and determin theta angle of ejected Compton electron
       G4double sol_select = G4UniformRand();

      if (sol_select < 0.5)
      {
           ThetaE = std::acos(X_p);
      }
      if (sol_select > 0.5)
      {
          ThetaE = std::acos(X_m);
      }

      cosThetaE = std::cos(ThetaE);
      sinThetaE = std::sin(ThetaE);
      G4double Theta = std::acos(cosTheta);

      //Calculate electron Phi
      G4double iSinThetaE = std::sqrt(1+std::tan((pi/2.0)-ThetaE)*std::tan((pi/2.0)-ThetaE));
      G4double iSinTheta = std::sqrt(1+std::tan((pi/2.0)-Theta)*std::tan((pi/2.0)-Theta));
      G4double ivar_A = iSinTheta/ (pERecoil*u_p_temp);
      // Trigs
      cosPhiE = (var_C - var_B*cosThetaE)*(ivar_A*iSinThetaE);

     // End of calculation of ejection Compton electron direction

      //Fix for floating point errors

    } while ( (iteration <= maxDopplerIterations) && (abs(cosPhiE) > 1));

   // Revert to original if maximum number of iterations threshold has been reached     
  if (iteration >= maxDopplerIterations)
    {
      pERecoil = photonEnergy0 ;
      bindingE = 0.;
      dirx=0.0;
      diry=0.0;
      dirz=1.0;
    }

  // Set "scattered" photon direction and energy

  G4ThreeVector photonDirection1(dirx,diry,dirz);
  SystemOfRefChange(photonDirection0,photonDirection1,
		    photonPolarization0,photonPolarization1);


  if (pERecoil > 0.)
    {
     fParticleChange->SetProposedKineticEnergy(pERecoil) ;
     fParticleChange->ProposeMomentumDirection(photonDirection1) ;
     fParticleChange->ProposePolarization(photonPolarization1);

     // Set ejected Compton electron direction and energy
     G4double PhiE = std::acos(cosPhiE);
     G4double eDirX = sinThetaE * std::cos(phi+PhiE);
     G4double eDirY = sinThetaE * std::sin(phi+PhiE);
     G4double eDirZ = cosThetaE;

     G4double eKineticEnergy = pEIncident - pERecoil - bindingE;

     G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
     SystemOfRefChangeElect(photonDirection0,eDirection,
                  photonPolarization0);

     G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
                                                   eDirection,eKineticEnergy) ;
     fvect->push_back(dp);

    }
  else
    {
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeTrackStatus(fStopAndKill);
    }

  // sample deexcitation
  //

  if (verboseLevel > 3) {
    G4cout << "Started atomic de-excitation " << fAtomDeexcitation << G4endl;
  }

  if(fAtomDeexcitation && iteration < maxDopplerIterations) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
        for (size_t i=nbefore; i<nafter; ++i) {
          //Check if there is enough residual energy 
          if (bindingE >= ((*fvect)[i])->GetKineticEnergy())
           {
             //Ok, this is a valid secondary: keep it
             bindingE -= ((*fvect)[i])->GetKineticEnergy();
           }
          else
           {
             //Invalid secondary: not enough energy to create it!
             //Keep its energy in the local deposit
             delete (*fvect)[i];
             (*fvect)[i]=0;
           }
        }
      }
    }
  }

  //This should never happen
  if(bindingE < 0.0)
     G4Exception("G4LowEPPolarizedComptonModel::SampleSecondaries()",
                 "em2051",FatalException,"Negative local energy deposit");

  fParticleChange->ProposeLocalEnergyDeposit(bindingE);

}

//****************************************************************************

G4double
G4LowEPPolarizedComptonModel::ComputeScatteringFunction(G4double x, G4int Z)
{
  G4double value = Z;
  if (x <= ScatFuncFitParam[Z][2]) {

    G4double lgq = G4Log(x)/ln10;

    if (lgq < ScatFuncFitParam[Z][1]) {
      value = ScatFuncFitParam[Z][3] + lgq*ScatFuncFitParam[Z][4];
    } else {
      value = ScatFuncFitParam[Z][5] + lgq*ScatFuncFitParam[Z][6] +
        lgq*lgq*ScatFuncFitParam[Z][7] + lgq*lgq*lgq*ScatFuncFitParam[Z][8];
    }
    value = G4Exp(value*ln10);
  }
  return value;
}


//****************************************************************************

#include "G4AutoLock.hh"
namespace { G4Mutex LowEPPolarizedComptonModelMutex = G4MUTEX_INITIALIZER; }

void
G4LowEPPolarizedComptonModel::InitialiseForElement(const G4ParticleDefinition*,
                                              G4int Z)
{
  G4AutoLock l(&LowEPPolarizedComptonModelMutex);
  if(!data[Z]) { ReadData(Z); }
  l.unlock();
}

//****************************************************************************

//Fitting data to compute scattering function 

const G4double G4LowEPPolarizedComptonModel::ScatFuncFitParam[101][9] = {
{  0,    0.,          0.,      0.,    0.,       0.,     0.,     0.,    0.},
{  1, 6.673, 1.49968E+08, -14.352, 1.999, -143.374, 50.787, -5.951, 0.2304418},
{  2, 6.500, 2.50035E+08, -14.215, 1.970, -53.649, 13.892, -0.948, 0.006996759},
{  3, 6.551, 3.99945E+08, -13.555, 1.993, -62.090, 21.462, -2.453, 0.093416},
{  4, 6.500, 5.00035E+08, -13.746, 1.998, -127.906, 46.491, -5.614, 0.2262103},
{  5, 6.500, 5.99791E+08, -13.800, 1.998, -131.153, 47.132, -5.619, 0.2233819},
{  6, 6.708, 6.99842E+08, -13.885, 1.999, -128.143, 45.379, -5.325, 0.2083009},
{  7, 6.685, 7.99834E+08, -13.885, 2.000, -131.048, 46.314, -5.421, 0.2114925},
{  8, 6.669, 7.99834E+08, -13.962, 2.001, -128.225, 44.818, -5.183, 0.1997155},
{  9, 6.711, 7.99834E+08, -13.999, 2.000, -122.112, 42.103, -4.796, 0.1819099},
{ 10, 6.702, 7.99834E+08, -14.044, 1.999, -110.143, 37.225, -4.143, 0.1532094},
{ 11, 6.425, 1.00000E+09, -13.423, 1.993, -41.137, 12.313, -1.152, 0.03384553},
{ 12, 6.542, 1.00000E+09, -13.389, 1.997, -53.549, 17.420, -1.840, 0.06431849},
{ 13, 6.570, 1.49968E+09, -13.401, 1.997, -66.243, 22.297, -2.460, 0.09045854},
{ 14, 6.364, 1.49968E+09, -13.452, 1.999, -78.271, 26.757, -3.008, 0.1128195},
{ 15, 6.500, 1.49968E+09, -13.488, 1.998, -85.069, 29.164, -3.291, 0.1239113},
{ 16, 6.500, 1.49968E+09, -13.532, 1.998, -93.640, 32.274, -3.665, 0.1388633},
{ 17, 6.500, 1.49968E+09, -13.584, 2.000, -98.534, 33.958, -3.857, 0.1461557},
{ 18, 6.500, 1.49968E+09, -13.618, 1.999, -100.077, 34.379, -3.891, 0.1468902},
{ 19, 6.500, 1.99986E+09, -13.185, 1.992, -53.819, 17.528, -1.851, 0.0648722},
{ 20, 6.490, 1.99986E+09, -13.123, 1.993, -52.221, 17.169, -1.832, 0.06502094},
{ 21, 6.498, 1.99986E+09, -13.157, 1.994, -55.365, 18.276, -1.961, 0.07002778},
{ 22, 6.495, 1.99986E+09, -13.183, 1.994, -57.412, 18.957, -2.036, 0.07278856},
{ 23, 6.487, 1.99986E+09, -13.216, 1.995, -58.478, 19.270, -2.065, 0.07362722},
{ 24, 6.500, 1.99986E+09, -13.330, 1.997, -62.192, 20.358, -2.167, 0.07666583},
{ 25, 6.488, 1.99986E+09, -13.277, 1.997, -58.007, 18.924, -2.003, 0.0704305},
{ 26, 6.500, 5.00035E+09, -13.292, 1.997, -61.176, 20.067, -2.141, 0.0760269},
{ 27, 6.500, 5.00035E+09, -13.321, 1.998, -61.909, 20.271, -2.159, 0.07653559},
{ 28, 6.500, 5.00035E+09, -13.340, 1.998, -62.402, 20.391, -2.167, 0.07664243},
{ 29, 6.500, 5.00035E+09, -13.439, 1.998, -67.305, 21.954, -2.331, 0.0823267},
{ 30, 6.500, 5.00035E+09, -13.383, 1.999, -62.064, 20.136, -2.122, 0.07437589},
{ 31, 6.500, 5.00035E+09, -13.349, 1.997, -61.068, 19.808, -2.086, 0.07307488},
{ 32, 6.500, 5.00035E+09, -13.373, 1.999, -63.126, 20.553, -2.175, 0.07660222},
{ 33, 6.500, 5.00035E+09, -13.395, 1.999, -65.674, 21.445, -2.278, 0.08054694},
{ 34, 6.500, 5.00035E+09, -13.417, 1.999, -69.457, 22.811, -2.442, 0.08709536},
{ 35, 6.500, 5.00035E+09, -13.442, 2.000, -72.283, 23.808, -2.558, 0.09156808},
{ 36, 6.500, 5.00035E+09, -13.451, 1.998, -74.696, 24.641, -2.653, 0.09516597},
{ 37, 6.500, 5.00035E+09, -13.082, 1.991, -46.235, 14.519, -1.458, 0.04837659},
{ 38, 6.465, 5.00035E+09, -13.022, 1.993, -41.784, 13.065, -1.300, 0.04267703},
{ 39, 6.492, 5.00035E+09, -13.043, 1.994, -44.609, 14.114, -1.429, 0.0479348},
{ 40, 6.499, 5.00035E+09, -13.064, 1.994, -47.142, 15.019, -1.536, 0.0521347},
{ 41, 6.384, 5.00035E+09, -13.156, 1.996, -53.114, 17.052, -1.766, 0.06079426},
{ 42, 6.500, 5.00035E+09, -13.176, 1.996, -54.590, 17.550, -1.822, 0.06290335},
{ 43, 6.500, 5.00035E+09, -13.133, 1.997, -51.272, 16.423, -1.694, 0.05806108},
{ 44, 6.500, 5.00035E+09, -13.220, 1.996, -58.314, 18.839, -1.969, 0.0684608},
{ 45, 6.500, 5.00035E+09, -13.246, 1.998, -59.674, 19.295, -2.020, 0.07037294},
{ 46, 6.500, 5.00035E+09, -13.407, 1.999, -72.228, 23.693, -2.532, 0.09017969},
{ 47, 6.500, 5.00035E+09, -13.277, 1.998, -60.890, 19.647, -2.053, 0.07138694},
{ 48, 6.500, 5.00035E+09, -13.222, 1.998, -56.152, 18.002, -1.863, 0.06410123},
{ 49, 6.500, 5.00035E+09, -13.199, 1.997, -56.208, 18.052, -1.872, 0.06456884},
{ 50, 6.500, 5.00035E+09, -13.215, 1.998, -58.478, 18.887, -1.973, 0.06860356},
{ 51, 6.500, 5.00035E+09, -13.230, 1.998, -60.708, 19.676, -2.066, 0.07225841},
{ 52, 6.500, 7.99834E+09, -13.246, 1.998, -63.341, 20.632, -2.180, 0.0767412},
{ 53, 6.500, 5.00035E+09, -13.262, 1.998, -66.339, 21.716, -2.310, 0.08191981},
{ 54, 6.500, 7.99834E+09, -13.279, 1.998, -67.649, 22.151, -2.357, 0.08357825},
{ 55, 6.500, 5.00035E+09, -12.951, 1.990, -45.302, 14.219, -1.423, 0.04712317},
{ 56, 6.425, 5.00035E+09, -12.882, 1.992, -39.825, 12.363, -1.214, 0.03931009},
{ 57, 6.466, 2.82488E+09, -12.903, 1.992, -38.952, 11.982, -1.160, 0.03681554},
{ 58, 6.451, 5.00035E+09, -12.915, 1.993, -41.959, 13.118, -1.302, 0.04271291},
{ 59, 6.434, 5.00035E+09, -12.914, 1.993, -40.528, 12.555, -1.230, 0.03971407},
{ 60, 6.444, 5.00035E+09, -12.922, 1.992, -39.986, 12.329, -1.200, 0.03843737},
{ 61, 6.414, 7.99834E+09, -12.930, 1.993, -42.756, 13.362, -1.327, 0.0436124},
{ 62, 6.420, 7.99834E+09, -12.938, 1.992, -42.682, 13.314, -1.319, 0.04322509},
{ 63, 6.416, 7.99834E+09, -12.946, 1.993, -42.399, 13.185, -1.301, 0.04243861},
{ 64, 6.443, 7.99834E+09, -12.963, 1.993, -43.226, 13.475, -1.335, 0.04377341},
{ 65, 6.449, 7.99834E+09, -12.973, 1.993, -43.232, 13.456, -1.330, 0.04347536},
{ 66, 6.419, 7.99834E+09, -12.966, 1.993, -42.047, 12.990, -1.270, 0.04095499},
{ 67, 6.406, 7.99834E+09, -12.976, 1.993, -42.405, 13.106, -1.283, 0.04146024},
{ 68, 6.424, 7.99834E+09, -12.986, 1.993, -41.974, 12.926, -1.259, 0.040435},
{ 69, 6.417, 7.99834E+09, -12.989, 1.993, -42.132, 12.967, -1.262, 0.04048908},
{ 70, 6.405, 7.99834E+09, -13.000, 1.994, -42.582, 13.122, -1.280, 0.04119599},
{ 71, 6.449, 7.99834E+09, -13.015, 1.994, -42.586, 13.115, -1.278, 0.04107587},
{ 72, 6.465, 7.99834E+09, -13.030, 1.994, -43.708, 13.509, -1.324, 0.04286491},
{ 73, 6.447, 7.99834E+09, -13.048, 1.996, -44.838, 13.902, -1.369, 0.04457132},
{ 74, 6.452, 7.99834E+09, -13.073, 1.997, -45.545, 14.137, -1.395, 0.04553459},
{ 75, 6.432, 7.99834E+09, -13.082, 1.997, -46.426, 14.431, -1.428, 0.04678218},
{ 76, 6.439, 7.99834E+09, -13.100, 1.997, -47.513, 14.806, -1.471, 0.04842566},
{ 77, 6.432, 7.99834E+09, -13.110, 1.997, -48.225, 15.042, -1.497, 0.04938364},
{ 78, 6.500, 7.99834E+09, -13.185, 1.997, -53.256, 16.739, -1.687, 0.05645173},
{ 79, 6.500, 7.99834E+09, -13.200, 1.997, -53.900, 16.946, -1.709, 0.05723134},
{ 80, 6.500, 7.99834E+09, -13.156, 1.998, -49.801, 15.536, -1.547, 0.05103522},
{ 81, 6.500, 7.99834E+09, -13.128, 1.997, -49.651, 15.512, -1.548, 0.05123203},
{ 82, 6.500, 7.99834E+09, -13.134, 1.997, -51.021, 16.018, -1.609, 0.05364831},
{ 83, 6.500, 7.99834E+09, -13.148, 1.998, -52.693, 16.612, -1.679, 0.05638698},
{ 84, 6.500, 7.99834E+09, -13.161, 1.998, -54.415, 17.238, -1.754, 0.05935566},
{ 85, 6.500, 7.99834E+09, -13.175, 1.998, -56.083, 17.834, -1.824, 0.06206068},
{ 86, 6.500, 7.99834E+09, -13.189, 1.998, -57.860, 18.463, -1.898, 0.0649633},
{ 87, 6.500, 7.99834E+09, -12.885, 1.990, -39.973, 12.164, -1.162, 0.0364598},
{ 88, 6.417, 7.99834E+09, -12.816, 1.991, -34.591, 10.338, -0.956, 0.0287409},
{ 89, 6.442, 7.99834E+09, -12.831, 1.992, -36.002, 10.867, -1.021, 0.03136835},
{ 90, 6.463, 7.99834E+09, -12.850, 1.993, -37.660, 11.475, -1.095, 0.03435334},
{ 91, 6.447, 7.99834E+09, -12.852, 1.993, -37.268, 11.301, -1.071, 0.0330539},
{ 92, 6.439, 7.99834E+09, -12.858, 1.993, -37.695, 11.438, -1.085, 0.03376669},
{ 93, 6.437, 1.00000E+10, -12.866, 1.993, -39.010, 11.927, -1.146, 0.03630848},
{ 94, 6.432, 7.99834E+09, -12.862, 1.993, -37.192, 11.229, -1.057, 0.0325621},
{ 95, 6.435, 7.99834E+09, -12.869, 1.993, -37.589, 11.363, -1.072, 0.03312393},
{ 96, 6.449, 1.00000E+10, -12.886, 1.993, -39.573, 12.095, -1.162, 0.03680527},
{ 97, 6.446, 1.00000E+10, -12.892, 1.993, -40.007, 12.242, -1.178, 0.03737377},
{ 98, 6.421, 1.00000E+10, -12.887, 1.993, -39.509, 12.041, -1.152, 0.03629023},
{ 99, 6.414, 1.00000E+10, -12.894, 1.993, -39.939, 12.183, -1.168, 0.03690464},
{100, 6.412, 1.00000E+10, -12.900, 1.993, -39.973, 12.180, -1.166, 0.036773}
  };

//****************************************************************************

//Supporting functions for photon polarisation effects

G4double G4LowEPPolarizedComptonModel::SetPhi(G4double energyRate,
                                             G4double sinT2)
{
  G4double rand1;
  G4double rand2;
  G4double phiProbability;
  G4double phi;
  G4double a, b;

  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      phiProbability=0.;
      phi = twopi*rand1;

      a = 2*sinT2;
      b = energyRate + 1/energyRate;

      phiProbability = 1 - (a/b)*(std::cos(phi)*std::cos(phi));



    }
  while ( rand2 > phiProbability );
  return phi;
}

//****************************************************************************

G4ThreeVector G4LowEPPolarizedComptonModel::SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}

//****************************************************************************

G4ThreeVector G4LowEPPolarizedComptonModel::GetRandomPolarization(G4ThreeVector& direction0)
{
  G4ThreeVector d0 = direction0.unit();
  G4ThreeVector a1 = SetPerpendicularVector(d0); //different orthogonal
  G4ThreeVector a0 = a1.unit(); // unit vector

  G4double rand1 = G4UniformRand();

  G4double angle = twopi*rand1; // random polar angle
  G4ThreeVector b0 = d0.cross(a0); // cross product

  G4ThreeVector c;

  c.setX(std::cos(angle)*(a0.x())+std::sin(angle)*b0.x());
  c.setY(std::cos(angle)*(a0.y())+std::sin(angle)*b0.y());
  c.setZ(std::cos(angle)*(a0.z())+std::sin(angle)*b0.z());

  G4ThreeVector c0 = c.unit();

  return c0;

}

//****************************************************************************

G4ThreeVector G4LowEPPolarizedComptonModel::GetPerpendicularPolarization
(const G4ThreeVector& photonDirection, const G4ThreeVector& photonPolarization) const
{

  // 
  // The polarization of a photon is always perpendicular to its momentum direction.
  // Therefore this function removes those vector component of photonPolarization, which
  // points in direction of photonDirection
  //
  // Mathematically we search the projection of the vector a on the plane E, where n is the
  // plains normal vector.
  // The basic equation can be found in each geometry book (e.g. Bronstein):
  // p = a - (a o n)/(n o n)*n

  return photonPolarization - photonPolarization.dot(photonDirection)/photonDirection.dot(photonDirection) * photonDirection;
}

//****************************************************************************

G4ThreeVector G4LowEPPolarizedComptonModel::SetNewPolarization(G4double LowEPPCepsilon,
                                                              G4double sinT2,
                                                              G4double phi,
                                                              G4double costheta)
{
  G4double rand1;
  G4double rand2;
  G4double cosPhi = std::cos(phi);
  G4double sinPhi = std::sin(phi);
  G4double sinTheta = std::sqrt(sinT2);
  G4double cosP2 = cosPhi*cosPhi;
  G4double normalisation = std::sqrt(1. - cosP2*sinT2);


  // Method based on:
  // D. Xu, Z. He and F. Zhang
  // "Detection of Gamma Ray Polarization Using a 3-D Position Sensitive CdZnTe Detector"
  // IEEE TNS, Vol. 52(4), 1160-1164, 2005.

  // Determination of Theta 

  G4double theta;

  rand1 = G4UniformRand();
  rand2 = G4UniformRand();

  if (rand1<(LowEPPCepsilon+1.0/LowEPPCepsilon-2)/(2.0*(LowEPPCepsilon+1.0/LowEPPCepsilon)-4.0*sinT2*cosP2))
    {
      if (rand2<0.5)
        theta = pi/2.0;
      else
        theta = 3.0*pi/2.0;
    }
  else
    {
      if (rand2<0.5)
        theta = 0;
      else
        theta = pi;
    }
  G4double cosBeta = std::cos(theta);
  G4double sinBeta = std::sqrt(1-cosBeta*cosBeta);

  G4ThreeVector photonPolarization1;

  G4double xParallel = normalisation*cosBeta;
  G4double yParallel = -(sinT2*cosPhi*sinPhi)*cosBeta/normalisation;
  G4double zParallel = -(costheta*sinTheta*cosPhi)*cosBeta/normalisation;
  G4double xPerpendicular = 0.;
  G4double yPerpendicular = (costheta)*sinBeta/normalisation;
  G4double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;

  G4double xTotal = (xParallel + xPerpendicular);
  G4double yTotal = (yParallel + yPerpendicular);
  G4double zTotal = (zParallel + zPerpendicular);

  photonPolarization1.setX(xTotal);
  photonPolarization1.setY(yTotal);
  photonPolarization1.setZ(zTotal);

  return photonPolarization1;

}
void G4LowEPPolarizedComptonModel::SystemOfRefChange(G4ThreeVector& direction0,
                                                    G4ThreeVector& direction1,
                                                    G4ThreeVector& polarization0,
                                                    G4ThreeVector& polarization1)
{
  // direction0 is the original photon direction ---> z
  // polarization0 is the original photon polarization ---> x
  // need to specify y axis in the real reference frame ---> y 
  G4ThreeVector Axis_Z0 = direction0.unit();
  G4ThreeVector Axis_X0 = polarization0.unit();
  G4ThreeVector Axis_Y0 = (Axis_Z0.cross(Axis_X0)).unit(); // to be confirmed;

  G4double direction_x = direction1.getX();
  G4double direction_y = direction1.getY();
  G4double direction_z = direction1.getZ();
  
  direction1 = (direction_x*Axis_X0 + direction_y*Axis_Y0 + direction_z*Axis_Z0).unit();
  G4double polarization_x = polarization1.getX();
  G4double polarization_y = polarization1.getY();
  G4double polarization_z = polarization1.getZ();

  polarization1 = (polarization_x*Axis_X0 + polarization_y*Axis_Y0 + polarization_z*Axis_Z0).unit();

}

void G4LowEPPolarizedComptonModel::SystemOfRefChangeElect(G4ThreeVector& pdirection,
                                                    G4ThreeVector& edirection,
                                                    G4ThreeVector& ppolarization)
{
  // direction0 is the original photon direction ---> z
  // polarization0 is the original photon polarization ---> x
  // need to specify y axis in the real reference frame ---> y 
  G4ThreeVector Axis_Z0 = pdirection.unit();
  G4ThreeVector Axis_X0 = ppolarization.unit();
  G4ThreeVector Axis_Y0 = (Axis_Z0.cross(Axis_X0)).unit(); // to be confirmed;

  G4double direction_x = edirection.getX();
  G4double direction_y = edirection.getY();
  G4double direction_z = edirection.getZ();

  edirection = (direction_x*Axis_X0 + direction_y*Axis_Y0 + direction_z*Axis_Z0).unit();

}


