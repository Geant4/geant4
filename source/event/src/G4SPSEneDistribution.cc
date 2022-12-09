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
// G4SPSEneDistribution class implementation
//
// Author: Fan Lei, QinetiQ ltd - 05/02/2004
// Customer: ESA/ESTEC
// Revisions: Andrew Green, Andrea Dotti
// --------------------------------------------------------------------
#include "G4SPSEneDistribution.hh"

#include "G4Exp.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"

G4SPSEneDistribution::G4SPSEneDistribution()
{
  G4MUTEXINIT(mutex);

  // Initialise all variables

  particle_energy = 1.0 * MeV;
  EnergyDisType = "Mono";
  weight=1.;
  MonoEnergy = 1 * MeV;
  Emin = 0.;
  Emax = 1.e30;
  alpha = 0.;
  biasalpha = 0.;
  prob_norm = 1.0;
  Ezero = 0.;
  SE = 0.;
  Temp = 0.;
  grad = 0.;
  cept = 0.;
  IntType = "NULL"; // Interpolation type

  ArbEmin = 0.;
  ArbEmax = 1.e30;

  verbosityLevel = 0;
    
  threadLocal_t& data = threadLocalData.Get();
  data.Emax = Emax;
  data.Emin = Emin;
  data.alpha =alpha;
  data.cept = cept;
  data.Ezero = Ezero;
  data.grad = grad;
  data.particle_energy = 0.;
  data.particle_definition = nullptr;
  data.weight = weight;
}

G4SPSEneDistribution::~G4SPSEneDistribution()
{
  G4MUTEXDESTROY(mutex);
  if(Arb_grad_cept_flag)
  {
    delete [] Arb_grad;
    delete [] Arb_cept;
  }
    
  if(Arb_alpha_Const_flag)
  {
    delete [] Arb_alpha;
    delete [] Arb_Const;
  }
    
  if(Arb_ezero_flag)
  {
    delete [] Arb_ezero;
  }
  delete Bbody_x;
  delete BBHist;
  delete CP_x;
  delete CPHist;
  for (auto & it : SplineInt)
  {
    delete it;
    it = nullptr;
  }
  SplineInt.clear();
}

void G4SPSEneDistribution::SetEnergyDisType(const G4String& DisType)
{
  G4AutoLock l(&mutex);
  EnergyDisType = DisType;
  if (EnergyDisType == "User")
  {
    UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
    IPDFEnergyExist = false;
  }
  else if (EnergyDisType == "Arb")
  {
    ArbEnergyH = IPDFArbEnergyH = ZeroPhysVector;
    IPDFArbExist = false;
  }
  else if (EnergyDisType == "Epn")
  {
    UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
    IPDFEnergyExist = false;
    EpnEnergyH = ZeroPhysVector;
  }
}

const G4String& G4SPSEneDistribution::GetEnergyDisType()
{
  G4AutoLock l(&mutex);
  return EnergyDisType;
}

void G4SPSEneDistribution::SetEmin(G4double emi)
{
  G4AutoLock l(&mutex);
  Emin = emi;
  threadLocalData.Get().Emin = Emin;
}

G4double G4SPSEneDistribution::GetEmin() const
{
  return threadLocalData.Get().Emin;
}

G4double G4SPSEneDistribution::GetArbEmin()
{
  G4AutoLock l(&mutex);
  return ArbEmin;
}

G4double G4SPSEneDistribution::GetArbEmax()
{
  G4AutoLock l(&mutex);
  return ArbEmax;
}

void G4SPSEneDistribution::SetEmax(G4double ema)
{
  G4AutoLock l(&mutex);
  Emax = ema;
  threadLocalData.Get().Emax = Emax;
}

G4double G4SPSEneDistribution::GetEmax() const
{
  return threadLocalData.Get().Emax;
}

void G4SPSEneDistribution::SetMonoEnergy(G4double menergy)
{
  G4AutoLock l(&mutex);
  MonoEnergy = menergy;
}

void G4SPSEneDistribution::SetBeamSigmaInE(G4double e)
{
  G4AutoLock l(&mutex);
  SE = e;
}
void G4SPSEneDistribution::SetAlpha(G4double alp)
{
  G4AutoLock l(&mutex);
  alpha = alp;
  threadLocalData.Get().alpha = alpha;
}

void G4SPSEneDistribution::SetBiasAlpha(G4double alp)
{
  G4AutoLock l(&mutex);
  biasalpha = alp;
  Biased = true;
}

void G4SPSEneDistribution::SetTemp(G4double tem)
{
  G4AutoLock l(&mutex);
  Temp = tem;
}

void G4SPSEneDistribution::SetEzero(G4double eze)
{
  G4AutoLock l(&mutex);
  Ezero = eze;
  threadLocalData.Get().Ezero = Ezero;
}

void G4SPSEneDistribution::SetGradient(G4double gr)
{
  G4AutoLock l(&mutex);
  grad = gr;
  threadLocalData.Get().grad = grad;
}

void G4SPSEneDistribution::SetInterCept(G4double c)
{
  G4AutoLock l(&mutex);
  cept = c;
  threadLocalData.Get().cept = cept;
}

const G4String& G4SPSEneDistribution::GetIntType()
{
  G4AutoLock l(&mutex);
  return IntType;
}

void G4SPSEneDistribution::SetBiasRndm(G4SPSRandomGenerator* a)
{
  G4AutoLock l(&mutex);
  eneRndm = a;
}

void G4SPSEneDistribution::SetVerbosity(G4int a)
{
  G4AutoLock l(&mutex);
  verbosityLevel = a;
}

G4double G4SPSEneDistribution::GetWeight() const
{
  return threadLocalData.Get().weight;
}

G4double G4SPSEneDistribution::GetMonoEnergy()
{
  G4AutoLock l(&mutex);
  return MonoEnergy;
}

G4double G4SPSEneDistribution::GetSE()
{
  G4AutoLock l(&mutex);
  return SE;
}

G4double G4SPSEneDistribution::Getalpha() const
{
  return threadLocalData.Get().alpha;
}

G4double G4SPSEneDistribution::GetEzero() const
{
  return threadLocalData.Get().Ezero;
}

G4double G4SPSEneDistribution::GetTemp()
{
  G4AutoLock l(&mutex);
  return Temp;
}

G4double G4SPSEneDistribution::Getgrad() const
{
  return threadLocalData.Get().grad;
}

G4double G4SPSEneDistribution::Getcept() const
{
  return threadLocalData.Get().cept;
}

G4PhysicsFreeVector G4SPSEneDistribution::GetUserDefinedEnergyHisto()
{
  G4AutoLock l(&mutex);
  return UDefEnergyH;
}

G4PhysicsFreeVector G4SPSEneDistribution::GetArbEnergyHisto()
{
  G4AutoLock l(&mutex);
  return ArbEnergyH;
}

void G4SPSEneDistribution::UserEnergyHisto(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  G4double ehi = input.x(),
           val = input.y();
  if (verbosityLevel > 1)
  {
    G4cout << "In UserEnergyHisto" << G4endl;
    G4cout << " " << ehi << " " << val << G4endl;
  }
  UDefEnergyH.InsertValues(ehi, val);
  Emax = ehi;
  threadLocalData.Get().Emax = Emax;
}

void G4SPSEneDistribution::ArbEnergyHisto(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  G4double ehi = input.x(),
           val = input.y();
  if (verbosityLevel > 1)
  {
    G4cout << "In ArbEnergyHisto" << G4endl;
    G4cout << " " << ehi << " " << val << G4endl;
  }
  ArbEnergyH.InsertValues(ehi, val);
}

void G4SPSEneDistribution::ArbEnergyHistoFile(const G4String& filename)
{
  G4AutoLock l(&mutex);
  std::ifstream infile(filename, std::ios::in);
  if (!infile)
  {
    G4Exception("G4SPSEneDistribution::ArbEnergyHistoFile", "Event0301",
                FatalException, "Unable to open the histo ASCII file");
  }
  G4double ehi, val;
  while (infile >> ehi >> val)
  {
    ArbEnergyH.InsertValues(ehi, val);
  }
}

void G4SPSEneDistribution::EpnEnergyHisto(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  G4double ehi = input.x(),
           val = input.y();
  if (verbosityLevel > 1)
  {
    G4cout << "In EpnEnergyHisto" << G4endl;
    G4cout << " " << ehi << " " << val << G4endl;
  }
  EpnEnergyH.InsertValues(ehi, val);
  Emax = ehi;
  threadLocalData.Get().Emax = Emax;
  Epnflag = true;
}

void G4SPSEneDistribution::Calculate()
{
  G4AutoLock l(&mutex);
  if (EnergyDisType == "Cdg")
  {
    CalculateCdgSpectrum();
  }
  else if (EnergyDisType == "Bbody")
  {
    if(!BBhistInit)
    {
      BBInitHists();
    }
    CalculateBbodySpectrum();
  }
  else if (EnergyDisType == "CPow")
  {
    if(!CPhistInit)
    {
      CPInitHists();
    }
    CalculateCPowSpectrum();
  }
}

void G4SPSEneDistribution::BBInitHists()  // MT: Lock in caller
{
  BBHist = new std::vector<G4double>(10001, 0.0);
  Bbody_x = new std::vector<G4double>(10001, 0.0);
  BBhistInit = true;
}

void G4SPSEneDistribution::CPInitHists()  // MT: Lock in caller
{
  CPHist = new std::vector<G4double>(10001, 0.0);
  CP_x = new std::vector<G4double>(10001, 0.0);
  CPhistInit = true;
}

void G4SPSEneDistribution::CalculateCdgSpectrum()  // MT: Lock in caller
{
  // This uses the spectrum from the INTEGRAL Mass Model (TIMM)
  // to generate a Cosmic Diffuse X/gamma ray spectrum.

  G4double pfact[2] = { 8.5, 112 };
  G4double spind[2] = { 1.4, 2.3 };
  G4double ene_line[3] = { 1. * keV, 18. * keV, 1E6 * keV };
  G4int n_par;

  ene_line[0] = threadLocalData.Get().Emin;
  if (threadLocalData.Get().Emin < 18 * keV)
  {
    n_par = 2;
    ene_line[2] = threadLocalData.Get().Emax;
    if (threadLocalData.Get().Emax < 18 * keV)
    {
      n_par = 1;
      ene_line[1] = threadLocalData.Get().Emax;
    }
  }
  else
  {
    n_par = 1;
    pfact[0] = 112.;
    spind[0] = 2.3;
    ene_line[1] = threadLocalData.Get().Emax;
  }

  // Create a cumulative histogram
  //
  CDGhist[0] = 0.;
  G4double omalpha;
  G4int i = 0;
  while (i < n_par)
  {
    omalpha = 1. - spind[i];
    CDGhist[i + 1] = CDGhist[i] + (pfact[i] / omalpha)
                                * (std::pow(ene_line[i + 1] / keV, omalpha)
                                - std::pow(ene_line[i] / keV,omalpha));
    ++i;
  }

  // Normalise histo and divide by 1000 to make MeV
  //
  i = 0;
  while (i < n_par)
  {
    CDGhist[i + 1] = CDGhist[i + 1] / CDGhist[n_par];
    ++i;
  }
}

void G4SPSEneDistribution::CalculateBbodySpectrum()  // MT: Lock in caller
{
  // Create bbody spectrum
  // Proved very hard to integrate indefinitely, so different method.
  // User inputs emin, emax and T. These are used to create a 10,000
  // bin histogram.
  // Use photon density spectrum = 2 nu**2/c**2 * (std::exp(h nu/kT)-1)
  // = 2 E**2/h**2c**2 times the exponential
    
  G4double erange = threadLocalData.Get().Emax - threadLocalData.Get().Emin;
  G4double steps = erange / 10000.;
    
  const G4double k = 8.6181e-11; //Boltzmann const in MeV/K
  const G4double h = 4.1362e-21; // Plancks const in MeV s
  const G4double c = 3e8; // Speed of light
  const G4double h2 = h * h;
  const G4double c2 = c * c;
  G4int count = 0;
  G4double sum = 0.;
  BBHist->at(0) = 0.;
    
  while (count < 10000)
  {
    Bbody_x->at(count) = threadLocalData.Get().Emin + G4double(count * steps);
    G4double Bbody_y = (2. * std::pow(Bbody_x->at(count), 2.))
                     / (h2*c2*(std::exp(Bbody_x->at(count) / (k*Temp)) - 1.));
    sum = sum + Bbody_y;
    BBHist->at(count + 1) = BBHist->at(count) + Bbody_y;
    ++count;
  }

  Bbody_x->at(10000) = threadLocalData.Get().Emax;

  // Normalise cumulative histo
  //
  count = 0;
  while (count < 10001)
  {
    BBHist->at(count) = BBHist->at(count) / sum;
    ++count;
  }
}

void G4SPSEneDistribution::CalculateCPowSpectrum()  // MT: Lock in caller
{
  // Create cutoff power-law spectrum, x^a exp(-x/b)
  // The integral of this function is an incomplete gamma function, which
  // is only available in the Boost library.
  //
  // User inputs are emin, emax and alpha and Ezero. These are used to
  // create a 10,000 bin histogram.
    
  G4double erange = threadLocalData.Get().Emax - threadLocalData.Get().Emin;
  G4double steps = erange / 10000.;
  alpha = threadLocalData.Get().alpha ;
  Ezero = threadLocalData.Get().Ezero ;
    
  G4int count = 0;
  G4double sum = 0.;
  CPHist->at(0) = 0.;
    
  while (count < 10000)
  {
    CP_x->at(count) = threadLocalData.Get().Emin + G4double(count * steps);
    G4double CP_y = std::pow(CP_x->at(count), alpha)
                  * std::exp(-CP_x->at(count) / Ezero);
    sum = sum + CP_y;
    CPHist->at(count + 1) = CPHist->at(count) + CP_y;
    ++count;
  }

  CP_x->at(10000) = threadLocalData.Get().Emax;

  // Normalise cumulative histo
  //
  count = 0;
  while (count < 10001)
  {
    CPHist->at(count) = CPHist->at(count) / sum;
    ++count;
  }
}

void G4SPSEneDistribution::InputEnergySpectra(G4bool value)
{
  G4AutoLock l(&mutex);

  // Allows user to specify spectrum is momentum
  //
  EnergySpec = value; // false if momentum
  if (verbosityLevel > 1)
  {
    G4cout << "EnergySpec has value " << EnergySpec << G4endl;
  }
}

void G4SPSEneDistribution::InputDifferentialSpectra(G4bool value)
{
  G4AutoLock l(&mutex);

  // Allows user to specify integral or differential spectra
  //
  DiffSpec = value; // true = differential, false = integral
  if (verbosityLevel > 1)
  {
    G4cout << "Diffspec has value " << DiffSpec << G4endl;
  }
}

void G4SPSEneDistribution::ArbInterpolate(const G4String& IType)
{
  G4AutoLock l(&mutex);

  IntType = IType;
  ArbEmax = ArbEnergyH.GetMaxEnergy();
  ArbEmin = ArbEnergyH.Energy(0);

  // Now interpolate points

  if (IntType == "Lin") LinearInterpolation();
  if (IntType == "Log") LogInterpolation();
  if (IntType == "Exp") ExpInterpolation();
  if (IntType == "Spline") SplineInterpolation();
}

void G4SPSEneDistribution::LinearInterpolation()  // MT: Lock in caller
{
  // Method to do linear interpolation on the Arb points
  // Calculate equation of each line segment, max 1024.
  // Calculate Area under each segment
  // Create a cumulative array which is then normalised Arb_Cum_Area

  G4double Area_seg[1024]; // Stores area under each segment
  G4double sum = 0., Arb_x[1024]={0.}, Arb_y[1024]={0.}, Arb_Cum_Area[1024]={0.};
  std::size_t i, count;
  std::size_t maxi = ArbEnergyH.GetVectorLength();
  for (i = 0; i < maxi; ++i)
  {
    Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(i);
    Arb_y[i] = ArbEnergyH(i);
  }

  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy

  if (!DiffSpec)
  {
    // Converts integral point-wise spectra to Differential
    //
    for (count = 0; count < maxi-1; ++count)
    {
      Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
                   / (Arb_x[count + 1] - Arb_x[count]);
    }
    --maxi;
  }

  if (!EnergySpec)
  {
    // change currently stored values (emin etc) which are actually momenta
    // to energies
    //
    G4ParticleDefinition* pdef = threadLocalData.Get().particle_definition;
    if (pdef == nullptr)
    {
      G4Exception("G4SPSEneDistribution::LinearInterpolation",
                  "Event0302", FatalException,
                  "Error: particle not defined");
    }
    else
    {
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent

      G4double mass = pdef->GetPDGMass();

      // Convert point to energy unit and its value to per energy unit
      //
      G4double total_energy;
      for (count = 0; count < maxi; ++count)
      {
        total_energy = std::sqrt((Arb_x[count] * Arb_x[count])
                     + (mass * mass)); // total energy
        Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
        Arb_x[count] = total_energy - mass; // kinetic energy
      }
    }
  }

  i = 1;

  Arb_grad = new G4double [1024];
  Arb_cept = new G4double [1024];
  Arb_grad_cept_flag = true;
    
  Arb_grad[0] = 0.;
  Arb_cept[0] = 0.;
  Area_seg[0] = 0.;
  Arb_Cum_Area[0] = 0.;
  while (i < maxi)
  {
    // calculate gradient and intercept for each segment
    //
    Arb_grad[i] = (Arb_y[i] - Arb_y[i - 1]) / (Arb_x[i] - Arb_x[i - 1]);
    if (verbosityLevel == 2)
    {
      G4cout << Arb_grad[i] << G4endl;
    }
    if (Arb_grad[i] > 0.)
    {
      if (verbosityLevel == 2)
      {
         G4cout << "Arb_grad is positive" << G4endl;
      }
      Arb_cept[i] = Arb_y[i] - (Arb_grad[i] * Arb_x[i]);
    }
    else if (Arb_grad[i] < 0.)
    {
      if (verbosityLevel == 2)
      {
        G4cout << "Arb_grad is negative" << G4endl;
      }
      Arb_cept[i] = Arb_y[i] + (-Arb_grad[i] * Arb_x[i]);
    }
    else
    {
      if (verbosityLevel == 2)
      {
        G4cout << "Arb_grad is 0." << G4endl;
      }
      Arb_cept[i] = Arb_y[i];
    }

    Area_seg[i] = ((Arb_grad[i] / 2)
                * (Arb_x[i] * Arb_x[i] - Arb_x[i - 1] * Arb_x[i - 1])
                + Arb_cept[i] * (Arb_x[i] - Arb_x[i - 1]));
    Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
    sum = sum + Area_seg[i];
    if (verbosityLevel == 2)
    {
      G4cout << Arb_x[i] << Arb_y[i] << Area_seg[i] << sum
             << Arb_grad[i] << G4endl;
    }
    ++i;
  }

  i = 0;
  while (i < maxi)
  {
    Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
    IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
    ++i;
  }

  // now scale the ArbEnergyH, needed by Probability()
  //
  ArbEnergyH.ScaleVector(1., 1./sum);

  if (verbosityLevel >= 1)
  {
    G4cout << "Leaving LinearInterpolation" << G4endl;
    ArbEnergyH.DumpValues();
    IPDFArbEnergyH.DumpValues();
  }
}

void G4SPSEneDistribution::LogInterpolation()  // MT: Lock in caller
{
  // Interpolation based on Logarithmic equations
  // Generate equations of line segments
  // y = Ax**alpha => log y = alpha*logx + logA
  // Find area under line segments
  // Create normalised, cumulative array Arb_Cum_Area

  G4double Area_seg[1024]; // Stores area under each segment
  G4double sum = 0., Arb_x[1024]={0.}, Arb_y[1024]={0.}, Arb_Cum_Area[1024]={0.};
  std::size_t i, count;
  std::size_t maxi = ArbEnergyH.GetVectorLength();
  for (i = 0; i < maxi; ++i)
  {
    Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(i);
    Arb_y[i] = ArbEnergyH(i);
  }

  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy

  if (!DiffSpec)
  {
    // Converts integral point-wise spectra to Differential
    //
    for (count = 0; count < maxi-1; ++count)
    {
      Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
                   / (Arb_x[count + 1] - Arb_x[count]);
    }
    --maxi;
  }

  if (!EnergySpec)
  {
    // Change currently stored values (emin etc) which are actually momenta
    // to energies

    G4ParticleDefinition* pdef = threadLocalData.Get().particle_definition;
    if (pdef == nullptr)
    {
      G4Exception("G4SPSEneDistribution::LogInterpolation",
                  "Event0302", FatalException,
                  "Error: particle not defined");
    }
    else
    {
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent

      G4double mass = pdef->GetPDGMass();

      // Convert point to energy unit and its value to per energy unit
      //
      G4double total_energy;
      for (count = 0; count < maxi; ++count)
      {
        total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass * mass));
        Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
        Arb_x[count] = total_energy - mass; // kinetic energy
      }
    }
  }

  i = 1;

  if ( Arb_ezero != nullptr ) { delete [] Arb_ezero; Arb_ezero = nullptr; }
  if ( Arb_Const != nullptr ) { delete [] Arb_Const; Arb_Const = nullptr; }
  Arb_alpha = new G4double [1024];
  Arb_Const = new G4double [1024];
  Arb_alpha_Const_flag = true;

  Arb_alpha[0] = 0.;
  Arb_Const[0] = 0.;
  Area_seg[0] = 0.;
  Arb_Cum_Area[0] = 0.;
  if (Arb_x[0] <= 0. || Arb_y[0] <= 0.)
  {
    G4cout << "You should not use log interpolation with points <= 0."
           << G4endl;
    G4cout << "These will be changed to 1e-20, which may cause problems"
           << G4endl;
    if (Arb_x[0] <= 0.)
    {
      Arb_x[0] = 1e-20;
    }
    if (Arb_y[0] <= 0.)
    {
      Arb_y[0] = 1e-20;
    }
  }

  G4double alp;
  while (i < maxi)
  {
    // In case points are negative or zero
    //
    if (Arb_x[i] <= 0. || Arb_y[i] <= 0.)
    {
      G4cout << "You should not use log interpolation with points <= 0."
             << G4endl;
      G4cout << "These will be changed to 1e-20, which may cause problems"
             << G4endl;
      if (Arb_x[i] <= 0.)
      {
        Arb_x[i] = 1e-20;
      }
      if (Arb_y[i] <= 0.)
      {
        Arb_y[i] = 1e-20;
      }
    }

    Arb_alpha[i] = (std::log10(Arb_y[i]) - std::log10(Arb_y[i - 1]))
                 / (std::log10(Arb_x[i]) - std::log10(Arb_x[i - 1]));
    Arb_Const[i] = Arb_y[i] / (std::pow(Arb_x[i], Arb_alpha[i]));
    alp = Arb_alpha[i] + 1;
    if (alp == 0.)
    {
      Area_seg[i] = Arb_Const[i]
                  * (std::log(Arb_x[i]) - std::log(Arb_x[i - 1])); 
    }
    else
    {
      Area_seg[i] = (Arb_Const[i] / alp)
                  * (std::pow(Arb_x[i], alp) - std::pow(Arb_x[i - 1], alp));
    }
    sum = sum + Area_seg[i];
    Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
    if (verbosityLevel == 2)
    {
      G4cout << Arb_alpha[i] << Arb_Const[i] << Area_seg[i] << G4endl;
    }
    ++i;
  }

  i = 0;
  while (i < maxi)
  {
    Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
    IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
    ++i;
  }

  // Now scale the ArbEnergyH, needed by Probability()
  //
  ArbEnergyH.ScaleVector(1., 1./sum);

  if (verbosityLevel >= 1)
  {
    G4cout << "Leaving LogInterpolation " << G4endl;
  }
}

void G4SPSEneDistribution::ExpInterpolation()  // MT: Lock in caller
{
  // Interpolation based on Exponential equations
  // Generate equations of line segments
  // y = Ae**-(x/e0) => ln y = -x/e0 + lnA
  // Find area under line segments
  // Create normalised, cumulative array Arb_Cum_Area

  G4double Area_seg[1024]; // Stores area under each segment
  G4double sum = 0., Arb_x[1024]={0.}, Arb_y[1024]={0.}, Arb_Cum_Area[1024]={0.};
  std::size_t i, count;
  std::size_t maxi = ArbEnergyH.GetVectorLength();
  for (i = 0; i < maxi; ++i)
  {
    Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(i);
    Arb_y[i] = ArbEnergyH(i);
  }

  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy

  if (!DiffSpec)
  {
    // Converts integral point-wise spectra to Differential
    //
    for (count = 0; count < maxi - 1; ++count)
    {
      Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
                   / (Arb_x[count + 1] - Arb_x[count]);
    }
    --maxi;
  }

  if (!EnergySpec)
  {
    // Change currently stored values (emin etc) which are actually momenta
    // to energies
    //
    G4ParticleDefinition* pdef = threadLocalData.Get().particle_definition;
    if (pdef == nullptr)
    {
      G4Exception("G4SPSEneDistribution::ExpInterpolation",
                  "Event0302", FatalException,
                  "Error: particle not defined");
    }
    else
    {
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent

      G4double mass = pdef->GetPDGMass();

      // Convert point to energy unit and its value to per energy unit
      //
      G4double total_energy;
      for (count = 0; count < maxi; ++count)
      {
        total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass * mass));
        Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
        Arb_x[count] = total_energy - mass; // kinetic energy
      }
    }
  }

  i = 1;

  if ( Arb_ezero != nullptr ) { delete[] Arb_ezero; Arb_ezero = nullptr; }
  if ( Arb_Const != nullptr ) { delete[] Arb_Const; Arb_Const = nullptr; }
  Arb_ezero = new G4double [1024];
  Arb_Const = new G4double [1024];
  Arb_ezero_flag = true;
    
  Arb_ezero[0] = 0.;
  Arb_Const[0] = 0.;
  Area_seg[0] = 0.;
  Arb_Cum_Area[0] = 0.;
  while (i < maxi)
  {
    G4double test = std::log(Arb_y[i]) - std::log(Arb_y[i - 1]);
    if (test > 0. || test < 0.)
    {
      Arb_ezero[i] = -(Arb_x[i] - Arb_x[i - 1])
                   / (std::log(Arb_y[i]) - std::log(Arb_y[i - 1]));
      Arb_Const[i] = Arb_y[i] / (std::exp(-Arb_x[i] / Arb_ezero[i]));
      Area_seg[i] = -(Arb_Const[i] * Arb_ezero[i])
                  * (std::exp(-Arb_x[i] / Arb_ezero[i])
                   - std::exp(-Arb_x[i - 1] / Arb_ezero[i]));
    }
    else
    {
      G4Exception("G4SPSEneDistribution::ExpInterpolation",
                  "Event0302", JustWarning,
                  "Flat line segment: problem, setting to zero parameters.");
      G4cout << "Flat line segment: problem" << G4endl;
      Arb_ezero[i] = 0.;
      Arb_Const[i] = 0.;
      Area_seg[i] = 0.;
    }
    sum = sum + Area_seg[i];
    Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
    if (verbosityLevel == 2)
    {
      G4cout << Arb_ezero[i] << Arb_Const[i] << Area_seg[i] << G4endl;
    }
    ++i;
  }

  i = 0;
  while (i < maxi)
  {
    Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
    IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
    ++i;
  }

  // Now scale the ArbEnergyH, needed by Probability()
  //
  ArbEnergyH.ScaleVector(1., 1./sum);

  if (verbosityLevel >= 1)
  {
    G4cout << "Leaving ExpInterpolation " << G4endl;
  }
}

void G4SPSEneDistribution::SplineInterpolation()  // MT: Lock in caller
{
  // Interpolation using Splines.
  // Create Normalised arrays, make x 0->1 and y hold the function (Energy)
  // 
  // Current method based on the above will not work in all cases. 
  // New method is implemented below.
  
  G4double sum, Arb_x[1024]={0.}, Arb_y[1024]={0.}, Arb_Cum_Area[1024]={0.};
  std::size_t i, count;
  std::size_t maxi = ArbEnergyH.GetVectorLength();

  for (i = 0; i < maxi; ++i)
  {
    Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(i);
    Arb_y[i] = ArbEnergyH(i);
  }

  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy

  if (!DiffSpec)
  {
    // Converts integral point-wise spectra to Differential
    //
    for (count = 0; count < maxi - 1; ++count)
    {
      Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
                   / (Arb_x[count + 1] - Arb_x[count]);
    }
    --maxi;
  }

  if (!EnergySpec)
  {
    // Change currently stored values (emin etc) which are actually momenta
    // to energies
    //
    G4ParticleDefinition* pdef = threadLocalData.Get().particle_definition;
    if (pdef == nullptr)
    {
      G4Exception("G4SPSEneDistribution::SplineInterpolation",
                  "Event0302", FatalException,
                  "Error: particle not defined");
    }
    else
    {
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent

      G4double mass = pdef->GetPDGMass();

      // Convert point to energy unit and its value to per energy unit
      //
      G4double total_energy;
      for (count = 0; count < maxi; ++count)
      {
        total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass * mass));
        Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
        Arb_x[count] = total_energy - mass; // kinetic energy
      }
    }
  }

  i = 1;
  Arb_Cum_Area[0] = 0.;
  sum = 0.;
  Splinetemp = new G4DataInterpolation(Arb_x, Arb_y, (G4int)maxi, 0., 0.);
  G4double ei[101], prob[101];
  for (auto & it : SplineInt)
  {
    delete it;
    it = 0;
  }
  SplineInt.clear();
  SplineInt.resize(1024,nullptr);
  while (i < maxi)
  {
    // 100 step per segment for the integration of area

    G4double de = (Arb_x[i] - Arb_x[i - 1])/100.;
    G4double area = 0.;

    for (count = 0; count < 101; ++count)
    {
      ei[count] = Arb_x[i - 1] + de*count ;
      prob[count] =  Splinetemp->CubicSplineInterpolation(ei[count]);
      if (prob[count] < 0.)
      { 
        G4ExceptionDescription ED;
        ED << "Warning: G4DataInterpolation returns value < 0  " << prob[count]
           << " " << ei[count] << G4endl;
        G4Exception("G4SPSEneDistribution::SplineInterpolation", "Event0303",
                    FatalException, ED);
      }
      area += prob[count]*de;
    }
    Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + area;
    sum += area; 

    prob[0] = prob[0]/(area/de);
    for (count = 1; count < 100; ++count)
    {
      prob[count] = prob[count-1] + prob[count]/(area/de);
    }

    SplineInt[i] = new G4DataInterpolation(prob, ei, 101, 0., 0.);

    // NOTE: i starts from 1!
    //
    ++i;
  }

  i = 0;
  while (i < maxi)
  {
    Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
    IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
    ++i;
  }

  // Now scale the ArbEnergyH, needed by Probability()
  //
  ArbEnergyH.ScaleVector(1., 1./sum);

  if (verbosityLevel > 0)
  {
    G4cout << "Leaving SplineInterpolation " << G4endl;
  }
}

void G4SPSEneDistribution::GenerateMonoEnergetic()
{
  // Method to generate MonoEnergetic particles

  threadLocalData.Get().particle_energy = MonoEnergy;
}

void G4SPSEneDistribution::GenerateGaussEnergies()
{
  // Method to generate Gaussian particles

  G4double ene = G4RandGauss::shoot(MonoEnergy,SE);
  if (ene < 0) ene = 0.;
  threadLocalData.Get().particle_energy = ene;
}

void G4SPSEneDistribution::GenerateLinearEnergies(G4bool bArb = false)
{
  G4double rndm;
  threadLocal_t& params = threadLocalData.Get();
  G4double emaxsq = std::pow(params.Emax, 2.); // Emax squared
  G4double eminsq = std::pow(params.Emin, 2.); // Emin squared
  G4double intersq = std::pow(params.cept, 2.); // cept squared

  if (bArb) rndm = G4UniformRand();
  else      rndm = eneRndm->GenRandEnergy();

  G4double bracket = ((params.grad / 2.)
                   * (emaxsq - eminsq)
                   + params.cept * (params.Emax - params.Emin));
  bracket = bracket * rndm;
  bracket = bracket + (params.grad / 2.) * eminsq + params.cept * params.Emin;

  // Now have a quad of form m/2 E**2 + cE - bracket = 0
  //
  bracket = -bracket;

  if (params.grad != 0.)
  {
    G4double sqbrack = (intersq - 4 * (params.grad / 2.) * (bracket));
    sqbrack = std::sqrt(sqbrack);
    G4double root1 = -params.cept + sqbrack;
    root1 = root1 / (2. * (params.grad / 2.));

    G4double root2 = -params.cept - sqbrack;
    root2 = root2 / (2. * (params.grad / 2.));

    if (root1 > params.Emin && root1 < params.Emax)
    {
      params.particle_energy = root1;
    }
    if (root2 > params.Emin && root2 < params.Emax)
    {
      params.particle_energy = root2;
    }
  }
  else if (params.grad == 0.)
  {
    // have equation of form cE - bracket =0
    //
    params.particle_energy = bracket / params.cept;
  }

  if (params.particle_energy < 0.)
  {
    params.particle_energy = -params.particle_energy;
  }

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GeneratePowEnergies(G4bool bArb = false)
{
  // Method to generate particle energies distributed as a power-law

  G4double rndm;
  G4double emina, emaxa;
    
  threadLocal_t& params = threadLocalData.Get();
    
  emina = std::pow(params.Emin, params.alpha + 1);
  emaxa = std::pow(params.Emax, params.alpha + 1);

  if (bArb) rndm = G4UniformRand();
  else      rndm = eneRndm->GenRandEnergy();

  if (params.alpha != -1.)
  {
    G4double ene = ((rndm * (emaxa - emina)) + emina);
    ene = std::pow(ene, (1. / (params.alpha + 1.)));
    params.particle_energy = ene;
  }
  else
  {
    G4double ene = (std::log(params.Emin) 
                 + rndm * (std::log(params.Emax) - std::log(params.Emin)));
    params.particle_energy = std::exp(ene);
  }
  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenerateCPowEnergies()
{
  // Method to generate particle energies distributed in
  // cutoff power-law distribution
  //
  // CP_x holds Energies, and CPHist holds the cumulative histo.
  // binary search to find correct bin then lin interpolation.
  // Use the earlier defined histogram + RandGeneral method to generate
  // random numbers following the histos distribution

  G4double rndm = eneRndm->GenRandEnergy();
  G4int nabove = 10001, nbelow = 0, middle;

  G4AutoLock l(&mutex);
  G4bool done = CPhistCalcd;
  l.unlock();

  if(!done)
  {
    Calculate(); //This is has a lock inside, risk is to do it twice
    l.lock();
    CPhistCalcd = true;
    l.unlock();
  }

  // Binary search to find bin that rndm is in
  //
  while (nabove - nbelow > 1)
  {
    middle = (nabove + nbelow) / 2;
    if (rndm == CPHist->at(middle))
    {
      break;
    }
    if (rndm < CPHist->at(middle))
    {
      nabove = middle;
    }
    else
    {
      nbelow = middle;
    }
  }

  // Now interpolate in that bin to find the correct output value
  //
  G4double x1, x2, y1, y2, t, q;
  x1 = CP_x->at(nbelow);
  if(nbelow+1 == static_cast<G4int>(CP_x->size()))
  {
    x2 = CP_x->back();
  }
  else
  {
    x2 = CP_x->at(nbelow + 1);
  }
  y1 = CPHist->at(nbelow);
  if(nbelow+1 == static_cast<G4int>(CPHist->size()))
  {
    G4cout << CPHist->back() << G4endl;
    y2 = CPHist->back();
  }
  else
  {
    y2 = CPHist->at(nbelow + 1);
  }
  t = (y2 - y1) / (x2 - x1);
  q = y1 - t * x1;

  threadLocalData.Get().particle_energy = (rndm - q) / t;

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << threadLocalData.Get().particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenerateBiasPowEnergies()
{
  // Method to generate particle energies distributed as
  // in biased power-law and calculate its weight
    
  threadLocal_t& params = threadLocalData.Get();
    
  G4double rndm;
  G4double emina, emaxa, emin, emax;

  G4double normal = 1.;

  emin = params.Emin;
  emax = params.Emax;

  rndm = eneRndm->GenRandEnergy();

  if (biasalpha != -1.)
  {
    emina = std::pow(emin, biasalpha + 1);
    emaxa = std::pow(emax, biasalpha + 1);
    G4double ee = ((rndm * (emaxa - emina)) + emina);
    params.particle_energy = std::pow(ee, (1. / (biasalpha + 1.)));
    normal = 1./(1+biasalpha) * (emaxa - emina);
  }
  else
  {
    G4double ee = (std::log(emin) + rndm * (std::log(emax) - std::log(emin)));
    params.particle_energy = std::exp(ee);
    normal = std::log(emax) - std::log(emin);
  }
  params.weight = GetProbability(params.particle_energy)
                / (std::pow(params.particle_energy,biasalpha)/normal);

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenerateExpEnergies(G4bool bArb = false)
{
  // Method to generate particle energies distributed according
  // to an exponential curve

  G4double rndm;

  if (bArb) rndm = G4UniformRand();
  else      rndm = eneRndm->GenRandEnergy();

  threadLocal_t& params = threadLocalData.Get();
  params.particle_energy = -params.Ezero
                         * (std::log(rndm * (std::exp(-params.Emax
                                                     / params.Ezero)
                                           - std::exp(-params.Emin
                                                     / params.Ezero))
                                   + std::exp(-params.Emin / params.Ezero)));
  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy  << G4endl;
  }
}

void G4SPSEneDistribution::GenerateBremEnergies()
{
  // Method to generate particle energies distributed according
  // to a Bremstrahlung equation of the form
  // I = const*((kT)**1/2)*E*(e**(-E/kT))

  G4double rndm = eneRndm->GenRandEnergy();
  G4double expmax, expmin, k;
    
  k = 8.6181e-11; // Boltzmann's const in MeV/K
  G4double ksq = std::pow(k, 2.); // k squared
  G4double Tsq = std::pow(Temp, 2.); // Temp squared

  threadLocal_t& params = threadLocalData.Get();

  expmax = std::exp(-params.Emax / (k * Temp));
  expmin = std::exp(-params.Emin / (k * Temp));

  // If either expmax or expmin are zero then this will cause problems
  // Most probably this will be because T is too low or E is too high

  if (expmax == 0.)
  {
    G4Exception("G4SPSEneDistribution::GenerateBremEnergies",
                "Event0302", FatalException,
                "*****EXPMAX=0. Choose different E's or Temp");
  }
  if (expmin == 0.)
  {
    G4Exception("G4SPSEneDistribution::GenerateBremEnergies",
                "Event0302", FatalException,
                "*****EXPMIN=0. Choose different E's or Temp");
  }

  G4double tempvar = rndm * ((-k) * Temp * (params.Emax * expmax
                                          - params.Emin * expmin)
                          - (ksq * Tsq * (expmax - expmin)));

  G4double bigc = (tempvar - k * Temp * params.Emin * expmin
                 - ksq * Tsq * expmin) / (-k * Temp);

  // This gives an equation of form: Ee(-E/kT) + kTe(-E/kT) - C =0
  // Solve this iteratively, step from Emin to Emax in 1000 steps
  // and take the best solution.

  G4double erange = params.Emax - params.Emin;
  G4double steps = erange / 1000.;
  G4int i;
  G4double etest, diff, err = 100000.;

  for (i = 1; i < 1000; ++i)
  {
    etest = params.Emin + (i - 1) * steps;
    diff = etest * (std::exp(-etest / (k * Temp)))
         + k * Temp * (std::exp(-etest / (k * Temp))) - bigc;

    if (diff < 0.)
    {
      diff = -diff;
    }

    if (diff < err)
    {
      err = diff;
      params.particle_energy = etest;
    }
  }
  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenerateBbodyEnergies()
{
  // BBody_x holds Energies, and BBHist holds the cumulative histo.
  // Binary search to find correct bin then lin interpolation.
  // Use the earlier defined histogram + RandGeneral method to generate
  // random numbers following the histos distribution

  G4double rndm = eneRndm->GenRandEnergy();
  G4int nabove = 10001, nbelow = 0, middle;

  G4AutoLock l(&mutex);
  G4bool done = BBhistCalcd;
  l.unlock();

  if(!done)
  {
    Calculate(); //This is has a lock inside, risk is to do it twice
    l.lock();
    BBhistCalcd = true;
    l.unlock();
  }

  // Binary search to find bin that rndm is in
  //
  while (nabove - nbelow > 1)
  {
    middle = (nabove + nbelow) / 2;
    if (rndm == BBHist->at(middle))
    {
      break;
    }
    if (rndm < BBHist->at(middle))
    {
      nabove = middle;
    }
    else
    {
      nbelow = middle;
    }
  }

  // Now interpolate in that bin to find the correct output value
  //
  G4double x1, x2, y1, y2, t, q;
  x1 = Bbody_x->at(nbelow);

  if(nbelow+1 == static_cast<G4int>(Bbody_x->size()))
  {
    x2 = Bbody_x->back();
  }
  else
  {
    x2 = Bbody_x->at(nbelow + 1);
  }
  y1 = BBHist->at(nbelow);
  if(nbelow+1 == static_cast<G4int>(BBHist->size()))
  {
    G4cout << BBHist->back() << G4endl;
    y2 = BBHist->back();
  }
  else
  {
    y2 = BBHist->at(nbelow + 1);
  }
  t = (y2 - y1) / (x2 - x1);
  q = y1 - t * x1;

  threadLocalData.Get().particle_energy = (rndm - q) / t;

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << threadLocalData.Get().particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenerateCdgEnergies()
{
  // Generate random numbers, compare with values in cumhist
  // to find appropriate part of spectrum and then
  // generate energy in the usual inversion way

  G4double rndm, rndm2;
  G4double ene_line[3]={0,0,0};
  G4double omalpha[2]={0,0};
  threadLocal_t& params = threadLocalData.Get();
  if (params.Emin < 18 * keV && params.Emax < 18 * keV)
  {
    omalpha[0] = 1. - 1.4;
    ene_line[0] = params.Emin;
    ene_line[1] = params.Emax;
  }
  if (params.Emin < 18 * keV && params.Emax > 18 * keV)
  {
    omalpha[0] = 1. - 1.4;
    omalpha[1] = 1. - 2.3;
    ene_line[0] = params.Emin;
    ene_line[1] = 18. * keV;
    ene_line[2] = params.Emax;
  }
  if (params.Emin > 18 * keV)
  {
    omalpha[0] = 1. - 2.3;
    ene_line[0] = params.Emin;
    ene_line[1] = params.Emax;
  }
  rndm = eneRndm->GenRandEnergy();
  rndm2 = eneRndm->GenRandEnergy();

  G4int i = 0;
  while (rndm >= CDGhist[i])
  {
    ++i;
  }

  // Generate final energy
  //
  G4double ene = (std::pow(ene_line[i - 1], omalpha[i - 1])
               + (std::pow(ene_line[i], omalpha[i - 1])
                - std::pow(ene_line[i - 1], omalpha[i- 1])) * rndm2);
  params.particle_energy = std::pow(ene, (1. / omalpha[i - 1]));

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << params.particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::GenUserHistEnergies()
{
  // Histograms are DIFFERENTIAL

  G4AutoLock l(&mutex);

  if (!IPDFEnergyExist)
  {
    std::size_t ii;
    std::size_t maxbin = UDefEnergyH.GetVectorLength();
    G4double bins[1024], vals[1024], sum;
    for ( ii = 0 ; ii<1024 ; ++ii ) { bins[ii]=0; vals[ii]=0; }
    sum = 0.;

    if ( (!EnergySpec)
      && (threadLocalData.Get().particle_definition == nullptr))
    {
      G4Exception("G4SPSEneDistribution::GenUserHistEnergies",
                  "Event0302", FatalException,
                  "Error: particle definition is NULL");
    }

    if (maxbin > 1024)
    {
      G4Exception("G4SPSEneDistribution::GenUserHistEnergies",
                  "Event0302", JustWarning,
                 "Maxbin>1024\n Setting maxbin to 1024, other bins are lost");
      maxbin = 1024;
    }

    if (!DiffSpec)
    {
      G4cout << "Histograms are Differential!!! " << G4endl;
    }
    else
    {
      bins[0] = UDefEnergyH.GetLowEdgeEnergy(0);
      vals[0] = UDefEnergyH(0);
      sum = vals[0];
      for (ii = 1; ii < maxbin; ++ii)
      {
        bins[ii] = UDefEnergyH.GetLowEdgeEnergy(ii);
        vals[ii] = UDefEnergyH(ii) + vals[ii - 1];
        sum = sum + UDefEnergyH(ii);
      }
    }

    if (!EnergySpec)
    {
      G4double mass = threadLocalData.Get().particle_definition->GetPDGMass();

      // Multiply the function (vals) up by the bin width
      // to make the function counts/s (i.e. get rid of momentum dependence)

      for (ii = 1; ii < maxbin; ++ii)
      {
        vals[ii] = vals[ii] * (bins[ii] - bins[ii - 1]);
      }

      // Put energy bins into new histo, plus divide by energy bin width
      // to make evals counts/s/energy
      //
      for (ii = 0; ii < maxbin; ++ii)
      {
        // kinetic energy
        //
        bins[ii] = std::sqrt((bins[ii]*bins[ii])+(mass*mass))-mass;
      }
      for (ii = 1; ii < maxbin; ++ii)
      {
        vals[ii] = vals[ii] / (bins[ii] - bins[ii - 1]);
      }
      sum = vals[maxbin - 1];
      vals[0] = 0.;
    }
    for (ii = 0; ii < maxbin; ++ii)
    {
      vals[ii] = vals[ii] / sum;
      IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
    }

    IPDFEnergyExist = true;
    if (verbosityLevel > 1)
    {
      IPDFEnergyH.DumpValues();
    }
  }
  l.unlock();
    
  // IPDF has been create so carry on
  //
  G4double rndm = eneRndm->GenRandEnergy();
  threadLocalData.Get().particle_energy= IPDFEnergyH.GetEnergy(rndm);

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << particle_energy << G4endl;
  }
}

G4double G4SPSEneDistribution::GetArbEneWeight(G4double ene)
{
  auto nbelow = IPDFArbEnergyH.FindBin(ene,(IPDFArbEnergyH.GetVectorLength())/2);
  G4double wei = 0.;
  if(IntType=="Lin")
  {
    // note: grad[i] and cept[i] are calculated with x[i-1] and x[i]
    auto gr = Arb_grad[nbelow + 1];
    auto ce = Arb_cept[nbelow + 1];
    wei = ene*gr + ce;
  }
  else if(IntType=="Log")
  {
    auto alp = Arb_alpha[nbelow + 1];
    auto cns = Arb_Const[nbelow + 1];
    wei = cns * std::pow(ene,alp);
  }
  else if(IntType=="Exp")
  {
    auto e0 = Arb_ezero[nbelow + 1];
    auto cns = Arb_Const[nbelow + 1];
    wei = cns * std::exp(-ene/e0);
  }
  else if(IntType=="Spline")
  {
    wei = SplineInt[nbelow+1]->CubicSplineInterpolation(ene);
  }
  return wei;
}

void G4SPSEneDistribution::GenArbPointEnergies()
{
  if (verbosityLevel > 0)
  {
    G4cout << "In GenArbPointEnergies" << G4endl;
  }

  G4double rndm = eneRndm->GenRandEnergy();

  // Find the Bin, have x, y, no of points, and cumulative area distribution
  //
  std::size_t nabove = IPDFArbEnergyH.GetVectorLength(), nbelow = 0, middle;

  // Binary search to find bin that rndm is in
  //
  while (nabove - nbelow > 1)
  {
    middle = (nabove + nbelow) / 2;
    if (rndm == IPDFArbEnergyH(middle))
    {
      break;
    }
    if (rndm < IPDFArbEnergyH(middle))
    {
      nabove = middle;
    }
    else
    {
      nbelow = middle;
    }
  }
  threadLocal_t& params = threadLocalData.Get();
  if (IntType == "Lin")
  {
    // Update thread-local copy of parameters
    //
    params.Emax = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow + 1);
    params.Emin = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow);
    params.grad = Arb_grad[nbelow + 1];
    params.cept = Arb_cept[nbelow + 1];
    GenerateLinearEnergies(true);
  }
  else if (IntType == "Log")
  {
    params.Emax = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow + 1);
    params.Emin = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow);
    params.alpha = Arb_alpha[nbelow + 1];
    GeneratePowEnergies(true);
  }
  else if (IntType == "Exp")
  {
    params.Emax = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow + 1);
    params.Emin = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow);
    params.Ezero = Arb_ezero[nbelow + 1];
    GenerateExpEnergies(true);
  }
  else if (IntType == "Spline")
  {
    params.Emax = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow + 1);
    params.Emin = IPDFArbEnergyH.GetLowEdgeEnergy(nbelow);
    params.particle_energy = -1e100;
    rndm = eneRndm->GenRandEnergy();
    while (params.particle_energy < params.Emin
        || params.particle_energy > params.Emax)
    {
      params.particle_energy =
        SplineInt[nbelow+1]->CubicSplineInterpolation(rndm);
      rndm = eneRndm->GenRandEnergy();
    }
    if (verbosityLevel >= 1)
    {
      G4cout << "Energy is " << params.particle_energy << G4endl;
    }
  }
  else
  {
    G4Exception("G4SPSEneDistribution::GenArbPointEnergies", "Event0302",
                FatalException, "Error: IntType unknown type");
  }
}

void G4SPSEneDistribution::GenEpnHistEnergies()
{
  // Firstly convert to energy if not already done

  G4AutoLock l(&mutex);

  if (Epnflag)  // true means spectrum is epn, false means e
  {
    // Convert to energy by multiplying by A number
    //
    ConvertEPNToEnergy();
  }
  if (!IPDFEnergyExist)
  {
    // IPDF has not been created, so create it
    //
    G4double bins[1024], vals[1024], sum;
    std::size_t ii;
    std::size_t maxbin = UDefEnergyH.GetVectorLength();
    bins[0] = UDefEnergyH.GetLowEdgeEnergy(0);
    vals[0] = UDefEnergyH(0);
    sum = vals[0];
    for (ii = 1; ii < maxbin; ++ii)
    {
      bins[ii] = UDefEnergyH.GetLowEdgeEnergy(ii);
      vals[ii] = UDefEnergyH(ii) + vals[ii - 1];
      sum = sum + UDefEnergyH(ii);
    }

    l.lock();
    for (ii = 0; ii < maxbin; ++ii)
    {
      vals[ii] = vals[ii] / sum;
      IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
    }
    IPDFEnergyExist = true;
       
  }
  l.unlock();

  // IPDF has been create so carry on
  //
  G4double rndm = eneRndm->GenRandEnergy();
  threadLocalData.Get().particle_energy = IPDFEnergyH.GetEnergy(rndm);

  if (verbosityLevel >= 1)
  {
    G4cout << "Energy is " << threadLocalData.Get().particle_energy << G4endl;
  }
}

void G4SPSEneDistribution::ConvertEPNToEnergy()  // MT: lock in caller
{
  // Use this before particle generation to convert the
  // currently stored histogram from energy/nucleon to energy.

  threadLocal_t& params = threadLocalData.Get();
  if (params.particle_definition == nullptr)
  {
    G4cout << "Error: particle not defined" << G4endl;
  }
  else
  {
    // Need to multiply histogram by the number of nucleons.
    // Baryon Number looks to hold the no. of nucleons
    //
    G4int Bary = params.particle_definition->GetBaryonNumber();

    // Change values in histogram, Read it out, delete it, re-create it
    //
    std::size_t count, maxcount;
    maxcount = EpnEnergyH.GetVectorLength();
    G4double ebins[1024], evals[1024];
    if (maxcount > 1024)
    {
      G4Exception("G4SPSEneDistribution::ConvertEPNToEnergy()",
                  "gps001", JustWarning,
                  "Histogram contains more than 1024 bins!\n\
                   Those above 1024 will be ignored");
      maxcount = 1024;
    }
    if (maxcount < 1)
    {
      G4Exception("G4SPSEneDistribution::ConvertEPNToEnergy()",
                 "gps001", FatalException,
                 "Histogram contains less than 1 bin!\nRedefine the histogram");
      return;
    }
    for (count = 0; count < maxcount; ++count)
    {
      // Read out
      ebins[count] = EpnEnergyH.GetLowEdgeEnergy(count);
      evals[count] = EpnEnergyH(count);
    }

    // Multiply the channels by the nucleon number to give energies
    //
    for (count = 0; count < maxcount; ++count)
    {
      ebins[count] = ebins[count] * Bary;
    }

    // Set Emin and Emax
    //
    params.Emin = ebins[0];
    if (maxcount > 1)
    {
      params.Emax = ebins[maxcount - 1];
    }
    else
    {
      params.Emax = ebins[0];
    }

    // Put energy bins into new histogram - UDefEnergyH
    //
    for (count = 0; count < maxcount; ++count)
    {
      UDefEnergyH.InsertValues(ebins[count], evals[count]);
    }
    Epnflag = false; // so that you dont repeat this method
  }
}

void G4SPSEneDistribution::ReSetHist(const G4String& atype)
{
  G4AutoLock l(&mutex);
  if (atype == "energy")
  {
    UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
    IPDFEnergyExist = false;
    Emin = 0.;
    Emax = 1e30;
  }
  else if (atype == "arb")
  {
    ArbEnergyH = IPDFArbEnergyH = ZeroPhysVector;
    IPDFArbExist = false;
  }
  else if (atype == "epn")
  {
    UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
    IPDFEnergyExist = false;
    EpnEnergyH = ZeroPhysVector;
  }
  else
  {
    G4cout << "Error, histtype not accepted " << G4endl;
  }
}

G4double G4SPSEneDistribution::GenerateOne(G4ParticleDefinition* a)
{
  // Copy global shared status to thread-local one
  //
  threadLocal_t& params = threadLocalData.Get();
  params.particle_definition=a;
  params.particle_energy=-1;
  if(applyEvergyWeight)
  {
    params.Emax = ArbEmax;
    params.Emin = ArbEmin;
  }
  else
  {
    params.Emax = Emax;
    params.Emin = Emin;
  }
  params.alpha = alpha;
  params.Ezero = Ezero;
  params.grad = grad;
  params.cept = cept;
  params.weight = weight;
  // particle_energy = -1.;

  if((EnergyDisType == "Mono") && ((MonoEnergy>Emax)||(MonoEnergy<Emin)))
  {
    G4ExceptionDescription ed;
    ed << "MonoEnergy " << G4BestUnit(MonoEnergy,"Energy")
       << " is outside of [Emin,Emax] = ["
       << G4BestUnit(Emin,"Energy") << ", "
       << G4BestUnit(Emax,"Energy") << ". MonoEnergy is used anyway.";
    G4Exception("G4SPSEneDistribution::GenerateOne()",
                "GPS0001", JustWarning, ed);
    params.particle_energy=MonoEnergy;
    return params.particle_energy;
  }
  while ( (EnergyDisType == "Arb")
        ?   (params.particle_energy < ArbEmin
          || params.particle_energy > ArbEmax)
        :   (params.particle_energy < params.Emin
          || params.particle_energy > params.Emax) )
  {
    if (Biased)
    {
      GenerateBiasPowEnergies();
    }
    else
    {
      if (EnergyDisType == "Mono")
      {
        GenerateMonoEnergetic();
      }
      else if (EnergyDisType == "Lin")
      {
        GenerateLinearEnergies(false);
      }
      else if (EnergyDisType == "Pow")
      {
        GeneratePowEnergies(false);
      }
      else if (EnergyDisType == "CPow")
      {
        GenerateCPowEnergies();
      }
      else if (EnergyDisType == "Exp")
      {
        GenerateExpEnergies(false);
      }
      else if (EnergyDisType == "Gauss")
      {
        GenerateGaussEnergies();
      }
      else if (EnergyDisType == "Brem")
      {
        GenerateBremEnergies();
      }
      else if (EnergyDisType == "Bbody")
      {
        GenerateBbodyEnergies();
      }
      else if (EnergyDisType == "Cdg")
      {
        GenerateCdgEnergies();
      }
      else if (EnergyDisType == "User")
      {
        GenUserHistEnergies();
      }
      else if (EnergyDisType == "Arb")
      {
        GenArbPointEnergies();
      }
      else if (EnergyDisType == "Epn")
      {
        GenEpnHistEnergies();
      }
      else
      {
        G4cout << "Error: EnergyDisType has unusual value" << G4endl;
      }
    }
  }
   return params.particle_energy;
}

G4double G4SPSEneDistribution::GetProbability(G4double ene)
{
  G4double prob = 1.;

  threadLocal_t& params = threadLocalData.Get();
  if (EnergyDisType == "Lin")
  {
    if (prob_norm == 1.)
    {
      prob_norm = 0.5*params.grad*params.Emax*params.Emax
                + params.cept*params.Emax
                - 0.5*params.grad*params.Emin*params.Emin
                - params.cept*params.Emin;
    }
    prob = params.cept + params.grad * ene;
    prob /= prob_norm;
  }
  else if (EnergyDisType == "Pow")
  {
    if (prob_norm == 1.)
    {
      if (alpha != -1.)
      {
        G4double emina = std::pow(params.Emin, params.alpha + 1);
        G4double emaxa = std::pow(params.Emax, params.alpha + 1);
        prob_norm = 1./(1.+alpha) * (emaxa - emina);
      }
      else
      {
        prob_norm = std::log(params.Emax) - std::log(params.Emin) ;
      }
    }
    prob = std::pow(ene, params.alpha)/prob_norm;
  }
  else if (EnergyDisType == "Exp")
  {
    if (prob_norm == 1.)
    {
      prob_norm = -params.Ezero*(std::exp(-params.Emax/params.Ezero)
                               - std::exp(params.Emin/params.Ezero));
    }  
    prob = std::exp(-ene / params.Ezero);
    prob /= prob_norm;
  }
  else if (EnergyDisType == "Arb")
  {
    prob = ArbEnergyH.Value(ene);

    if (prob <= 0.)
    {
      G4cout << " Warning:G4SPSEneDistribution::GetProbability: prob<= 0. "
             << prob << " " << ene << G4endl;
      prob = 1e-30;
    }
  }
  else
  {
    G4cout << "Error: EnergyDisType not supported" << G4endl;
  }

  return prob;
}
