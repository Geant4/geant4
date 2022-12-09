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
// History:
// 2001-2002 R&D by V.Grichine
// 19.06.03 V. Grichine, modifications in BuildTable for the integration
//                       in respect of angle: range is increased, accuracy is
//                       improved
// 28.07.05, P.Gumplinger add G4ProcessType to constructor
// 28.09.07, V.Ivanchenko general cleanup without change of algorithms
//

#include "G4VXTRenergyLoss.hh"

#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4EmProcessSubType.hh"
#include "G4Integrator.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SandiaTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Timer.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"
#include "G4PhysicsModelCatalog.hh"

////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4VXTRenergyLoss::G4VXTRenergyLoss(G4LogicalVolume* anEnvelope,
                                   G4Material* foilMat, G4Material* gasMat,
                                   G4double a, G4double b, G4int n,
                                   const G4String& processName,
                                   G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
  , fGammaCutInKineticEnergy(nullptr)
  , fAngleDistrTable(nullptr)
  , fEnergyDistrTable(nullptr)
  , fAngleForEnergyTable(nullptr)
  , fPlatePhotoAbsCof(nullptr)
  , fGasPhotoAbsCof(nullptr)
  , fGammaTkinCut(0.0)
{
  verboseLevel = 1;
  secID = G4PhysicsModelCatalog::GetModelID("model_XTRenergyLoss");
  SetProcessSubType(fTransitionRadiation);

  fPtrGamma    = nullptr;
  fMinEnergyTR = fMaxEnergyTR = fMaxThetaTR = fGamma = fEnergy = 0.0;
  fVarAngle = fLambda = fTotalDist = fPlateThick = fGasThick = 0.0;
  fAlphaPlate = 100.;
  fAlphaGas = 40.;

  fTheMinEnergyTR = CLHEP::keV * 1.; //  1.; // 
  fTheMaxEnergyTR = CLHEP::keV * 100.; // 40.; //

  fTheMinAngle    = 1.e-8;  //
  fTheMaxAngle    = 4.e-4;

  fTotBin = 50;  //  number of bins in log scale 
  fBinTR  =  100; //   number of bins in TR vectors

  // min/max angle2 in log-vectors

  fMinThetaTR = 3.0e-9; 
  fMaxThetaTR = 1.0e-4;

  
  // Proton energy vector initialization
  fProtonEnergyVector =
    new G4PhysicsLogVector(fMinProtonTkin, fMaxProtonTkin, fTotBin);

  fXTREnergyVector =
    new G4PhysicsLogVector(fTheMinEnergyTR, fTheMaxEnergyTR, fBinTR);

  fEnvelope = anEnvelope;

  fPlateNumber = n;
  if(verboseLevel > 0)
    G4cout << "### G4VXTRenergyLoss: the number of TR radiator plates = "
           << fPlateNumber << G4endl;
  if(fPlateNumber == 0)
  {
    G4Exception("G4VXTRenergyLoss::G4VXTRenergyLoss()", "VXTRELoss01",
                FatalException, "No plates in X-ray TR radiator");
  }
  // default is XTR dEdx, not flux after radiator
  fExitFlux      = false;
  // default angle distribution according numerical integration
  fFastAngle     = false; // no angle according sum of delta-functions by default
  fAngleRadDistr = true;
  fCompton       = false;

  fLambda = DBL_MAX;

  // Mean thicknesses of plates and gas gaps
  fPlateThick = a;
  fGasThick   = b;
  fTotalDist  = fPlateNumber * (fPlateThick + fGasThick);
  if(verboseLevel > 0)
    G4cout << "total radiator thickness = " << fTotalDist / cm << " cm"
           << G4endl;

  // index of plate material
  fMatIndex1 = (G4int)foilMat->GetIndex();
  if(verboseLevel > 0)
    G4cout << "plate material = " << foilMat->GetName() << G4endl;

  // index of gas material
  fMatIndex2 = (G4int)gasMat->GetIndex();
  if(verboseLevel > 0)
    G4cout << "gas material = " << gasMat->GetName() << G4endl;

  // plasma energy squared for plate material
  fSigma1 = fPlasmaCof * foilMat->GetElectronDensity();
  if(verboseLevel > 0)
    G4cout << "plate plasma energy = " << std::sqrt(fSigma1) / eV << " eV"
           << G4endl;

  // plasma energy squared for gas material
  fSigma2 = fPlasmaCof * gasMat->GetElectronDensity();
  if(verboseLevel > 0)
    G4cout << "gas plasma energy = " << std::sqrt(fSigma2) / eV << " eV"
           << G4endl;

  // Compute cofs for preparation of linear photo absorption
  ComputePlatePhotoAbsCof();
  ComputeGasPhotoAbsCof();

  pParticleChange = &fParticleChange;
}

///////////////////////////////////////////////////////////////////////////
G4VXTRenergyLoss::~G4VXTRenergyLoss()
{
  delete fProtonEnergyVector;
  delete fXTREnergyVector;
  if(fEnergyDistrTable)
  {
    fEnergyDistrTable->clearAndDestroy();
    delete fEnergyDistrTable;
  }
  if(fAngleRadDistr)
  {
    fAngleDistrTable->clearAndDestroy();
    delete fAngleDistrTable;
  }
  if(fAngleForEnergyTable)
  {
    fAngleForEnergyTable->clearAndDestroy();
    delete fAngleForEnergyTable;
  }
}

void G4VXTRenergyLoss::ProcessDescription(std::ostream& out) const
{
  out << "Base class for 'fast' parameterisation model describing X-ray "
         "transition\n"
         "radiation. Angular distribution is very rough.\n";
}

///////////////////////////////////////////////////////////////////////////////
// Returns condition for application of the model depending on particle type
G4bool G4VXTRenergyLoss::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.0);
}

/////////////////////////////////////////////////////////////////////////////////
// Calculate step size for XTR process inside raaditor
G4double G4VXTRenergyLoss::GetMeanFreePath(const G4Track& aTrack, G4double,
                                           G4ForceCondition* condition)
{
  G4int iTkin, iPlace;
  G4double lambda, sigma, kinEnergy, mass, gamma;
  G4double charge, chargeSq, massRatio, TkinScaled;
  G4double E1, E2, W, W1, W2;

  *condition = NotForced;

  if(aTrack.GetVolume()->GetLogicalVolume() != fEnvelope)
    lambda = DBL_MAX;
  else
  {
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    kinEnergy                          = aParticle->GetKineticEnergy();
    mass  = aParticle->GetDefinition()->GetPDGMass();
    gamma = 1.0 + kinEnergy / mass;
    if(verboseLevel > 1)
    {
      G4cout << " gamma = " << gamma << ";   fGamma = " << fGamma << G4endl;
    }

    if(std::fabs(gamma - fGamma) < 0.05 * gamma)
      lambda = fLambda;
    else
    {
      charge     = aParticle->GetDefinition()->GetPDGCharge();
      chargeSq   = charge * charge;
      massRatio  = proton_mass_c2 / mass;
      TkinScaled = kinEnergy * massRatio;

      for(iTkin = 0; iTkin < fTotBin; ++iTkin)
      {
        if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))
          break;
      }
      iPlace = iTkin - 1;

      if(iTkin == 0)
        lambda = DBL_MAX;  // Tkin is too small, neglect of TR photon generation
      else  // general case: Tkin between two vectors of the material
      {
        if(iTkin == fTotBin)
        {
          sigma = (*(*fEnergyDistrTable)(iPlace))(0) * chargeSq;
        }
        else
        {
          E1    = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1);
          E2    = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
          W     = 1.0 / (E2 - E1);
          W1    = (E2 - TkinScaled) * W;
          W2    = (TkinScaled - E1) * W;
          sigma = ((*(*fEnergyDistrTable)(iPlace))(0) * W1 +
                   (*(*fEnergyDistrTable)(iPlace + 1))(0) * W2) *
                  chargeSq;
        }
        if(sigma < DBL_MIN)
          lambda = DBL_MAX;
        else
          lambda = 1. / sigma;
        fLambda = lambda;
        fGamma  = gamma;
        if(verboseLevel > 1)
        {
          G4cout << " lambda = " << lambda / mm << " mm" << G4endl;
        }
      }
    }
  }
  return lambda;
}

//////////////////////////////////////////////////////////////////////////
// Interface for build table from physics list
void G4VXTRenergyLoss::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  if(pd.GetPDGCharge() == 0.)
  {
    G4Exception("G4VXTRenergyLoss::BuildPhysicsTable", "Notification",
                JustWarning, "XTR initialisation for neutral particle ?!");
  }
  BuildEnergyTable();

  if(fAngleRadDistr)
  {
    if(verboseLevel > 0)
    {
      G4cout
        << "Build angle for energy distribution according the current radiator"
        << G4endl;
    }
    BuildAngleForEnergyBank();
  }
}

//////////////////////////////////////////////////////////////////////////
// Build integral energy distribution of XTR photons
void G4VXTRenergyLoss::BuildEnergyTable()
{
  G4int iTkin, iTR, iPlace;
  G4double radiatorCof = 1.0;  // for tuning of XTR yield
  G4double energySum   = 0.0;

  fEnergyDistrTable = new G4PhysicsTable(fTotBin);
  if(fAngleRadDistr)
    fAngleDistrTable = new G4PhysicsTable(fTotBin);

  fGammaTkinCut = 0.0;

  // setting of min/max TR energies
  if(fGammaTkinCut > fTheMinEnergyTR)
    fMinEnergyTR = fGammaTkinCut;
  else
    fMinEnergyTR = fTheMinEnergyTR;

  if(fGammaTkinCut > fTheMaxEnergyTR)
    fMaxEnergyTR = 2.0 * fGammaTkinCut;
  else
    fMaxEnergyTR = fTheMaxEnergyTR;

  G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
    integral;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();

  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "Lorentz Factor"
           << "\t"
           << "XTR photon number" << G4endl;
    G4cout << G4endl;
  }
  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    auto energyVector =
      new G4PhysicsLogVector(fMinEnergyTR, fMaxEnergyTR, fBinTR);

    fGamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);

    // if(fMaxThetaTR > fTheMaxAngle)     fMaxThetaTR = fTheMaxAngle;
    // else if(fMaxThetaTR < fTheMinAngle)     fMaxThetaTR = fTheMinAngle;

    energySum = 0.0;

    energyVector->PutValue(fBinTR - 1, energySum);

    for(iTR = fBinTR - 2; iTR >= 0; --iTR)
    {
      // Legendre96 or Legendre10

      energySum += radiatorCof * fCofTR *
	
	// integral.Legendre10(this, &G4VXTRenergyLoss::SpectralXTRdEdx,
	
                   integral.Legendre96(this, &G4VXTRenergyLoss::SpectralXTRdEdx,
				       
                                       energyVector->GetLowEdgeEnergy(iTR),
                                       energyVector->GetLowEdgeEnergy(iTR + 1));

      energyVector->PutValue(iTR, energySum / fTotalDist);
    }
    iPlace = iTkin;
    fEnergyDistrTable->insertAt(iPlace, energyVector);

    if(verboseLevel > 0)
    {
      G4cout << fGamma << "\t" << energySum << G4endl;
    }
  }
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "total time for build X-ray TR energy loss tables = "
           << timer.GetUserElapsed() << " s" << G4endl;
  }
  fGamma = 0.;
  return;
}

//////////////////////////////////////////////////////////////////////////
// Bank of angle distributions for given energies (slow!)

void G4VXTRenergyLoss::BuildAngleForEnergyBank()
{
  
  if( ( this->GetProcessName() == "TranspRegXTRadiator" ||
        this->GetProcessName() == "TranspRegXTRmodel" ||
        this->GetProcessName() == "RegularXTRadiator" ||
	this->GetProcessName() == "RegularXTRmodel"  )       && fFastAngle    ) // ffastAngle=true!
  {
    BuildAngleTable(); // by sum of delta-functions
    return;
  }
  G4int i, iTkin, iTR;
  G4double angleSum = 0.0;

  fGammaTkinCut = 0.0;

  // setting of min/max TR energies
  if(fGammaTkinCut > fTheMinEnergyTR)
    fMinEnergyTR = fGammaTkinCut;
  else
    fMinEnergyTR = fTheMinEnergyTR;

  if(fGammaTkinCut > fTheMaxEnergyTR)
    fMaxEnergyTR = 2.0 * fGammaTkinCut;
  else
    fMaxEnergyTR = fTheMaxEnergyTR;

  auto energyVector =
    new G4PhysicsLogVector(fMinEnergyTR, fMaxEnergyTR, fBinTR);

  G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
    integral;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();

  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    fGamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);

    if(fMaxThetaTR > fTheMaxAngle)
      fMaxThetaTR = fTheMaxAngle;
    else if(fMaxThetaTR < fTheMinAngle)
      fMaxThetaTR = fTheMinAngle;

    fAngleForEnergyTable = new G4PhysicsTable(fBinTR);

    for(iTR = 0; iTR < fBinTR; ++iTR)
    {
      angleSum = 0.0;
      fEnergy  = energyVector->GetLowEdgeEnergy(iTR);
      
     // log-vector to increase number of thin bins for small angles
      auto angleVector = new G4PhysicsLogVector(fMinThetaTR, fMaxThetaTR, fBinTR);
 
      

      angleVector->PutValue(fBinTR - 1, angleSum);

      for(i = fBinTR - 2; i >= 0; --i)
      {
        // Legendre96 or Legendre10

        angleSum +=
          integral.Legendre10(this, &G4VXTRenergyLoss::SpectralAngleXTRdEdx,
                              angleVector->GetLowEdgeEnergy(i),
                              angleVector->GetLowEdgeEnergy(i + 1));

        angleVector->PutValue(i, angleSum);
      }
      fAngleForEnergyTable->insertAt(iTR, angleVector);
    }
    fAngleBank.push_back(fAngleForEnergyTable);
  }
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "total time for build X-ray TR angle for energy loss tables = "
           << timer.GetUserElapsed() << " s" << G4endl;
  }
  fGamma = 0.;
  delete energyVector;
}

////////////////////////////////////////////////////////////////////////
// Build XTR angular distribution at given energy based on the model
// of transparent regular radiator
void G4VXTRenergyLoss::BuildAngleTable()
{
  G4int iTkin, iTR;
  G4double energy;

  fGammaTkinCut = 0.0;

  // setting of min/max TR energies
  if(fGammaTkinCut > fTheMinEnergyTR)
    fMinEnergyTR = fGammaTkinCut;
  else
    fMinEnergyTR = fTheMinEnergyTR;

  if(fGammaTkinCut > fTheMaxEnergyTR)
    fMaxEnergyTR = 2.0 * fGammaTkinCut;
  else
    fMaxEnergyTR = fTheMaxEnergyTR;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();
  if(verboseLevel > 0)
  {
    G4cout << G4endl << "Lorentz Factor" << "\t"
           << "XTR photon number" << G4endl << G4endl;
  }
  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    fGamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);

    // fMaxThetaTR = 25. * 2500.0 / (fGamma * fGamma);  // theta^2

    if(fMaxThetaTR > fTheMaxAngle)
      fMaxThetaTR = fTheMaxAngle;
    else
    {
      if(fMaxThetaTR < fTheMinAngle)
        fMaxThetaTR = fTheMinAngle;
    }

    fAngleForEnergyTable = new G4PhysicsTable(fBinTR);

    for(iTR = 0; iTR < fBinTR; ++iTR)
    {
      energy = fXTREnergyVector->GetLowEdgeEnergy(iTR);

      G4PhysicsFreeVector* angleVector = GetAngleVector(energy, fBinTR);

      fAngleForEnergyTable->insertAt(iTR, angleVector);
    }
    fAngleBank.push_back(fAngleForEnergyTable);
  }
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "total time for build XTR angle for given energy tables = "
           << timer.GetUserElapsed() << " s" << G4endl;
  }
  fGamma = 0.;

  return;
}

/////////////////////////////////////////////////////////////////////////
// Vector of angles and angle integral distributions
G4PhysicsFreeVector* G4VXTRenergyLoss::GetAngleVector(G4double energy, G4int n)
{
  G4double theta = 0., result, tmp = 0., cof1, cof2, cofMin, cofPHC,
           angleSum = 0.;
  G4int iTheta, k, kMin;

  auto angleVector = new G4PhysicsFreeVector(n);

  cofPHC = 4. * pi * hbarc;
  tmp    = (fSigma1 - fSigma2) / cofPHC / energy;
  cof1   = fPlateThick * tmp;
  cof2   = fGasThick * tmp;

  cofMin = energy * (fPlateThick + fGasThick) / fGamma / fGamma;
  cofMin += (fPlateThick * fSigma1 + fGasThick * fSigma2) / energy;
  cofMin /= cofPHC;

  kMin = G4int(cofMin);
  if(cofMin > kMin)
    kMin++;

  if(verboseLevel > 2)
  {
    G4cout << "n-1 = " << n - 1
           << "; theta = " << std::sqrt(fMaxThetaTR) * fGamma
           << "; tmp = " << 0. << ";    angleSum = " << angleSum << G4endl;
  }

  for(iTheta = n - 1; iTheta >= 1; --iTheta)
  {
    k      = iTheta - 1 + kMin;
    tmp    = pi * fPlateThick * (k + cof2) / (fPlateThick + fGasThick);
    result = (k - cof1) * (k - cof1) * (k + cof2) * (k + cof2);
    tmp    = std::sin(tmp) * std::sin(tmp) * std::abs(k - cofMin) / result;

    if(k == kMin && kMin == G4int(cofMin))
    {
      // angleSum += 0.5 * tmp;
      angleSum += tmp; // ATLAS TB 
    }
    else if(iTheta == n - 1)
      ;
    else
    {
      angleSum += tmp;
    }
    theta = std::abs(k - cofMin) * cofPHC / energy / (fPlateThick + fGasThick);

    if(verboseLevel > 2)
    {
      G4cout << "iTheta = " << iTheta << "; k = " << k
             << "; theta = " << std::sqrt(theta) * fGamma << "; tmp = " << tmp
             << ";    angleSum = " << angleSum << G4endl;
    }
    angleVector->PutValue(iTheta, theta, angleSum);
  }
  if(theta > 0.)
  {
    // angleSum += 0.5 * tmp;
    angleSum += 0.;  // ATLAS TB
    theta     = 0.;
  }
  if(verboseLevel > 2)
  {
    G4cout << "iTheta = " << iTheta << "; theta = " << std::sqrt(theta) * fGamma
           << "; tmp = " << tmp << ";    angleSum = " << angleSum << G4endl;
  }
  angleVector->PutValue(iTheta, theta, angleSum);

  return angleVector;
}

////////////////////////////////////////////////////////////////////////
// Build XTR angular distribution based on the model of transparent regular
// radiator
void G4VXTRenergyLoss::BuildGlobalAngleTable()
{
  G4int iTkin, iTR, iPlace;
  G4double radiatorCof = 1.0;  // for tuning of XTR yield
  G4double angleSum;
  fAngleDistrTable = new G4PhysicsTable(fTotBin);

  fGammaTkinCut = 0.0;

  // setting of min/max TR energies
  if(fGammaTkinCut > fTheMinEnergyTR)
    fMinEnergyTR = fGammaTkinCut;
  else
    fMinEnergyTR = fTheMinEnergyTR;

  if(fGammaTkinCut > fTheMaxEnergyTR)
    fMaxEnergyTR = 2.0 * fGammaTkinCut;
  else
    fMaxEnergyTR = fTheMaxEnergyTR;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "Lorentz Factor"
           << "\t"
           << "XTR photon number" << G4endl;
    G4cout << G4endl;
  }
  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    fGamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);

    // fMaxThetaTR = 25.0 / (fGamma * fGamma);  // theta^2
    // fMaxThetaTR = 1.e-4;  // theta^2

    if(fMaxThetaTR > fTheMaxAngle)
      fMaxThetaTR = fTheMaxAngle;
    else
    {
      if(fMaxThetaTR < fTheMinAngle)
        fMaxThetaTR = fTheMinAngle;
    }
    auto angleVector =
    // G4PhysicsLogVector* angleVector =
      new G4PhysicsLinearVector(0.0, fMaxThetaTR, fBinTR);
    //  new G4PhysicsLogVector(1.e-8, fMaxThetaTR, fBinTR);

    angleSum = 0.0;

    G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
      integral;

    angleVector->PutValue(fBinTR - 1, angleSum);

    for(iTR = fBinTR - 2; iTR >= 0; --iTR)
    {
      angleSum += radiatorCof * fCofTR *
                  integral.Legendre96(this, &G4VXTRenergyLoss::AngleXTRdEdx,
                                      angleVector->GetLowEdgeEnergy(iTR),
                                      angleVector->GetLowEdgeEnergy(iTR + 1));

      angleVector->PutValue(iTR, angleSum);
    }
    if(verboseLevel > 1)
    {
      G4cout << fGamma << "\t" << angleSum << G4endl;
    }
    iPlace = iTkin;
    fAngleDistrTable->insertAt(iPlace, angleVector);
  }
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "total time for build X-ray TR angle tables = "
           << timer.GetUserElapsed() << " s" << G4endl;
  }
  fGamma = 0.;

  return;
}

//////////////////////////////////////////////////////////////////////////////
// The main function which is responsible for the treatment of a particle
// passage through G4Envelope with discrete generation of G4Gamma
G4VParticleChange* G4VXTRenergyLoss::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step& aStep)
{
  G4int iTkin;
  G4double energyTR, theta, theta2, phi, dirX, dirY, dirZ;

  fParticleChange.Initialize(aTrack);

  if(verboseLevel > 1)
  {
    G4cout << "Start of G4VXTRenergyLoss::PostStepDoIt " << G4endl;
    G4cout << "name of current material =  "
           << aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()
           << G4endl;
  }
  if(aTrack.GetVolume()->GetLogicalVolume() != fEnvelope)
  {
    if(verboseLevel > 0)
    {
      G4cout << "Go out from G4VXTRenergyLoss::PostStepDoIt: wrong volume "
             << G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  else
  {
    G4StepPoint* pPostStepPoint        = aStep.GetPostStepPoint();
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    // Now we are ready to Generate one TR photon
    G4double kinEnergy = aParticle->GetKineticEnergy();
    G4double mass      = aParticle->GetDefinition()->GetPDGMass();
    G4double gamma     = 1.0 + kinEnergy / mass;

    if(verboseLevel > 1)
    {
      G4cout << "gamma = " << gamma << G4endl;
    }
    G4double massRatio           = proton_mass_c2 / mass;
    G4double TkinScaled          = kinEnergy * massRatio;
    G4ThreeVector position       = pPostStepPoint->GetPosition();
    G4ParticleMomentum direction = aParticle->GetMomentumDirection();
    G4double startTime           = pPostStepPoint->GetGlobalTime();

    for(iTkin = 0; iTkin < fTotBin; ++iTkin)
    {
      if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))
        break;
    }

    if(iTkin == 0)  // Tkin is too small, neglect of TR photon generation
    {
      if(verboseLevel > 0)
      {
        G4cout << "Go out from G4VXTRenergyLoss::PostStepDoIt:iTkin = " << iTkin
               << G4endl;
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    else  // general case: Tkin between two vectors of the material
    {
      fParticleChange.SetNumberOfSecondaries(1);

      energyTR = GetXTRrandomEnergy(TkinScaled, iTkin);

      if(verboseLevel > 1)
      {
        G4cout << "energyTR = " << energyTR / keV << " keV" << G4endl;
      }
      if(fAngleRadDistr)
      {
        theta2 = GetRandomAngle(energyTR, iTkin);
        if(theta2 > 0.)
          theta = std::sqrt(theta2);
        else
          theta = 0.;
      }
      else
        theta = std::fabs(G4RandGauss::shoot(0.0, pi / gamma));

      if(theta >= 0.1)
        theta = 0.1;

      phi = twopi * G4UniformRand();

      dirX = std::sin(theta) * std::cos(phi);
      dirY = std::sin(theta) * std::sin(phi);
      dirZ = std::cos(theta);

      G4ThreeVector directionTR(dirX, dirY, dirZ);
      directionTR.rotateUz(direction);
      directionTR.unit();

      auto aPhotonTR =
        new G4DynamicParticle(G4Gamma::Gamma(), directionTR, energyTR);

      // A XTR photon is set on the particle track inside the radiator
      // and is moved to the G4Envelope surface for standard X-ray TR models
      // only. The case of fExitFlux=true

      if(fExitFlux)
      {
        const G4RotationMatrix* rotM =
          pPostStepPoint->GetTouchable()->GetRotation();
        G4ThreeVector transl = pPostStepPoint->GetTouchable()->GetTranslation();
        G4AffineTransform transform = G4AffineTransform(rotM, transl);
        transform.Invert();
        G4ThreeVector localP = transform.TransformPoint(position);
        G4ThreeVector localV = transform.TransformAxis(directionTR);

        G4double distance =
          fEnvelope->GetSolid()->DistanceToOut(localP, localV);
        if(verboseLevel > 1)
        {
          G4cout << "distance to exit = " << distance / mm << " mm" << G4endl;
        }
        position += distance * directionTR;
        startTime += distance / c_light;
      }
      G4Track* aSecondaryTrack = new G4Track(aPhotonTR, startTime, position);
      aSecondaryTrack->SetTouchableHandle(
        aStep.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID(aTrack.GetTrackID());

      fParticleChange.AddSecondary(aSecondaryTrack);
      fParticleChange.ProposeEnergy(kinEnergy);
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

///////////////////////////////////////////////////////////////////////
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2, or 2->1)
// varAngle =2* (1 - std::cos(theta)) or approximately = theta*theta
G4complex G4VXTRenergyLoss::OneInterfaceXTRdEdx(G4double energy, G4double gamma,
                                                G4double varAngle)
{
  G4complex Z1 = GetPlateComplexFZ(energy, gamma, varAngle);
  G4complex Z2 = GetGasComplexFZ(energy, gamma, varAngle);

  G4complex zOut = (Z1 - Z2) * (Z1 - Z2) * (varAngle * energy / hbarc / hbarc);
  return zOut;
}

//////////////////////////////////////////////////////////////////////////////
// For photon energy distribution tables. Integrate first over angle
G4double G4VXTRenergyLoss::SpectralAngleXTRdEdx(G4double varAngle)
{
  G4double result = GetStackFactor(fEnergy, fGamma, varAngle);
  if(result < 0.0)
    result = 0.0;
  return result;
}

/////////////////////////////////////////////////////////////////////////
// For second integration over energy
G4double G4VXTRenergyLoss::SpectralXTRdEdx(G4double energy)
{
  G4int i;
  static constexpr G4int iMax = 8;
  G4double angleSum           = 0.0;

  G4double lim[iMax] = { 0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 };

  for(i = 0; i < iMax; ++i)
    lim[i] *= fMaxThetaTR;

  G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
    integral;

  fEnergy = energy;
  {
    for(i = 0; i < iMax - 1; ++i)
    {
      angleSum += integral.Legendre96(
        this, &G4VXTRenergyLoss::SpectralAngleXTRdEdx, lim[i], lim[i + 1]);
    }
  }
  return angleSum;
}

//////////////////////////////////////////////////////////////////////////
// for photon angle distribution tables
G4double G4VXTRenergyLoss::AngleSpectralXTRdEdx(G4double energy)
{
  G4double result = GetStackFactor(energy, fGamma, fVarAngle);
  if(result < 0)
    result = 0.0;
  return result;
}

///////////////////////////////////////////////////////////////////////////
// The XTR angular distribution based on transparent regular radiator
G4double G4VXTRenergyLoss::AngleXTRdEdx(G4double varAngle)
{
  G4double result;
  G4double sum = 0., tmp1, tmp2, tmp = 0., cof1, cof2, cofMin, cofPHC, energy1,
           energy2;
  G4int k, kMax, kMin, i;

  cofPHC = twopi * hbarc;

  cof1 = (fPlateThick + fGasThick) * (1. / fGamma / fGamma + varAngle);
  cof2 = fPlateThick * fSigma1 + fGasThick * fSigma2;

  cofMin = std::sqrt(cof1 * cof2);
  cofMin /= cofPHC;

  kMin = G4int(cofMin);
  if(cofMin > kMin)
    kMin++;

  kMax = kMin + 9;

  for(k = kMin; k <= kMax; ++k)
  {
    tmp1    = cofPHC * k;
    tmp2    = std::sqrt(tmp1 * tmp1 - cof1 * cof2);
    energy1 = (tmp1 + tmp2) / cof1;
    energy2 = (tmp1 - tmp2) / cof1;

    for(i = 0; i < 2; ++i)
    {
      if(i == 0)
      {
        if(energy1 > fTheMaxEnergyTR || energy1 < fTheMinEnergyTR)
          continue;

        tmp1 =
          (energy1 * energy1 * (1. / fGamma / fGamma + varAngle) + fSigma1) *
          fPlateThick / (4 * hbarc * energy1);
        tmp2 = std::sin(tmp1);
        tmp  = energy1 * tmp2 * tmp2;
        tmp2 = fPlateThick / (4. * tmp1);
        tmp1 =
          hbarc * energy1 /
          (energy1 * energy1 * (1. / fGamma / fGamma + varAngle) + fSigma2);
        tmp *= (tmp1 - tmp2) * (tmp1 - tmp2);
        tmp1 = cof1 / (4. * hbarc) - cof2 / (4. * hbarc * energy1 * energy1);
        tmp2 = std::abs(tmp1);

        if(tmp2 > 0.)
          tmp /= tmp2;
        else
          continue;
      }
      else
      {
        if(energy2 > fTheMaxEnergyTR || energy2 < fTheMinEnergyTR)
          continue;

        tmp1 =
          (energy2 * energy2 * (1. / fGamma / fGamma + varAngle) + fSigma1) *
          fPlateThick / (4. * hbarc * energy2);
        tmp2 = std::sin(tmp1);
        tmp  = energy2 * tmp2 * tmp2;
        tmp2 = fPlateThick / (4. * tmp1);
        tmp1 =
          hbarc * energy2 /
          (energy2 * energy2 * (1. / fGamma / fGamma + varAngle) + fSigma2);
        tmp *= (tmp1 - tmp2) * (tmp1 - tmp2);
        tmp1 = cof1 / (4. * hbarc) - cof2 / (4. * hbarc * energy2 * energy2);
        tmp2 = std::abs(tmp1);

        if(tmp2 > 0.)
          tmp /= tmp2;
        else
          continue;
      }
      sum += tmp;
    }
  }
  result = 4. * pi * fPlateNumber * sum * varAngle;
  result /= hbarc * hbarc;

  return result;
}

//////////////////////////////////////////////////////////////////////
// Calculates formation zone for plates. Omega is energy !!!
G4double G4VXTRenergyLoss::GetPlateFormationZone(G4double omega, G4double gamma,
                                                 G4double varAngle)
{
  G4double cof, lambda;
  lambda = 1.0 / gamma / gamma + varAngle + fSigma1 / omega / omega;
  cof    = 2.0 * hbarc / omega / lambda;
  return cof;
}

//////////////////////////////////////////////////////////////////////
// Calculates complex formation zone for plates. Omega is energy !!!
G4complex G4VXTRenergyLoss::GetPlateComplexFZ(G4double omega, G4double gamma,
                                              G4double varAngle)
{
  G4double cof, length, delta, real_v, image_v;

  length = 0.5 * GetPlateFormationZone(omega, gamma, varAngle);
  delta  = length * GetPlateLinearPhotoAbs(omega);
  cof    = 1.0 / (1.0 + delta * delta);

  real_v  = length * cof;
  image_v = real_v * delta;

  G4complex zone(real_v, image_v);
  return zone;
}

////////////////////////////////////////////////////////////////////////
// Computes matrix of Sandia photo absorption cross section coefficients for
// plate material
void G4VXTRenergyLoss::ComputePlatePhotoAbsCof()
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4Material* mat                   = (*theMaterialTable)[fMatIndex1];
  fPlatePhotoAbsCof                       = mat->GetSandiaTable();

  return;
}

//////////////////////////////////////////////////////////////////////
// Returns the value of linear photo absorption coefficient (in reciprocal
// length) for plate for given energy of X-ray photon omega
G4double G4VXTRenergyLoss::GetPlateLinearPhotoAbs(G4double omega)
{
  G4double omega2, omega3, omega4;

  omega2 = omega * omega;
  omega3 = omega2 * omega;
  omega4 = omega2 * omega2;

  const G4double* SandiaCof = fPlatePhotoAbsCof->GetSandiaCofForMaterial(omega);
  G4double cross            = SandiaCof[0] / omega + SandiaCof[1] / omega2 +
                   SandiaCof[2] / omega3 + SandiaCof[3] / omega4;
  return cross;
}

//////////////////////////////////////////////////////////////////////
// Calculates formation zone for gas. Omega is energy !!!
G4double G4VXTRenergyLoss::GetGasFormationZone(G4double omega, G4double gamma,
                                               G4double varAngle)
{
  G4double cof, lambda;
  lambda = 1.0 / gamma / gamma + varAngle + fSigma2 / omega / omega;
  cof    = 2.0 * hbarc / omega / lambda;
  return cof;
}

//////////////////////////////////////////////////////////////////////
// Calculates complex formation zone for gas gaps. Omega is energy !!!
G4complex G4VXTRenergyLoss::GetGasComplexFZ(G4double omega, G4double gamma,
                                            G4double varAngle)
{
  G4double cof, length, delta, real_v, image_v;

  length = 0.5 * GetGasFormationZone(omega, gamma, varAngle);
  delta  = length * GetGasLinearPhotoAbs(omega);
  cof    = 1.0 / (1.0 + delta * delta);

  real_v  = length * cof;
  image_v = real_v * delta;

  G4complex zone(real_v, image_v);
  return zone;
}

////////////////////////////////////////////////////////////////////////
// Computes matrix of Sandia photo absorption cross section coefficients for
// gas material
void G4VXTRenergyLoss::ComputeGasPhotoAbsCof()
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4Material* mat                   = (*theMaterialTable)[fMatIndex2];
  fGasPhotoAbsCof                         = mat->GetSandiaTable();
  return;
}

//////////////////////////////////////////////////////////////////////
// Returns the value of linear photo absorption coefficient (in reciprocal
// length) for gas
G4double G4VXTRenergyLoss::GetGasLinearPhotoAbs(G4double omega)
{
  G4double omega2, omega3, omega4;

  omega2 = omega * omega;
  omega3 = omega2 * omega;
  omega4 = omega2 * omega2;

  const G4double* SandiaCof = fGasPhotoAbsCof->GetSandiaCofForMaterial(omega);
  G4double cross            = SandiaCof[0] / omega + SandiaCof[1] / omega2 +
                   SandiaCof[2] / omega3 + SandiaCof[3] / omega4;
  return cross;
}

//////////////////////////////////////////////////////////////////////
// Calculates the product of linear cof by formation zone for plate.
// Omega is energy !!!
G4double G4VXTRenergyLoss::GetPlateZmuProduct(G4double omega, G4double gamma,
                                              G4double varAngle)
{
  return GetPlateFormationZone(omega, gamma, varAngle) *
         GetPlateLinearPhotoAbs(omega);
}
//////////////////////////////////////////////////////////////////////
// Calculates the product of linear cof by formation zone for plate.
// G4cout and output in file in some energy range.
void G4VXTRenergyLoss::GetPlateZmuProduct()
{
  std::ofstream outPlate("plateZmu.dat", std::ios::out);
  outPlate.setf(std::ios::scientific, std::ios::floatfield);

  G4int i;
  G4double omega, varAngle, gamma;
  gamma    = 10000.;
  varAngle = 1 / gamma / gamma;
  if(verboseLevel > 0)
    G4cout << "energy, keV" << "\t" << "Zmu for plate" << G4endl;
  for(i = 0; i < 100; ++i)
  {
    omega = (1.0 + i) * keV;
    if(verboseLevel > 1)
      G4cout << omega / keV << "\t"
             << GetPlateZmuProduct(omega, gamma, varAngle) << "\t";
    if(verboseLevel > 0)
      outPlate << omega / keV << "\t\t"
               << GetPlateZmuProduct(omega, gamma, varAngle) << G4endl;
  }
  return;
}

//////////////////////////////////////////////////////////////////////
// Calculates the product of linear cof by formation zone for gas.
// Omega is energy !!!
G4double G4VXTRenergyLoss::GetGasZmuProduct(G4double omega, G4double gamma,
                                            G4double varAngle)
{
  return GetGasFormationZone(omega, gamma, varAngle) *
         GetGasLinearPhotoAbs(omega);
}

//////////////////////////////////////////////////////////////////////
// Calculates the product of linear cof by formation zone for gas.
// G4cout and output in file in some energy range.
void G4VXTRenergyLoss::GetGasZmuProduct()
{
  std::ofstream outGas("gasZmu.dat", std::ios::out);
  outGas.setf(std::ios::scientific, std::ios::floatfield);
  G4int i;
  G4double omega, varAngle, gamma;
  gamma    = 10000.;
  varAngle = 1 / gamma / gamma;
  if(verboseLevel > 0)
    G4cout << "energy, keV" << "\t" << "Zmu for gas" << G4endl;
  for(i = 0; i < 100; ++i)
  {
    omega = (1.0 + i) * keV;
    if(verboseLevel > 1)
      G4cout << omega / keV << "\t" << GetGasZmuProduct(omega, gamma, varAngle)
             << "\t";
    if(verboseLevel > 0)
      outGas << omega / keV << "\t\t"
             << GetGasZmuProduct(omega, gamma, varAngle) << G4endl;
  }
  return;
}

////////////////////////////////////////////////////////////////////////
// Computes Compton cross section for plate material in 1/mm
G4double G4VXTRenergyLoss::GetPlateCompton(G4double omega)
{
  G4int i, numberOfElements;
  G4double xSection = 0., nowZ, sumZ = 0.;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  numberOfElements = (G4int)(*theMaterialTable)[fMatIndex1]->GetNumberOfElements();

  for(i = 0; i < numberOfElements; ++i)
  {
    nowZ = (*theMaterialTable)[fMatIndex1]->GetElement(i)->GetZ();
    sumZ += nowZ;
    xSection += GetComptonPerAtom(omega, nowZ);
  }
  xSection /= sumZ;
  xSection *= (*theMaterialTable)[fMatIndex1]->GetElectronDensity();
  return xSection;
}

////////////////////////////////////////////////////////////////////////
// Computes Compton cross section for gas material in 1/mm
G4double G4VXTRenergyLoss::GetGasCompton(G4double omega)
{
  G4int i, numberOfElements;
  G4double xSection = 0., nowZ, sumZ = 0.;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  numberOfElements = (G4int)(*theMaterialTable)[fMatIndex2]->GetNumberOfElements();

  for(i = 0; i < numberOfElements; ++i)
  {
    nowZ = (*theMaterialTable)[fMatIndex2]->GetElement(i)->GetZ();
    sumZ += nowZ;
    xSection += GetComptonPerAtom(omega, nowZ);
  }
  xSection /= sumZ;
  xSection *= (*theMaterialTable)[fMatIndex2]->GetElectronDensity();
  return xSection;
}

////////////////////////////////////////////////////////////////////////
// Computes Compton cross section per atom with Z electrons for gamma with
// the energy GammaEnergy
G4double G4VXTRenergyLoss::GetComptonPerAtom(G4double GammaEnergy, G4double Z)
{
  G4double CrossSection = 0.0;
  if(Z < 0.9999)
    return CrossSection;
  if(GammaEnergy < 0.1 * keV)
    return CrossSection;
  if(GammaEnergy > (100. * GeV / Z))
    return CrossSection;

  static constexpr G4double a = 20.0;
  static constexpr G4double b = 230.0;
  static constexpr G4double c = 440.0;

  static constexpr G4double d1 = 2.7965e-1 * barn, d2 = -1.8300e-1 * barn,
                            d3 = 6.7527 * barn, d4 = -1.9798e+1 * barn,
                            e1 = 1.9756e-5 * barn, e2 = -1.0205e-2 * barn,
                            e3 = -7.3913e-2 * barn, e4 = 2.7079e-2 * barn,
                            f1 = -3.9178e-7 * barn, f2 = 6.8241e-5 * barn,
                            f3 = 6.0480e-5 * barn, f4 = 3.0274e-4 * barn;

  G4double p1Z = Z * (d1 + e1 * Z + f1 * Z * Z);
  G4double p2Z = Z * (d2 + e2 * Z + f2 * Z * Z);
  G4double p3Z = Z * (d3 + e3 * Z + f3 * Z * Z);
  G4double p4Z = Z * (d4 + e4 * Z + f4 * Z * Z);

  G4double T0 = 15.0 * keV;
  if(Z < 1.5)
    T0 = 40.0 * keV;

  G4double X = std::max(GammaEnergy, T0) / electron_mass_c2;
  CrossSection =
    p1Z * std::log(1. + 2. * X) / X +
    (p2Z + p3Z * X + p4Z * X * X) / (1. + a * X + b * X * X + c * X * X * X);

  //  modification for low energy. (special case for Hydrogen)
  if(GammaEnergy < T0)
  {
    G4double dT0 = 1. * keV;
    X            = (T0 + dT0) / electron_mass_c2;
    G4double sigma =
      p1Z * std::log(1. + 2. * X) / X +
      (p2Z + p3Z * X + p4Z * X * X) / (1. + a * X + b * X * X + c * X * X * X);
    G4double c1 = -T0 * (sigma - CrossSection) / (CrossSection * dT0);
    G4double c2 = 0.150;
    if(Z > 1.5)
      c2 = 0.375 - 0.0556 * std::log(Z);
    G4double y = std::log(GammaEnergy / T0);
    CrossSection *= std::exp(-y * (c1 + c2 * y));
  }
  return CrossSection;
}

///////////////////////////////////////////////////////////////////////
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2, or 2->1)
// varAngle =2* (1 - std::cos(theta)) or approximately = theta*theta
G4double G4VXTRenergyLoss::OneBoundaryXTRNdensity(G4double energy,
                                                  G4double gamma,
                                                  G4double varAngle) const
{
  G4double formationLength1, formationLength2;
  formationLength1 =
    1.0 / (1.0 / (gamma * gamma) + fSigma1 / (energy * energy) + varAngle);
  formationLength2 =
    1.0 / (1.0 / (gamma * gamma) + fSigma2 / (energy * energy) + varAngle);
  return (varAngle / energy) * (formationLength1 - formationLength2) *
         (formationLength1 - formationLength2);
}

G4double G4VXTRenergyLoss::GetStackFactor(G4double energy, G4double gamma,
                                          G4double varAngle)
{
  // return stack factor corresponding to one interface
  return std::real(OneInterfaceXTRdEdx(energy, gamma, varAngle));
}

//////////////////////////////////////////////////////////////////////////////
// For photon energy distribution tables. Integrate first over angle
G4double G4VXTRenergyLoss::XTRNSpectralAngleDensity(G4double varAngle)
{
  return OneBoundaryXTRNdensity(fEnergy, fGamma, varAngle) *
         GetStackFactor(fEnergy, fGamma, varAngle);
}

/////////////////////////////////////////////////////////////////////////
// For second integration over energy
G4double G4VXTRenergyLoss::XTRNSpectralDensity(G4double energy)
{
  fEnergy = energy;
  G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
    integral;
  return integral.Legendre96(this, &G4VXTRenergyLoss::XTRNSpectralAngleDensity,
                             0.0, 0.2 * fMaxThetaTR) +
         integral.Legendre10(this, &G4VXTRenergyLoss::XTRNSpectralAngleDensity,
                             0.2 * fMaxThetaTR, fMaxThetaTR);
}

//////////////////////////////////////////////////////////////////////////
// for photon angle distribution tables
G4double G4VXTRenergyLoss::XTRNAngleSpectralDensity(G4double energy)
{
  return OneBoundaryXTRNdensity(energy, fGamma, fVarAngle) *
         GetStackFactor(energy, fGamma, fVarAngle);
}

///////////////////////////////////////////////////////////////////////////
G4double G4VXTRenergyLoss::XTRNAngleDensity(G4double varAngle)
{
  fVarAngle = varAngle;
  G4Integrator<G4VXTRenergyLoss, G4double (G4VXTRenergyLoss::*)(G4double)>
    integral;
  return integral.Legendre96(this, &G4VXTRenergyLoss::XTRNAngleSpectralDensity,
                             fMinEnergyTR, fMaxEnergyTR);
}

//////////////////////////////////////////////////////////////////////////////
// Check number of photons for a range of Lorentz factors from both energy
// and angular tables
void G4VXTRenergyLoss::GetNumberOfPhotons()
{
  G4int iTkin;
  G4double gamma, numberE;

  std::ofstream outEn("numberE.dat", std::ios::out);
  outEn.setf(std::ios::scientific, std::ios::floatfield);

  std::ofstream outAng("numberAng.dat", std::ios::out);
  outAng.setf(std::ios::scientific, std::ios::floatfield);

  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    gamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);
    numberE = (*(*fEnergyDistrTable)(iTkin))(0);
    if(verboseLevel > 1)
      G4cout << gamma << "\t\t" << numberE << "\t" << G4endl;
    if(verboseLevel > 0)
      outEn << gamma << "\t\t" << numberE << G4endl;
  }
  return;
}

/////////////////////////////////////////////////////////////////////////
// Returns random energy of a X-ray TR photon for given scaled kinetic energy
// of a charged particle
G4double G4VXTRenergyLoss::GetXTRrandomEnergy(G4double scaledTkin, G4int iTkin)
{
  G4int iTransfer, iPlace;
  G4double transfer = 0.0, position, E1, E2, W1, W2, W;

  iPlace = iTkin - 1;

  if(iTkin == fTotBin)  // relativistic plato, try from left
  {
    position = (*(*fEnergyDistrTable)(iPlace))(0) * G4UniformRand();

    for(iTransfer = 0;; ++iTransfer)
    {
      if(position >= (*(*fEnergyDistrTable)(iPlace))(iTransfer))
        break;
    }
    transfer = GetXTRenergy(iPlace, position, iTransfer);
  }
  else
  {
    E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1);
    E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
    W  = 1.0 / (E2 - E1);
    W1 = (E2 - scaledTkin) * W;
    W2 = (scaledTkin - E1) * W;

    position = ((*(*fEnergyDistrTable)(iPlace))(0) * W1 +
                (*(*fEnergyDistrTable)(iPlace + 1))(0) * W2) *
               G4UniformRand();

    for(iTransfer = 0;; ++iTransfer)
    {
      if(position >= ((*(*fEnergyDistrTable)(iPlace))(iTransfer) *W1 +
                      (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer) *W2))
        break;
    }
    transfer = GetXTRenergy(iPlace, position, iTransfer);
  }
  if(transfer < 0.0)
    transfer = 0.0;
  return transfer;
}

////////////////////////////////////////////////////////////////////////
// Returns approximate position of X-ray photon energy during random sampling
// over integral energy distribution
G4double G4VXTRenergyLoss::GetXTRenergy(G4int iPlace, G4double, G4int iTransfer)
{
  G4double x1, x2, y1, y2, result;

  if(iTransfer == 0)
  {
    result = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  }
  else
  {
    y1 = (*(*fEnergyDistrTable)(iPlace))(iTransfer - 1);
    y2 = (*(*fEnergyDistrTable)(iPlace))(iTransfer);

    x1 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer - 1);
    x2 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    if(x1 == x2)
      result = x2;
    else
    {
      if(y1 == y2)
        result = x1 + (x2 - x1) * G4UniformRand();
      else
      {
        result = x1 + (x2 - x1) * G4UniformRand();
      }
    }
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////
//  Get XTR photon angle at given energy and Tkin

G4double G4VXTRenergyLoss::GetRandomAngle(G4double energyXTR, G4int iTkin)
{
  G4int iTR, iAngle;
  G4double position, angle;

  if(iTkin == fTotBin)
    --iTkin;

  fAngleForEnergyTable = fAngleBank[iTkin];

  for(iTR = 0; iTR < fBinTR; ++iTR)
  {
    if(energyXTR < fXTREnergyVector->GetLowEdgeEnergy(iTR))
      break;
  }
  if(iTR == fBinTR)
    --iTR;

  position = (*(*fAngleForEnergyTable)(iTR))(0) * G4UniformRand();
  // position = (*(*fAngleForEnergyTable)(iTR))(1) * G4UniformRand(); // ATLAS TB

  for(iAngle = 0;; ++iAngle)
  // for(iAngle = 1;; ++iAngle) // ATLAS TB
  {
    if(position >= (*(*fAngleForEnergyTable)(iTR))(iAngle))
      break;
  }
  angle = GetAngleXTR(iTR, position, iAngle);
  return angle;
}

////////////////////////////////////////////////////////////////////////
// Returns approximate position of X-ray photon angle at given energy during
// random sampling over integral energy distribution

G4double G4VXTRenergyLoss::GetAngleXTR(G4int iPlace, G4double position,
                                       G4int iTransfer)
{
  G4double x1, x2, y1, y2, result;

  if( iTransfer == 0 )
  // if( iTransfer == 1 ) // ATLAS TB
  {
    result = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  }
  else
  {
    y1 = (*(*fAngleForEnergyTable)(iPlace))(iTransfer - 1);
    y2 = (*(*fAngleForEnergyTable)(iPlace))(iTransfer);

    x1 = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer - 1);
    x2 = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    if(x1 == x2) result = x2;
    else
    {
      if( y1 == y2 )  result = x1 + (x2 - x1) * G4UniformRand();
      else
      {
        result = x1 + (position - y1) * (x2 - x1) / (y2 - y1);
        // result = x1 + 0.1*(position - y1) * (x2 - x1) / (y2 - y1); // ATLAS TB
        // result = x1 + 0.05*(position - y1) * (x2 - x1) / (y2 - y1); // ATLAS TB
      }
    }
  }
  return result;
}
