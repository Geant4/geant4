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

#include "G4VEmAdjointModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4AdjointInterpolator.hh"
#include "G4AdjointPositron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Integrator.hh"
#include "G4ParticleChange.hh"
#include "G4ProductionCutsTable.hh"
#include "G4TrackStatus.hh"

////////////////////////////////////////////////////////////////////////////////
G4VEmAdjointModel::G4VEmAdjointModel(const G4String& nam)
  : fName(nam)
{
  fCSManager = G4AdjointCSManager::GetAdjointCSManager();
  fCSManager->RegisterEmAdjointModel(this);
}

///////////////////////////////////e////////////////////////////////////////////
G4VEmAdjointModel::~G4VEmAdjointModel()
{
  delete fCSMatrixProdToProjBackScat;
  fCSMatrixProdToProjBackScat = nullptr;

  delete fCSMatrixProjToProjBackScat;
  fCSMatrixProjToProjBackScat = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::AdjointCrossSection(
  const G4MaterialCutsCouple* aCouple, G4double primEnergy,
  G4bool isScatProjToProj)
{
  DefineCurrentMaterial(aCouple);
  fPreStepEnergy = primEnergy;

  std::vector<G4double>* CS_Vs_Element = &fElementCSProdToProj;
  if(isScatProjToProj)
    CS_Vs_Element = &fElementCSScatProjToProj;
  fLastCS =
    fCSManager->ComputeAdjointCS(fCurrentMaterial, this, primEnergy,
                                 fTcutSecond, isScatProjToProj, *CS_Vs_Element);
  if(isScatProjToProj)
    fLastAdjointCSForScatProjToProj = fLastCS;
  else
    fLastAdjointCSForProdToProj = fLastCS;

  return fLastCS;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionPerAtomPrimToSecond(
  G4double kinEnergyProj, G4double kinEnergyProd, G4double Z, G4double A)
{
  G4double dSigmadEprod = 0.;
  G4double Emax_proj    = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double Emin_proj    = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);

  // the produced particle should have a kinetic energy less than the projectile
  if(kinEnergyProj > Emin_proj && kinEnergyProj <= Emax_proj)
  {
    G4double E1     = kinEnergyProd;
    G4double E2     = kinEnergyProd * 1.000001;
    G4double sigma1 = fDirectModel->ComputeCrossSectionPerAtom(
      fDirectPrimaryPart, kinEnergyProj, Z, A, E1, 1.e20);
    G4double sigma2 = fDirectModel->ComputeCrossSectionPerAtom(
      fDirectPrimaryPart, kinEnergyProj, Z, A, E2, 1.e20);

    dSigmadEprod = (sigma1 - sigma2) / (E2 - E1);
  }
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionPerAtomPrimToScatPrim(
  G4double kinEnergyProj, G4double kinEnergyScatProj, G4double Z, G4double A)
{
  G4double kinEnergyProd = kinEnergyProj - kinEnergyScatProj;
  G4double dSigmadEprod;
  if(kinEnergyProd <= 0.)
    dSigmadEprod = 0.;
  else
    dSigmadEprod =
      DiffCrossSectionPerAtomPrimToSecond(kinEnergyProj, kinEnergyProd, Z, A);
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToSecond(
  const G4Material* aMaterial, G4double kinEnergyProj, G4double kinEnergyProd)
{
  G4double dSigmadEprod = 0.;
  G4double Emax_proj    = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double Emin_proj    = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);

  if(kinEnergyProj > Emin_proj && kinEnergyProj <= Emax_proj)
  {
    G4double E1     = kinEnergyProd;
    G4double E2     = kinEnergyProd * 1.0001;
    G4double sigma1 = fDirectModel->CrossSectionPerVolume(
      aMaterial, fDirectPrimaryPart, kinEnergyProj, E1, 1.e20);
    G4double sigma2 = fDirectModel->CrossSectionPerVolume(
      aMaterial, fDirectPrimaryPart, kinEnergyProj, E2, 1.e20);
    dSigmadEprod = (sigma1 - sigma2) / (E2 - E1);
  }
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToScatPrim(
  const G4Material* aMaterial, G4double kinEnergyProj,
  G4double kinEnergyScatProj)
{
  G4double kinEnergyProd = kinEnergyProj - kinEnergyScatProj;
  G4double dSigmadEprod;
  if(kinEnergyProd <= 0.)
    dSigmadEprod = 0.;
  else
    dSigmadEprod = DiffCrossSectionPerVolumePrimToSecond(
      aMaterial, kinEnergyProj, kinEnergyProd);
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionFunction1(G4double kinEnergyProj)
{
  G4double bias_factor =
    fCsBiasingFactor * fKinEnergyProdForIntegration / kinEnergyProj;

  if(fUseMatrixPerElement)
  {
    return DiffCrossSectionPerAtomPrimToSecond(
             kinEnergyProj, fKinEnergyProdForIntegration, fZSelectedNucleus,
             fASelectedNucleus) *
           bias_factor;
  }
  else
  {
    return DiffCrossSectionPerVolumePrimToSecond(
             fSelectedMaterial, kinEnergyProj, fKinEnergyProdForIntegration) *
           bias_factor;
  }
}

////////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::DiffCrossSectionFunction2(G4double kinEnergyProj)
{
  G4double bias_factor =
    fCsBiasingFactor * fKinEnergyScatProjForIntegration / kinEnergyProj;
  if(fUseMatrixPerElement)
  {
    return DiffCrossSectionPerAtomPrimToScatPrim(
             kinEnergyProj, fKinEnergyScatProjForIntegration, fZSelectedNucleus,
             fASelectedNucleus) *
           bias_factor;
  }
  else
  {
    return DiffCrossSectionPerVolumePrimToScatPrim(
             fSelectedMaterial, kinEnergyProj,
             fKinEnergyScatProjForIntegration) *
           bias_factor;
  }
}

////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<G4double>*>
G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerAtomForSecond(
  G4double kinEnergyProd, G4double Z, G4double A,
  G4int nbin_pro_decade)  // nb bins per order of magnitude of energy
{
  G4Integrator<G4VEmAdjointModel, G4double (G4VEmAdjointModel::*)(G4double)>
    integral;
  fASelectedNucleus            = G4lrint(A);
  fZSelectedNucleus            = G4lrint(Z);
  fKinEnergyProdForIntegration = kinEnergyProd;

  // compute the vector of integrated cross sections
  G4double minEProj = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);
  G4double maxEProj = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double E1       = minEProj;
  std::vector<G4double>* log_ESec_vector = new std::vector<G4double>();
  std::vector<G4double>* log_Prob_vector = new std::vector<G4double>();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);

  G4double E2 =
    std::pow(10., G4double(G4int(std::log10(minEProj) * nbin_pro_decade) + 1) /
                    nbin_pro_decade);
  G4double fE = std::pow(10., 1. / nbin_pro_decade);

  if(std::pow(fE, 5.) > (maxEProj / minEProj))
    fE = std::pow(maxEProj / minEProj, 0.2);

  G4double int_cross_section = 0.;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(E1 < maxEProj * 0.9999999)
  {
    int_cross_section +=
      integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction1, E1,
                       std::min(E2, maxEProj * 0.99999999), 5);
    log_ESec_vector->push_back(std::log(std::min(E2, maxEProj)));
    log_Prob_vector->push_back(std::log(int_cross_section));
    E1 = E2;
    E2 *= fE;
  }
  std::vector<std::vector<G4double>*> res_mat;
  if(int_cross_section > 0.)
  {
    res_mat.push_back(log_ESec_vector);
    res_mat.push_back(log_Prob_vector);
  }
  else {
	delete  log_ESec_vector;
	delete  log_Prob_vector;
  }
  return res_mat;
}

/////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<G4double>*>
G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerAtomForScatProj(
  G4double kinEnergyScatProj, G4double Z, G4double A,
  G4int nbin_pro_decade)  // nb bins pro order of magnitude of energy
{
  G4Integrator<G4VEmAdjointModel, G4double (G4VEmAdjointModel::*)(G4double)>
    integral;
  fASelectedNucleus                = G4lrint(A);
  fZSelectedNucleus                = G4lrint(Z);
  fKinEnergyScatProjForIntegration = kinEnergyScatProj;

  // compute the vector of integrated cross sections
  G4double minEProj = GetSecondAdjEnergyMinForScatProjToProj(kinEnergyScatProj);
  G4double maxEProj = GetSecondAdjEnergyMaxForScatProjToProj(kinEnergyScatProj);
  G4double dEmax    = maxEProj - kinEnergyScatProj;
  G4double dEmin    = GetLowEnergyLimit();
  G4double dE1      = dEmin;
  G4double dE2      = dEmin;

  std::vector<G4double>* log_ESec_vector = new std::vector<G4double>();
  std::vector<G4double>* log_Prob_vector = new std::vector<G4double>();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins = std::max(G4int(std::log10(dEmax / dEmin)) * nbin_pro_decade, 5);
  G4double fE = std::pow(dEmax / dEmin, 1. / nbins);

  G4double int_cross_section = 0.;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(dE1 < dEmax * 0.9999999999999)
  {
    dE2 = dE1 * fE;
    int_cross_section +=
      integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction2,
                       minEProj + dE1, std::min(minEProj + dE2, maxEProj), 5);
    log_ESec_vector->push_back(std::log(std::min(dE2, maxEProj - minEProj)));
    log_Prob_vector->push_back(std::log(int_cross_section));
    dE1 = dE2;
  }

  std::vector<std::vector<G4double>*> res_mat;
  if(int_cross_section > 0.)
  {
    res_mat.push_back(log_ESec_vector);
    res_mat.push_back(log_Prob_vector);
  }
  else {
  	delete  log_ESec_vector;
  	delete  log_Prob_vector;
  }

  return res_mat;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<G4double>*>
G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForSecond(
  G4Material* aMaterial, G4double kinEnergyProd,
  G4int nbin_pro_decade)  // nb bins pro order of magnitude of energy
{
  G4Integrator<G4VEmAdjointModel, G4double (G4VEmAdjointModel::*)(G4double)>
    integral;
  fSelectedMaterial            = aMaterial;
  fKinEnergyProdForIntegration = kinEnergyProd;

  // compute the vector of integrated cross sections
  G4double minEProj = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);
  G4double maxEProj = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double E1       = minEProj;
  std::vector<G4double>* log_ESec_vector = new std::vector<G4double>();
  std::vector<G4double>* log_Prob_vector = new std::vector<G4double>();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);

  G4double E2 =
    std::pow(10., G4double(G4int(std::log10(minEProj) * nbin_pro_decade) + 1) /
                    nbin_pro_decade);
  G4double fE = std::pow(10., 1. / nbin_pro_decade);

  if(std::pow(fE, 5.) > (maxEProj / minEProj))
    fE = std::pow(maxEProj / minEProj, 0.2);

  G4double int_cross_section = 0.;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(E1 < maxEProj * 0.9999999)
  {
    int_cross_section +=
      integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction1, E1,
                       std::min(E2, maxEProj * 0.99999999), 5);
    log_ESec_vector->push_back(std::log(std::min(E2, maxEProj)));
    log_Prob_vector->push_back(std::log(int_cross_section));
    E1 = E2;
    E2 *= fE;
  }
  std::vector<std::vector<G4double>*> res_mat;

  if(int_cross_section > 0.)
  {
    res_mat.push_back(log_ESec_vector);
    res_mat.push_back(log_Prob_vector);
  }
  else {
  	delete  log_ESec_vector;
  	delete  log_Prob_vector;
  }
  return res_mat;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<G4double>*>
G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForScatProj(
  G4Material* aMaterial, G4double kinEnergyScatProj,
  G4int nbin_pro_decade)  // nb bins pro order of magnitude of energy
{
  G4Integrator<G4VEmAdjointModel, G4double (G4VEmAdjointModel::*)(G4double)>
    integral;
  fSelectedMaterial                = aMaterial;
  fKinEnergyScatProjForIntegration = kinEnergyScatProj;

  // compute the vector of integrated cross sections
  G4double minEProj = GetSecondAdjEnergyMinForScatProjToProj(kinEnergyScatProj);
  G4double maxEProj = GetSecondAdjEnergyMaxForScatProjToProj(kinEnergyScatProj);

  G4double dEmax = maxEProj - kinEnergyScatProj;
  G4double dEmin = GetLowEnergyLimit();
  G4double dE1   = dEmin;
  G4double dE2   = dEmin;

  std::vector<G4double>* log_ESec_vector = new std::vector<G4double>();
  std::vector<G4double>* log_Prob_vector = new std::vector<G4double>();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins = std::max(int(std::log10(dEmax / dEmin)) * nbin_pro_decade, 5);
  G4double fE = std::pow(dEmax / dEmin, 1. / nbins);

  G4double int_cross_section = 0.;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(dE1 < dEmax * 0.9999999999999)
  {
    dE2 = dE1 * fE;
    int_cross_section +=
      integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction2,
                       minEProj + dE1, std::min(minEProj + dE2, maxEProj), 5);
    log_ESec_vector->push_back(std::log(std::min(dE2, maxEProj - minEProj)));
    log_Prob_vector->push_back(std::log(int_cross_section));
    dE1 = dE2;
  }

  std::vector<std::vector<G4double>*> res_mat;
  if(int_cross_section > 0.)
  {
    res_mat.push_back(log_ESec_vector);
    res_mat.push_back(log_Prob_vector);
  }
  else {
  	delete  log_ESec_vector;
  	delete  log_Prob_vector;
  }

  return res_mat;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromCSMatrix(
  std::size_t MatrixIndex, G4double aPrimEnergy, G4bool isScatProjToProj)
{
  G4AdjointCSMatrix* theMatrix = (*fCSMatrixProdToProjBackScat)[MatrixIndex];
  if(isScatProjToProj)
    theMatrix = (*fCSMatrixProjToProjBackScat)[MatrixIndex];
  std::vector<G4double>* theLogPrimEnergyVector =
    theMatrix->GetLogPrimEnergyVector();

  if(theLogPrimEnergyVector->empty())
  {
    G4cout << "No data are contained in the given AdjointCSMatrix!" << G4endl;
    G4cout << "The sampling procedure will be stopped." << G4endl;
    return 0.;
  }

  G4AdjointInterpolator* theInterpolator = G4AdjointInterpolator::GetInstance();
  G4double aLogPrimEnergy                = std::log(aPrimEnergy);
  G4int ind = (G4int)theInterpolator->FindPositionForLogVector(
    aLogPrimEnergy, *theLogPrimEnergyVector);

  G4double aLogPrimEnergy1, aLogPrimEnergy2;
  G4double aLogCS1, aLogCS2;
  G4double log01, log02;
  std::vector<G4double>* aLogSecondEnergyVector1 = nullptr;
  std::vector<G4double>* aLogSecondEnergyVector2 = nullptr;
  std::vector<G4double>* aLogProbVector1         = nullptr;
  std::vector<G4double>* aLogProbVector2         = nullptr;
  std::vector<std::size_t>* aLogProbVectorIndex1    = nullptr;
  std::vector<std::size_t>* aLogProbVectorIndex2    = nullptr;

  theMatrix->GetData(ind, aLogPrimEnergy1, aLogCS1, log01,
                     aLogSecondEnergyVector1, aLogProbVector1,
					 aLogProbVectorIndex1 );
  theMatrix->GetData(ind + 1, aLogPrimEnergy2, aLogCS2, log02,
                     aLogSecondEnergyVector2, aLogProbVector2,
                     aLogProbVectorIndex2);

  if (! (aLogProbVector1 && aLogProbVector2 &&
		       aLogSecondEnergyVector1 && aLogSecondEnergyVector2)){
	 return  0.;
  }

  G4double rand_var     = G4UniformRand();
  G4double log_rand_var = std::log(rand_var);
  G4double log_Tcut     = std::log(fTcutSecond);
  G4double Esec         = 0.;
  G4double log_dE1, log_dE2;
  G4double log_rand_var1, log_rand_var2;
  G4double log_E1, log_E2;
  log_rand_var1 = log_rand_var;
  log_rand_var2 = log_rand_var;

  G4double Emin = 0.;
  G4double Emax = 0.;
  if(theMatrix->IsScatProjToProj())
  {  // case where Tcut plays a role
    Emin = GetSecondAdjEnergyMinForScatProjToProj(aPrimEnergy, fTcutSecond);
    Emax = GetSecondAdjEnergyMaxForScatProjToProj(aPrimEnergy);
    G4double dE = 0.;
    if(Emin < Emax)
    {
      if(fApplyCutInRange)
      {
        if(fSecondPartSameType && fTcutSecond > aPrimEnergy)
          return aPrimEnergy;

        log_rand_var1 = log_rand_var +
                        theInterpolator->InterpolateForLogVector(
                          log_Tcut, *aLogSecondEnergyVector1, *aLogProbVector1);
        log_rand_var2 = log_rand_var +
                        theInterpolator->InterpolateForLogVector(
                          log_Tcut, *aLogSecondEnergyVector2, *aLogProbVector2);
      }
      log_dE1 = theInterpolator->Interpolate(log_rand_var1, *aLogProbVector1,
                                             *aLogSecondEnergyVector1, "Lin");
      log_dE2 = theInterpolator->Interpolate(log_rand_var2, *aLogProbVector2,
                                             *aLogSecondEnergyVector2, "Lin");
      dE      = std::exp(theInterpolator->LinearInterpolation(
        aLogPrimEnergy, aLogPrimEnergy1, aLogPrimEnergy2, log_dE1, log_dE2));
    }

    Esec = aPrimEnergy + dE;
    Esec = std::max(Esec, Emin);
    Esec = std::min(Esec, Emax);
  }
  else
  {  // Tcut condition is already full-filled

    log_E1 = theInterpolator->Interpolate(log_rand_var, *aLogProbVector1,
                                          *aLogSecondEnergyVector1, "Lin");
    log_E2 = theInterpolator->Interpolate(log_rand_var, *aLogProbVector2,
                                          *aLogSecondEnergyVector2, "Lin");

    Esec = std::exp(theInterpolator->LinearInterpolation(
      aLogPrimEnergy, aLogPrimEnergy1, aLogPrimEnergy2, log_E1, log_E2));
    Emin = GetSecondAdjEnergyMinForProdToProj(aPrimEnergy);
    Emax = GetSecondAdjEnergyMaxForProdToProj(aPrimEnergy);
    Esec = std::max(Esec, Emin);
    Esec = std::min(Esec, Emax);
  }
  return Esec;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromCSMatrix(
  G4double aPrimEnergy, G4bool isScatProjToProj)
{
  SelectCSMatrix(isScatProjToProj);
  return SampleAdjSecEnergyFromCSMatrix(fCSMatrixUsed, aPrimEnergy,
                                        isScatProjToProj);
}

//////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::SelectCSMatrix(G4bool isScatProjToProj)
{
  fCSMatrixUsed = 0;
  if(!fUseMatrixPerElement)
    fCSMatrixUsed = fCurrentMaterial->GetIndex();
  else if(!fOneMatrixForAllElements)
  {  // Select Material
    std::vector<G4double>* CS_Vs_Element = &fElementCSScatProjToProj;
    fLastCS                              = fLastAdjointCSForScatProjToProj;
    if(!isScatProjToProj)
    {
      CS_Vs_Element = &fElementCSProdToProj;
      fLastCS       = fLastAdjointCSForProdToProj;
    }
    G4double SumCS = 0.;
    std::size_t ind = 0;
    for(std::size_t i = 0; i < CS_Vs_Element->size(); ++i)
    {
      SumCS += (*CS_Vs_Element)[i];
      if(G4UniformRand() <= SumCS / fLastCS)
      {
        ind = i;
        break;
      }
    }
    fCSMatrixUsed = fCurrentMaterial->GetElement((G4int)ind)->GetIndex();
  }
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromDiffCrossSectionPerAtom(
  G4double prim_energy, G4bool isScatProjToProj)
{
  // here we try to use the rejection method
  constexpr G4int iimax = 1000;
  G4double E            = 0.;
  G4double x, xmin, greject;
  if(isScatProjToProj)
  {
    G4double Emax = GetSecondAdjEnergyMaxForScatProjToProj(prim_energy);
    G4double Emin = prim_energy + fTcutSecond;
    xmin          = Emin / Emax;
    G4double grejmax =
      DiffCrossSectionPerAtomPrimToScatPrim(Emin, prim_energy, 1) * prim_energy;

    G4int ii = 0;
    do
    {
      // q = G4UniformRand();
      x = 1. / (G4UniformRand() * (1. / xmin - 1.) + 1.);
      E = x * Emax;
      greject =
        DiffCrossSectionPerAtomPrimToScatPrim(E, prim_energy, 1) * prim_energy;
      ++ii;
      if(ii >= iimax)
      {
        break;
      }
    }
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    while(greject < G4UniformRand() * grejmax);
  }
  else
  {
    G4double Emax = GetSecondAdjEnergyMaxForProdToProj(prim_energy);
    G4double Emin = GetSecondAdjEnergyMinForProdToProj(prim_energy);
    xmin          = Emin / Emax;
    G4double grejmax =
      DiffCrossSectionPerAtomPrimToSecond(Emin, prim_energy, 1);
    G4int ii = 0;
    do
    {
      x       = std::pow(xmin, G4UniformRand());
      E       = x * Emax;
      greject = DiffCrossSectionPerAtomPrimToSecond(E, prim_energy, 1);
      ++ii;
      if(ii >= iimax)
      {
        break;
      }
    }
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    while(greject < G4UniformRand() * grejmax);
  }

  return E;
}

////////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::CorrectPostStepWeight(G4ParticleChange* fParticleChange,
                                              G4double old_weight,
                                              G4double adjointPrimKinEnergy,
                                              G4double projectileKinEnergy,
                                              G4bool isScatProjToProj)
{
  G4double new_weight = old_weight;
  G4double w_corr =
    fCSManager->GetPostStepWeightCorrection() / fCsBiasingFactor;

  fLastCS = fLastAdjointCSForScatProjToProj;
  if(!isScatProjToProj)
    fLastCS = fLastAdjointCSForProdToProj;
  if((adjointPrimKinEnergy - fPreStepEnergy) / fPreStepEnergy > 0.001)
  {
    G4double post_stepCS = AdjointCrossSection(
      fCurrentCouple, adjointPrimKinEnergy, isScatProjToProj);
    if(post_stepCS > 0. && fLastCS > 0.)
      w_corr *= post_stepCS / fLastCS;
  }

  new_weight *= w_corr;
  new_weight *= projectileKinEnergy / adjointPrimKinEnergy;
  // This is needed due to the biasing of diff CS
  // by the factor adjointPrimKinEnergy/projectileKinEnergy

  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::GetSecondAdjEnergyMaxForScatProjToProj(
  G4double kinEnergyScatProj)
{
  G4double maxEProj = GetHighEnergyLimit();
  if(fSecondPartSameType)
    maxEProj = std::min(kinEnergyScatProj * 2., maxEProj);
  return maxEProj;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::GetSecondAdjEnergyMinForScatProjToProj(
  G4double primAdjEnergy, G4double tcut)
{
  G4double Emin = primAdjEnergy;
  if(fApplyCutInRange)
    Emin += tcut;
  return Emin;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::GetSecondAdjEnergyMaxForProdToProj(G4double)
{
  return fHighEnergyLimit;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VEmAdjointModel::GetSecondAdjEnergyMinForProdToProj(
  G4double primAdjEnergy)
{
  G4double minEProj = primAdjEnergy;
  if(fSecondPartSameType)
    minEProj = primAdjEnergy * 2.;
  return minEProj;
}

////////////////////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::DefineCurrentMaterial(
  const G4MaterialCutsCouple* couple)
{
  if(couple != fCurrentCouple)
  {
    fCurrentCouple   = const_cast<G4MaterialCutsCouple*>(couple);
    fCurrentMaterial = const_cast<G4Material*>(couple->GetMaterial());
    std::size_t idx       = 56;
    fTcutSecond      = 1.e-11;
    if(fAdjEquivDirectSecondPart)
    {
      if(fAdjEquivDirectSecondPart == G4AdjointGamma::AdjointGamma())
        idx = 0;
      else if(fAdjEquivDirectSecondPart == G4AdjointElectron::AdjointElectron())
        idx = 1;
      else if(fAdjEquivDirectSecondPart == G4AdjointPositron::AdjointPositron())
        idx = 2;
      if(idx < 56)
      {
        const std::vector<G4double>* aVec =
          G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(
            idx);
        fTcutSecond = (*aVec)[couple->GetIndex()];
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::SetHighEnergyLimit(G4double aVal)
{
  fHighEnergyLimit = aVal;
  if(fDirectModel)
    fDirectModel->SetHighEnergyLimit(aVal);
}

////////////////////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::SetLowEnergyLimit(G4double aVal)
{
  fLowEnergyLimit = aVal;
  if(fDirectModel)
    fDirectModel->SetLowEnergyLimit(aVal);
}

////////////////////////////////////////////////////////////////////////////////////////////
void G4VEmAdjointModel::SetAdjointEquivalentOfDirectPrimaryParticleDefinition(
  G4ParticleDefinition* aPart)
{
  fAdjEquivDirectPrimPart = aPart;
  if(fAdjEquivDirectPrimPart->GetParticleName() == "adj_e-")
    fDirectPrimaryPart = G4Electron::Electron();
  else if(fAdjEquivDirectPrimPart->GetParticleName() == "adj_gamma")
    fDirectPrimaryPart = G4Gamma::Gamma();
}
