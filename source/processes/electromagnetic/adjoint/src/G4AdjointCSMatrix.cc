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

#include "G4AdjointCSMatrix.hh"

#include "G4AdjointInterpolator.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>
#include <fstream>

///////////////////////////////////////////////////////
G4AdjointCSMatrix::G4AdjointCSMatrix(G4bool aBool) { fScatProjToProj = aBool; }

///////////////////////////////////////////////////////
G4AdjointCSMatrix::~G4AdjointCSMatrix()
{
  fLogPrimEnergyVector.clear();
  fLogCrossSectionVector.clear();

  for (auto p : fLogSecondEnergyMatrix) {
    p->clear();
    delete p;
    p = nullptr;
  }
  fLogSecondEnergyMatrix.clear();

  for (auto p : fLogProbMatrix) {
    p->clear();
    delete p;
    p = nullptr;
  }
  fLogProbMatrix.clear();

  for (auto p : fLogProbMatrixIndex) {
    if (p) {
      p->clear();
      delete p;
      p = nullptr;
    }
  }
  fLogProbMatrixIndex.clear();
}

///////////////////////////////////////////////////////
void G4AdjointCSMatrix::Clear()
{
  fLogPrimEnergyVector.clear();
  fLogCrossSectionVector.clear();
  fLogSecondEnergyMatrix.clear();
  fLogProbMatrix.clear();
  fLogProbMatrixIndex.clear();
  fLog0Vector.clear();
  fNbPrimEnergy = 0;
}

///////////////////////////////////////////////////////
void G4AdjointCSMatrix::AddData(G4double aLogPrimEnergy, G4double aLogCS,
                                std::vector<double>* aLogSecondEnergyVector,
                                std::vector<double>* aLogProbVector,
                                size_t n_pro_decade)
{
  G4AdjointInterpolator* theInterpolator = G4AdjointInterpolator::GetInstance();

  // At this time we consider that the energy is increasing monotically
  fLogPrimEnergyVector.push_back(aLogPrimEnergy);
  fLogCrossSectionVector.push_back(aLogCS);
  fLogSecondEnergyMatrix.push_back(aLogSecondEnergyVector);
  fLogProbMatrix.push_back(aLogProbVector);

  std::vector<size_t>* aLogProbVectorIndex = nullptr;

  if(n_pro_decade > 0 && !aLogProbVector->empty())
  {
    aLogProbVectorIndex = new std::vector<size_t>();
    G4double dlog       = std::log(10.) / n_pro_decade;
    G4double log_val =
      int(std::min((*aLogProbVector)[0], aLogProbVector->back()) / dlog) * dlog;
    fLog0Vector.push_back(log_val);

    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    while(log_val < 0.)
    {
      aLogProbVectorIndex->push_back(
        theInterpolator->FindPosition(log_val, (*aLogProbVector)));
      log_val += dlog;
    }
  }
  else
  {
    fLog0Vector.push_back(0.);
  }
  fLogProbMatrixIndex.push_back(aLogProbVectorIndex);

  ++fNbPrimEnergy;
}

///////////////////////////////////////////////////////
G4bool G4AdjointCSMatrix::GetData(unsigned int i, G4double& aLogPrimEnergy,
                                  G4double& aLogCS, G4double& log0,
                                  std::vector<double>*& aLogSecondEnergyVector,
                                  std::vector<double>*& aLogProbVector,
                                  std::vector<size_t>*& aLogProbVectorIndex)
{
  if(i >= fNbPrimEnergy)
    return false;
  aLogPrimEnergy         = fLogPrimEnergyVector[i];
  aLogCS                 = fLogCrossSectionVector[i];
  aLogSecondEnergyVector = fLogSecondEnergyMatrix[i];
  aLogProbVector         = fLogProbMatrix[i];
  aLogProbVectorIndex    = fLogProbMatrixIndex[i];
  log0                   = fLog0Vector[i];
  return true;
}

///////////////////////////////////////////////////////
void G4AdjointCSMatrix::Write(G4String file_name)
{
  std::fstream FileOutput(file_name, std::ios::out);
  FileOutput << std::setiosflags(std::ios::scientific);
  FileOutput << std::setprecision(6);
  FileOutput << fLogPrimEnergyVector.size() << G4endl;
  for(size_t i = 0; i < fLogPrimEnergyVector.size(); ++i)
  {
    FileOutput << std::exp(fLogPrimEnergyVector[i]) / MeV << '\t'
               << std::exp(fLogCrossSectionVector[i]) << G4endl;
    size_t j1 = 0;
    FileOutput << fLogSecondEnergyMatrix[i]->size() << G4endl;
    for(size_t j = 0; j < fLogSecondEnergyMatrix[i]->size(); ++j)
    {
      FileOutput << std::exp((*fLogSecondEnergyMatrix[i])[j]);
      ++j1;
      if(j1 < 10)
        FileOutput << '\t';
      else
      {
        FileOutput << G4endl;
        j1 = 0;
      }
    }
    if(j1 > 0)
      FileOutput << G4endl;
    j1 = 0;
    FileOutput << fLogProbMatrix[i]->size() << G4endl;
    for(size_t j = 0; j < fLogProbMatrix[i]->size(); ++j)
    {
      FileOutput << std::exp((*fLogProbMatrix[i])[j]);
      ++j1;
      if(j1 < 10)
        FileOutput << '\t';
      else
      {
        FileOutput << G4endl;
        j1 = 0;
      }
    }
    if(j1 > 0)
      FileOutput << G4endl;
  }
}

///////////////////////////////////////////////////////
void G4AdjointCSMatrix::Read(G4String file_name)
{
  std::fstream FileOutput(file_name, std::ios::in);
  size_t n1, n2;

  fLogPrimEnergyVector.clear();
  fLogCrossSectionVector.clear();
  fLogSecondEnergyMatrix.clear();
  fLogProbMatrix.clear();
  FileOutput >> n1;
  for(size_t i = 0; i < n1; ++i)
  {
    G4double E, CS;
    FileOutput >> E >> CS;
    fLogPrimEnergyVector.push_back(E);
    fLogCrossSectionVector.push_back(CS);
    FileOutput >> n2;
    fLogSecondEnergyMatrix.push_back(new std::vector<G4double>());
    fLogProbMatrix.push_back(new std::vector<G4double>());

    for(size_t j = 0; j < n2; ++j)
    {
      G4double E1;
      FileOutput >> E1;
      fLogSecondEnergyMatrix[i]->push_back(E1);
    }
    FileOutput >> n2;
    for(size_t j = 0; j < n2; ++j)
    {
      G4double prob;
      FileOutput >> prob;
      fLogProbMatrix[i]->push_back(prob);
    }
  }
}
