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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointCSMatrix
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  An adjoint CS matrix is used by the model of a reverse process to sample
//  an adjoint secondary (being equivalent to a forward primary). It represents
//  the integration over the energy of the adjoint secondary (therefore the
//  forward primary) of the differential cross section of the equivalent forward
//  discrete process (Ionisation, Brem, PE effect, Compton,..). Each reverse
//  model has its own cross section matrix for a given cut, material couple. It
//  is therefore recomputed after a modification of the cuts by the user.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointCSMatrix_h
#define G4AdjointCSMatrix_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include <vector>

class G4AdjointCSMatrix
{
 public:
  G4AdjointCSMatrix(G4bool aBool);
  ~G4AdjointCSMatrix();

  void Clear();

  void AddData(G4double aPrimEnergy, G4double aCS,
               std::vector<double>* aLogSecondEnergyVector,
               std::vector<double>* aLogProbVector, size_t n_pro_decade = 0);

  G4bool GetData(unsigned int i, G4double& aPrimEnergy, G4double& aCS,
                 G4double& log0, std::vector<double>*& aLogSecondEnergyVector,
                 std::vector<double>*& aLogProbVector,
                 std::vector<size_t>*& aLogProbVectorIndex);

  inline std::vector<double>* GetLogPrimEnergyVector()
  {
    return &fLogPrimEnergyVector;
  }

  inline std::vector<double>* GetLogCrossSectionvector()
  {
    return &fLogCrossSectionVector;
  }

  inline G4bool IsScatProjToProj() { return fScatProjToProj; }

  void Write(G4String file_name);

  void Read(G4String file_name);

 private:
  std::vector<double> fLogPrimEnergyVector;
  // Adjoint Cross sections as functions of primary energy
  std::vector<double> fLogCrossSectionVector;

  std::vector<std::vector<double>*> fLogSecondEnergyMatrix;
  std::vector<std::vector<double>*> fLogProbMatrix;
  // Each column represents the integrated probability of
  // getting a secondary

  // index of equidistant LogProb
  std::vector<std::vector<size_t>*> fLogProbMatrixIndex;
  std::vector<double> fLog0Vector;

  size_t fNbPrimEnergy = 0;

  G4bool fScatProjToProj;
};
#endif
