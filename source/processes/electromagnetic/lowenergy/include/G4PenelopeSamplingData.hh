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
//
// Author: Luciano Pandola
//
// History:
// -----------
// 09 Dec 2009   L. Pandola   1st implementation. 
//
// -------------------------------------------------------------------
//
// Class description:
// This is a container of data that are used for sampling algorithm
// of Penelope08 Rayleigh scattering
// -------------------------------------------------------------------

#ifndef G4PENELOPESAMPLINGDATA_HH
#define G4PENELOPESAMPLINGDATA_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
//
//This is a container of data that are used for sampling algoritm
//
class G4PenelopeSamplingData
{
public:
  explicit G4PenelopeSamplingData(G4int npoints=150);
  ~G4PenelopeSamplingData();

  void AddPoint(G4double x0,G4double pac0,G4double a0,G4double b0,size_t ITTL0,
		size_t ITTU0);
  size_t GetNumberOfStoredPoints();
  void Clear();
  void DumpTable();
  
  G4double GetX(size_t index);
  G4double GetPAC(size_t index);
  G4double GetA(size_t index);
  G4double GetB(size_t index);

  G4double SampleValue(G4double rndm);

  G4PenelopeSamplingData & operator=(const G4PenelopeSamplingData &right) = delete;
  G4PenelopeSamplingData(const G4PenelopeSamplingData&) = delete;
  
private:  
  G4DataVector* fX; //grid points, in increasing order
  G4DataVector* fPAC; //value of the cumulative pdf at x_i
  G4DataVector* fA; // rational inverse cumulative inverse distribution parameters
  G4DataVector* fB;
  
  std::vector<size_t> *fITTL; //largest j for which pac(j) < (i-1)/(np-1)
  std::vector<size_t> *fITTU; //smallest k for which pac(k) > i/(np-1)

  G4int fNP; //number of grid points
};

#endif

