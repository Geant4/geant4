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
// --------
// 09 Dec 2009   L Pandola    First implementation
//
#include "G4PenelopeSamplingData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
G4PenelopeSamplingData::G4PenelopeSamplingData(G4int nPoints) : 
  fNP(nPoints)
{
  //create vectors
  fX = new G4DataVector();
  fPAC = new G4DataVector();
  fA = new G4DataVector();
  fB = new G4DataVector();
  fITTL = new std::vector<size_t>;
  fITTU = new std::vector<size_t>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
G4PenelopeSamplingData::~G4PenelopeSamplingData()
{
  if (fX) delete fX;
  if (fPAC) delete fPAC;
  if (fA) delete fA;
  if (fB) delete fB;
  if (fITTL) delete fITTL;
  if (fITTU) delete fITTU;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.......oooOO0OOooo...

size_t G4PenelopeSamplingData::GetNumberOfStoredPoints()
{
  size_t points = fX->size();

  //check everything is all right
  if (fPAC->size() != points || fA->size() != points || 
      fB->size() != points || fITTL->size() != points ||
      fITTU->size() != points)
    {
      G4ExceptionDescription ed;
      ed << "Data vectors look to have different dimensions !" << G4endl;
      G4Exception("G4PenelopeSamplingData::GetNumberOfStoredPoints()","em2040",
		  FatalException,ed);      
    }
  return points;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
void G4PenelopeSamplingData::Clear()
{
  if (fX) delete fX;
  if (fPAC) delete fPAC;
  if (fA) delete fA;
  if (fB) delete fB;
  if (fITTL) delete fITTL;
  if (fITTU) delete fITTU;
  //create vectors
  fX = new G4DataVector();
  fPAC = new G4DataVector();
  fA = new G4DataVector();
  fB = new G4DataVector();
  fITTL = new std::vector<size_t>;
  fITTU = new std::vector<size_t>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
void G4PenelopeSamplingData::AddPoint(G4double x0,G4double pac0,G4double a0,G4double b0,
					size_t ITTL0,size_t ITTU0)
{
  fX->push_back(x0);
  fPAC->push_back(pac0);
  fA->push_back(a0);
  fB->push_back(b0);
  fITTL->push_back(ITTL0);
  fITTU->push_back(ITTU0);

  //check how many points we do have now
  size_t nOfPoints = GetNumberOfStoredPoints();

  if (nOfPoints > ((size_t)fNP))
    {
      G4cout << "G4PenelopeSamplingData::AddPoint() " << G4endl;
      G4cout << "WARNING: Up to now there are " << nOfPoints << " points in the table" << G4endl;
      G4cout << "while the anticipated (declared) number is " << fNP << G4endl;
    }
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
void G4PenelopeSamplingData::DumpTable()
{
  
  G4cout << "*************************************************************************" << G4endl;
  G4cout << GetNumberOfStoredPoints() << " points" << G4endl;
  G4cout << "*************************************************************************" << G4endl;
  for (size_t i=0;i<GetNumberOfStoredPoints();i++)
    {
      G4cout << i << " " << (*fX)[i] << " " << (*fPAC)[i] << " " << (*fA)[i] << " " << 
	(*fB)[i] << " " << (*fITTL)[i] << " " << (*fITTU)[i] << G4endl;
    }
  G4cout << "*************************************************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetX(size_t index)
{
  if (index < fX->size())
    return (*fX)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetPAC(size_t index)
{
  if (index < fPAC->size())
    return (*fPAC)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetA(size_t index)
{
  if (index < fA->size())
    return (*fA)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetB(size_t index)
{
  if (index < fB->size())
    return (*fB)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::SampleValue(G4double maxRand)
{
  //One passes here a random number in (0,1).
  //Notice: it possible that is between (0,b) with b<1
  size_t points = GetNumberOfStoredPoints();
 
  size_t itn = (size_t) (maxRand*(points-1)); 
  size_t i = (*fITTL)[itn];
  size_t j = (*fITTU)[itn];

  while ((j-i) > 1)
    {
      size_t k = (i+j)/2;
      if (maxRand > (*fPAC)[k])
	i = k;
      else
	j = k;
    }

  //Sampling from the rational inverse cumulative distribution
  G4double result = 0;

  G4double rr = maxRand - (*fPAC)[i];
  if (rr > 1e-16)
    {
      G4double d = (*fPAC)[i+1]-(*fPAC)[i];
      result = (*fX)[i]+
	((1.0+(*fA)[i]+(*fB)[i])*d*rr/
	 (d*d+((*fA)[i]*d+(*fB)[i]*rr)*rr))*((*fX)[i+1]-(*fX)[i]);      
    }
  else
    result = (*fX)[i]; 
  
  return result;
}
