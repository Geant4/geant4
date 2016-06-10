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
// $Id: G4PenelopeSamplingData.cc 66241 2012-12-13 18:34:42Z gunter $
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
  np(nPoints)
{
  //create vectors
  x = new G4DataVector();
  pac = new G4DataVector();
  a = new G4DataVector();
  b = new G4DataVector();
  ITTL = new std::vector<size_t>;
  ITTU = new std::vector<size_t>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
G4PenelopeSamplingData::~G4PenelopeSamplingData()
{
  if (x) delete x;
  if (pac) delete pac;
  if (a) delete a;
  if (b) delete b;
  if (ITTL) delete ITTL;
  if (ITTU) delete ITTU;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.......oooOO0OOooo...
size_t G4PenelopeSamplingData::GetNumberOfStoredPoints()
{
  size_t points = x->size();

  //check everything is all right
  if (pac->size() != points || a->size() != points || 
      b->size() != points || ITTL->size() != points ||
      ITTU->size() != points)
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
  if (x) delete x;
  if (pac) delete pac;
  if (a) delete a;
  if (b) delete b;
  if (ITTL) delete ITTL;
  if (ITTU) delete ITTU;
  //create vectors
  x = new G4DataVector();
  pac = new G4DataVector();
  a = new G4DataVector();
  b = new G4DataVector();
  ITTL = new std::vector<size_t>;
  ITTU = new std::vector<size_t>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...
void G4PenelopeSamplingData::AddPoint(G4double x0,G4double pac0,G4double a0,G4double b0,
					size_t ITTL0,size_t ITTU0)
{
  x->push_back(x0);
  pac->push_back(pac0);
  a->push_back(a0);
  b->push_back(b0);
  ITTL->push_back(ITTL0);
  ITTU->push_back(ITTU0);

  //check how many points we do have now
  size_t nOfPoints = GetNumberOfStoredPoints();

  if (nOfPoints > ((size_t)np))
    {
      G4cout << "G4PenelopeSamplingData::AddPoint() " << G4endl;
      G4cout << "WARNING: Up to now there are " << nOfPoints << " points in the table" << G4endl;
      G4cout << "while the anticipated (declared) number is " << np << G4endl;
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
      G4cout << i << " " << (*x)[i] << " " << (*pac)[i] << " " << (*a)[i] << " " << 
	(*b)[i] << " " << (*ITTL)[i] << " " << (*ITTU)[i] << G4endl;
    }
  G4cout << "*************************************************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetX(size_t index)
{
  if (index < x->size())
    return (*x)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetPAC(size_t index)
{
  if (index < pac->size())
    return (*pac)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetA(size_t index)
{
  if (index < a->size())
    return (*a)[index];
  else
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
G4double G4PenelopeSamplingData::GetB(size_t index)
{
  if (index < b->size())
    return (*b)[index];
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
  size_t i = (*ITTL)[itn];
  size_t j = (*ITTU)[itn];

  while ((j-i) > 1)
    {
      size_t k = (i+j)/2;
      if (maxRand > (*pac)[k])
	i = k;
      else
	j = k;
    }

  //Sampling from the rational inverse cumulative distribution
  G4double result = 0;

  G4double rr = maxRand - (*pac)[i];
  if (rr > 1e-16)
    {
      G4double d = (*pac)[i+1]-(*pac)[i];
      result = (*x)[i]+
	((1.0+(*a)[i]+(*b)[i])*d*rr/
	 (d*d+((*a)[i]*d+(*b)[i]*rr)*rr))*((*x)[i+1]-(*x)[i]);      
    }
  else
    result = (*x)[i]; 
  
  return result;
}
