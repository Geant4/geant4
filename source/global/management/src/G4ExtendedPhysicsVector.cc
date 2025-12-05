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
// G4ExtendedPhysicsVector class implementation
//
// Author:  V.Ivanchenko 09.09.2025
//
// --------------------------------------------------------------------

#include "G4ExtendedPhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include <iomanip>

// --------------------------------------------------------------
G4ExtendedPhysicsVector::G4ExtendedPhysicsVector(G4PhysicsVector* ptr, G4int nxsec)
  : totalData(ptr)
{
  if (nullptr == ptr)
  {
    if (0 < verboseLevel)
    {
      G4cout << "### G4ExtendedPhysicsVector with undefined G4PhysicsVector " << G4endl;
    }
    totalData = new G4PhysicsFreeVector(false);
  }
  if (nxsec > 1)
  {
    nPartialXS = nxsec - 1;
    dataPartialXS = new std::vector<std::vector<G4float>* >;
    dataPartialXS->resize((std::size_t)nPartialXS, nullptr);
  }
}

// --------------------------------------------------------------
G4ExtendedPhysicsVector::~G4ExtendedPhysicsVector()
{
  if (nPartialXS > 0)
  {
    for (auto const & p : *dataPartialXS) { delete p; }
    delete dataPartialXS;
  }
  delete totalData;
}

// --------------------------------------------------------------------
void G4ExtendedPhysicsVector::SetDataLength(G4int dlength)
{
  // this method may be applied for empty vector only
  totalData->SetDataLength(dlength);
  numberOfNodes = totalData->GetVectorLength();
  if (1 < numberOfNodes)
  {
    idxmax = numberOfNodes - 2;
  }
  else
  {
   if (0 < verboseLevel)
    {
      G4cout << "### G4ExtendedPhysicsVector::SetDataLength numberOfNodes="
             << numberOfNodes << " data structure left empty." << G4endl;
    }
    return;
  }

  if (nPartialXS > 0)
  {
    for (G4int i = 0; i < nPartialXS; ++i)
    {
      (*dataPartialXS)[i] = new std::vector<G4float>(numberOfNodes, 0.f); 
    }
  }
}

// --------------------------------------------------------------------
G4double G4ExtendedPhysicsVector::LogLogValue(const G4double e, std::size_t& idx) const
{
  G4bool interpolation = totalData->CheckIndex(e, idx);

  G4double y1 = (*totalData)[idx];
  if (interpolation)
  {
    G4double e1 = totalData->Energy(idx);
    G4double e2 = totalData->Energy(idx + 1);
    G4double y2 = (*totalData)[idx + 1];
    if (e1 > 0.0 && e2 > e1 && y1 > 0.0 && y2 > 0.0)
    {
      y1 *= G4Exp(G4Log(e/e1) * G4Log(y2/y1) / G4Log(e2/e1));
    }
  }
  return y1;
}

// --------------------------------------------------------------------
void G4ExtendedPhysicsVector::PutPartialXSData(const std::size_t idx, const G4double* y)
{
  if (idx >= numberOfNodes)
  {
    if (0 < verboseLevel)
    {
      G4cout << "### G4ExtendedPhysicsVector::PutPartialXSData(..) idx=" << idx
             << " is out of range " << numberOfNodes 
             << G4endl;
    }
    return;
  }

  // prepare data vector, the last vector is not created
  // because of normalisation of the sum to 1.
  G4float sum = 0.f;
  for (G4int i = 0; i < nPartialXS; ++i)
  {
    sum += (G4float)y[i];
    (*((*dataPartialXS)[i]))[idx] = sum;
  }
  sum += (G4float)y[nPartialXS];
  if (sum > 0.f)
  {
    sum = 1.f/sum;
    for (G4int i = 0; i < nPartialXS; ++i)
    {
      (*((*dataPartialXS)[i]))[idx] *= sum;
    }
  }
}

// --------------------------------------------------------------------
G4int G4ExtendedPhysicsVector::SampleReactionChannel(const G4double e,
                                                     const G4double rand,
                                                     std::size_t& idx) const
{
  if (nPartialXS <= 1) { return 0; }
  G4bool interpolation = totalData->CheckIndex(e, idx);

  G4double e1 = totalData->Energy(idx);
  G4double e2 = e1;
  if (interpolation)
  {
    e2 = totalData->Energy(idx + 1) - e1;
    if (e2 <= 0.0) { interpolation = false; }
  }
  for (G4int i=0; i < nPartialXS; ++i)
  {
    G4double xs = (G4double)(*((*dataPartialXS)[i]))[idx];
    if (interpolation)
    {
      G4double xs2 = (G4double)((*((*dataPartialXS)[i]))[idx + 1]) - xs;
      xs += (e - e1) * xs2 / e2;
    }
    if (xs >= rand) { return i; }  
  }
  return nPartialXS;
}

// --------------------------------------------------------------------
G4int
G4ExtendedPhysicsVector::SampleReactionChannelLogLog(const G4double e,
                                                     const G4double rand,
                                                     std::size_t& idx) const
{
  if (nPartialXS <= 1) { return 0; }
  G4bool interpolation = totalData->CheckIndex(e, idx);

  G4double e1 = totalData->Energy(idx);
  if (e1 <= 0.0) { interpolation = false; }
  G4double e2 = e1;
  if (interpolation)
  {
    e2 = totalData->Energy(idx + 1);
    if (e2 <= e1) { interpolation = false; }
    e2 = G4Log(e2 / e1);
  }
  for (G4int i = 0; i < nPartialXS; ++i)
  {
    G4double xs = (G4double)(*((*dataPartialXS)[i]))[idx];
    if (interpolation && xs > 0.0)
    {
      G4double xs2 = (G4double)((*((*dataPartialXS)[i]))[idx + 1]) - xs;
      if (xs2 > 0.0)
      {
        xs *= G4Exp(G4Log(e / e1) * G4Log(xs2 / xs) / e2);
      }
    }
    if (xs >= rand) { return i; }
  }
  return nPartialXS;
}

// --------------------------------------------------------------
void G4ExtendedPhysicsVector::DumpValues(G4double unitE, G4double unitV) const
{
  G4cout << "====== Data length " << numberOfNodes << " =====" << G4endl;
  // e, partial cumulative normalized x-sections, total x-section
  for (std::size_t i = 0; i < numberOfNodes; ++i)
  {
    G4cout << i << ". " << totalData->Energy(i) / unitE;
    for (G4int j = 0; j < nPartialXS; ++j)
    {
      G4cout << "  " << (*((*dataPartialXS)[i]))[j];
    }
    G4cout << "  " << (*totalData)[i] / unitV << G4endl;
  }
}

//---------------------------------------------------------------
