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
// G4PhysicsVector class implementation
//
// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 03 Mar. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// --------------------------------------------------------------------

#include "G4PhysicsVector.hh"
#include <iomanip>

// --------------------------------------------------------------
G4PhysicsVector::G4PhysicsVector(G4bool val)
  : useSpline(val)
{}

// --------------------------------------------------------------------
void G4PhysicsVector::Initialise()
{
  idxmax = numberOfNodes - 2;
  if(0 < numberOfNodes)
  {
    edgeMin = binVector[0];
    edgeMax = binVector[numberOfNodes - 1];
  }
}

// --------------------------------------------------------------
G4bool G4PhysicsVector::Store(std::ofstream& fOut, G4bool ascii) const
{
  // Ascii mode
  if(ascii)
  {
    fOut << *this;
    return true;
  }
  // Binary Mode

  // binning
  fOut.write((char*) (&edgeMin), sizeof edgeMin);
  fOut.write((char*) (&edgeMax), sizeof edgeMax);
  fOut.write((char*) (&numberOfNodes), sizeof numberOfNodes);

  // contents
  std::size_t size = dataVector.size();
  fOut.write((char*) (&size), sizeof size);

  G4double* value = new G4double[2 * size];
  for(std::size_t i = 0; i < size; ++i)
  {
    value[2 * i]     = binVector[i];
    value[2 * i + 1] = dataVector[i];
  }
  fOut.write((char*) (value), 2 * size * (sizeof(G4double)));
  delete[] value;

  return true;
}

// --------------------------------------------------------------
G4bool G4PhysicsVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  // clear properties;
  dataVector.clear();
  binVector.clear();
  secDerivative.clear();

  // retrieve in ascii mode
  if(ascii)
  {
    // binning
    fIn >> edgeMin >> edgeMax >> numberOfNodes;
    if(fIn.fail() || numberOfNodes < 2)
    {
      return false;
    }
    // contents
    G4int siz = 0;
    fIn >> siz;
    if(fIn.fail() || siz != G4int(numberOfNodes))
    {
      return false;
    }

    binVector.reserve(siz);
    dataVector.reserve(siz);
    G4double vBin, vData;

    for(G4int i = 0; i < siz; ++i)
    {
      vBin  = 0.;
      vData = 0.;
      fIn >> vBin >> vData;
      if(fIn.fail())
      {
        return false;
      }
      binVector.push_back(vBin);
      dataVector.push_back(vData);
    }
    Initialise();
    return true;
  }

  // retrieve in binary mode
  // binning
  fIn.read((char*) (&edgeMin), sizeof edgeMin);
  fIn.read((char*) (&edgeMax), sizeof edgeMax);
  fIn.read((char*) (&numberOfNodes), sizeof numberOfNodes);

  // contents
  std::size_t size;
  fIn.read((char*) (&size), sizeof size);

  G4double* value = new G4double[2 * size];
  fIn.read((char*) (value), 2 * size * (sizeof(G4double)));
  if(G4int(fIn.gcount()) != G4int(2 * size * (sizeof(G4double))))
  {
    delete[] value;
    return false;
  }

  binVector.reserve(size);
  dataVector.reserve(size);
  for(std::size_t i = 0; i < size; ++i)
  {
    binVector.push_back(value[2 * i]);
    dataVector.push_back(value[2 * i + 1]);
  }
  delete[] value;

  Initialise();
  return true;
}

// --------------------------------------------------------------
void G4PhysicsVector::DumpValues(G4double unitE, G4double unitV) const
{
  for(std::size_t i = 0; i < numberOfNodes; ++i)
  {
    G4cout << binVector[i] / unitE << "   " << dataVector[i] / unitV 
           << G4endl;
  }
}

// --------------------------------------------------------------------
std::size_t G4PhysicsVector::FindBin(const G4double energy, 
                                     std::size_t idx) const
{
  if(idx + 1 < numberOfNodes && 
     energy >= binVector[idx] && energy <= binVector[idx])
  {
    return idx;
  } 
  if(energy <= binVector[1])
  {
    return 0;
  }
  if(energy >= binVector[idxmax])
  {
    return idxmax;
  }
  return GetBin(energy); 
}

// --------------------------------------------------------------------
void G4PhysicsVector::ScaleVector(const G4double factorE, 
                                  const G4double factorV)
{
  for(std::size_t i = 0; i < numberOfNodes; ++i)
  {
    binVector[i] *= factorE;
    dataVector[i] *= factorV;
  }
  Initialise();
}

// --------------------------------------------------------------------
void G4PhysicsVector::FillSecondDerivatives(const G4SplineType stype,
					    const G4double dir1,
					    const G4double dir2)
{
  if(!useSpline) { return; }
  // cannot compute derivatives for less than 5 points
  const std::size_t nmin = (stype == G4SplineType::Base) ? 5 : 4;
  if(nmin > numberOfNodes) 
  {
    if(0 < verboseLevel)
    { 
      G4cout << "### G4PhysicsVector: spline cannot be used for "
	     << numberOfNodes << " points - spline disabled" 
	     << G4endl;
      DumpValues();
    }
    useSpline = false;
    return;
  }
  // check energies of free vector
  if(type == T_G4PhysicsFreeVector)
  {
    for(std::size_t i=0; i<=idxmax; ++i) 
    {
      if(binVector[i + 1] <= binVector[i])
      {
        if(0 < verboseLevel) 
        {
	  G4cout << "### G4PhysicsVector: spline cannot be used, because "
		 << " E[" << i << "]=" << binVector[i]
		 << " >= E[" << i+1 << "]=" << binVector[i + 1]
		 << G4endl;
	  DumpValues();
        }
        useSpline = false;
        return;
      }
    }
  }

  // spline is possible
  Initialise();
  secDerivative.resize(numberOfNodes);

  if(1 < verboseLevel)
  {
    G4cout << "### G4PhysicsVector:: FillSecondDerivatives N=" 
           << numberOfNodes << G4endl;
    DumpValues();
  }

  switch(stype) 
  {
    case G4SplineType::Base:
      ComputeSecDerivative1();
      break;

    case G4SplineType::FixedEdges:
      ComputeSecDerivative2(dir1, dir2);
      break;

    default:
      ComputeSecDerivative0();
  }
}

// --------------------------------------------------------------
void G4PhysicsVector::ComputeSecDerivative0()
//  A simplified method of computation of second derivatives
{
  std::size_t n = numberOfNodes - 1;

  for(std::size_t i = 1; i < n; ++i)
  {
    secDerivative[i] = 3.0 *
      ((dataVector[i + 1] - dataVector[i]) / (binVector[i + 1] - binVector[i]) -
       (dataVector[i] - dataVector[i - 1]) /
         (binVector[i] - binVector[i - 1])) /
      (binVector[i + 1] - binVector[i - 1]);
  }
  secDerivative[n] = secDerivative[n - 1];
  secDerivative[0] = secDerivative[1];
}

// --------------------------------------------------------------
void G4PhysicsVector::ComputeSecDerivative1()
// Computation of second derivatives using "Not-a-knot" endpoint conditions
// B.I. Kvasov "Methods of shape-preserving spline approximation"
// World Scientific, 2000
{
  std::size_t n = numberOfNodes - 1;
  G4double* u = new G4double[n];
  G4double p, sig;

  u[1] = ((dataVector[2] - dataVector[1]) / (binVector[2] - binVector[1]) -
          (dataVector[1] - dataVector[0]) / (binVector[1] - binVector[0]));
  u[1] = 6.0 * u[1] * (binVector[2] - binVector[1]) /
         ((binVector[2] - binVector[0]) * (binVector[2] - binVector[0]));

  // Decomposition loop for tridiagonal algorithm. secDerivative[i]
  // and u[i] are used for temporary storage of the decomposed factors.

  secDerivative[1] = (2.0 * binVector[1] - binVector[0] - binVector[2]) /
                     (2.0 * binVector[2] - binVector[0] - binVector[1]);

  for(std::size_t i = 2; i < n - 1; ++i)
  {
    sig =
      (binVector[i] - binVector[i - 1]) / (binVector[i + 1] - binVector[i - 1]);
    p                = sig * secDerivative[i - 1] + 2.0;
    secDerivative[i] = (sig - 1.0) / p;
    u[i] =
      (dataVector[i + 1] - dataVector[i]) / (binVector[i + 1] - binVector[i]) -
      (dataVector[i] - dataVector[i - 1]) / (binVector[i] - binVector[i - 1]);
    u[i] =
      (6.0 * u[i] / (binVector[i + 1] - binVector[i - 1])) - sig * u[i - 1] / p;
  }

  sig =
    (binVector[n - 1] - binVector[n - 2]) / (binVector[n] - binVector[n - 2]);
  p = sig * secDerivative[n - 3] + 2.0;
  u[n - 1] =
    (dataVector[n] - dataVector[n - 1]) / (binVector[n] - binVector[n - 1]) -
    (dataVector[n - 1] - dataVector[n - 2]) /
      (binVector[n - 1] - binVector[n - 2]);
  u[n - 1] = 6.0 * sig * u[n - 1] / (binVector[n] - binVector[n - 2]) -
             (2.0 * sig - 1.0) * u[n - 2] / p;

  p = (1.0 + sig) + (2.0 * sig - 1.0) * secDerivative[n - 2];
  secDerivative[n - 1] = u[n - 1] / p;

  // The back-substitution loop for the triagonal algorithm of solving
  // a linear system of equations.

  for(std::size_t k = n - 2; k > 1; --k)
  {
    secDerivative[k] *=
      (secDerivative[k + 1] - u[k] * (binVector[k + 1] - binVector[k - 1]) /
                                (binVector[k + 1] - binVector[k]));
  }
  secDerivative[n] =
    (secDerivative[n - 1] - (1.0 - sig) * secDerivative[n - 2]) / sig;
  sig = 1.0 - ((binVector[2] - binVector[1]) / (binVector[2] - binVector[0]));
  secDerivative[1] *= (secDerivative[2] - u[1] / (1.0 - sig));
  secDerivative[0] = (secDerivative[1] - sig * secDerivative[2]) / (1.0 - sig);

  delete[] u;
}

// --------------------------------------------------------------
void G4PhysicsVector::ComputeSecDerivative2(G4double firstPointDerivative,
                                            G4double endPointDerivative)
// A standard method of computation of second derivatives
// First derivatives at the first and the last point should be provided
// See for example W.H. Press et al. "Numerical recipes in C"
// Cambridge University Press, 1997.
{
  std::size_t n = numberOfNodes - 1;
  G4double* u = new G4double[n];
  G4double p, sig, un;

  u[0] = (6.0 / (binVector[1] - binVector[0])) *
         ((dataVector[1] - dataVector[0]) / (binVector[1] - binVector[0]) -
          firstPointDerivative);

  secDerivative[0] = -0.5;

  // Decomposition loop for tridiagonal algorithm. secDerivative[i]
  // and u[i] are used for temporary storage of the decomposed factors.

  for(std::size_t i = 1; i < n; ++i)
  {
    sig =
      (binVector[i] - binVector[i - 1]) / (binVector[i + 1] - binVector[i - 1]);
    p                = sig * (secDerivative[i - 1]) + 2.0;
    secDerivative[i] = (sig - 1.0) / p;
    u[i] =
      (dataVector[i + 1] - dataVector[i]) / (binVector[i + 1] - binVector[i]) -
      (dataVector[i] - dataVector[i - 1]) / (binVector[i] - binVector[i - 1]);
    u[i] =
      6.0 * u[i] / (binVector[i + 1] - binVector[i - 1]) - sig * u[i - 1] / p;
  }

  sig =
    (binVector[n - 1] - binVector[n - 2]) / (binVector[n] - binVector[n - 2]);
  p  = sig * secDerivative[n - 2] + 2.0;
  un = (6.0 / (binVector[n] - binVector[n - 1])) *
         (endPointDerivative - (dataVector[n] - dataVector[n - 1]) /
                                 (binVector[n] - binVector[n - 1])) -
       u[n - 1] / p;
  secDerivative[n] = un / (secDerivative[n - 1] + 2.0);

  // The back-substitution loop for the triagonal algorithm of solving
  // a linear system of equations.

  for(std::size_t k = n - 1; k > 0; --k)
  {
    secDerivative[k] *=
      (secDerivative[k + 1] - u[k] * (binVector[k + 1] - binVector[k - 1]) /
                                (binVector[k + 1] - binVector[k]));
  }
  secDerivative[0] = 0.5 * (u[0] - secDerivative[1]);

  delete[] u;
}

// --------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4PhysicsVector& pv)
{
  // binning
  G4long prec = out.precision();
  out << std::setprecision(12) << pv.edgeMin << " " << pv.edgeMax << " "
      << pv.numberOfNodes << G4endl;

  // contents
  out << pv.dataVector.size() << G4endl;
  for(std::size_t i = 0; i < pv.dataVector.size(); ++i)
  {
    out << pv.binVector[i] << "  " << pv.dataVector[i] << G4endl;
  }
  out.precision(prec);

  return out;
}

//---------------------------------------------------------------
G4double G4PhysicsVector::GetEnergy(const G4double val) const
{
  if(0 == numberOfNodes)
  {
    return 0.0;
  }
  if(1 == numberOfNodes || val <= dataVector[0])
  {
    return edgeMin;
  }
  if(val >= dataVector[numberOfNodes - 1])
  {
    return edgeMax;
  }
  std::size_t bin = std::lower_bound(dataVector.cbegin(), dataVector.cend(), val)
                  - dataVector.cbegin() - 1;
  if(bin > idxmax) { bin = idxmax; } 
  G4double res = binVector[bin];
  G4double del = dataVector[bin + 1] - dataVector[bin];
  if(del > 0.0)
  {
    res += (val - dataVector[bin]) * (binVector[bin + 1] - res) / del;
  }
  return res;
}

//---------------------------------------------------------------
void G4PhysicsVector::PrintPutValueError(std::size_t index, 
                                         G4double val, 
                                         const G4String& text)
{
  G4ExceptionDescription ed;
  ed << "Vector type: " << type << " length= " << numberOfNodes
     << "; an attempt to put data at index= " << index
     << " value= " << val << " in " << text;
  G4Exception("G4PhysicsVector:", "gl0005", 
              FatalException, ed, "Wrong operation");
}

//---------------------------------------------------------------
