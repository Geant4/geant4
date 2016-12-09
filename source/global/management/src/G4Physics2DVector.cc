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
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4Physics2DVector.cc
//
//  Author:        Vladimir Ivanchenko
//
//  Creation date: 25.09.2011
//
// --------------------------------------------------------------

#include <iomanip>

#include "G4Physics2DVector.hh"

// --------------------------------------------------------------

G4Physics2DVector::G4Physics2DVector()
 : type(T_G4PhysicsFreeVector),
   numberOfXNodes(0), numberOfYNodes(0),
   verboseLevel(0), useBicubic(false)
{}

// --------------------------------------------------------------

G4Physics2DVector::G4Physics2DVector(size_t nx, size_t ny)
 : type(T_G4PhysicsFreeVector),
   numberOfXNodes(nx), numberOfYNodes(ny),
   verboseLevel(0), useBicubic(false)
{
  PrepareVectors();
}

// --------------------------------------------------------------

G4Physics2DVector::~G4Physics2DVector() 
{
  ClearVectors();
}

// --------------------------------------------------------------

G4Physics2DVector::G4Physics2DVector(const G4Physics2DVector& right)
{
  type         = right.type;

  numberOfXNodes = right.numberOfXNodes;
  numberOfYNodes = right.numberOfYNodes;

  verboseLevel = right.verboseLevel;
  useBicubic   = right.useBicubic;

  xVector      = right.xVector;
  yVector      = right.yVector;

  PrepareVectors();
  CopyData(right);
}

// --------------------------------------------------------------

G4Physics2DVector& G4Physics2DVector::operator=(const G4Physics2DVector& right)
{
  if (&right==this)  { return *this; }
  ClearVectors();

  type         = right.type;

  numberOfXNodes = right.numberOfXNodes;
  numberOfYNodes = right.numberOfYNodes;

  verboseLevel = right.verboseLevel;
  useBicubic   = right.useBicubic;

  PrepareVectors();
  CopyData(right);

  return *this;
}

// --------------------------------------------------------------

void G4Physics2DVector::PrepareVectors()
{
  xVector.resize(numberOfXNodes,0.);
  yVector.resize(numberOfYNodes,0.);
  value.resize(numberOfYNodes,0);
  for(size_t j=0; j<numberOfYNodes; ++j) {
    G4PV2DDataVector* v = new G4PV2DDataVector();
    v->resize(numberOfXNodes,0.);
    value[j] = v;
  }
}

// --------------------------------------------------------------

void G4Physics2DVector::ClearVectors()
{
  for(size_t j=0; j<numberOfYNodes; ++j) {
    delete value[j];
  }
}

// --------------------------------------------------------------

void G4Physics2DVector::CopyData(const G4Physics2DVector &right)
{
  for(size_t i=0; i<numberOfXNodes; ++i) {
    xVector[i] = right.xVector[i];
  }
  for(size_t j=0; j<numberOfYNodes; ++j) {
    yVector[j] = right.yVector[j];
    G4PV2DDataVector* v0 = right.value[j];
    for(size_t i=0; i<numberOfXNodes; ++i) { 
      PutValue(i,j,(*v0)[i]); 
    }
  }
}

// --------------------------------------------------------------

G4double G4Physics2DVector::Value(G4double xx, G4double yy, 
				  size_t& idx, size_t& idy) const
{
  G4double x = xx;
  G4double y = yy;

  // no interpolation outside the table
  if(x < xVector[0]) { 
    x = xVector[0]; 
  } else if(x > xVector[numberOfXNodes - 1]) { 
    x = xVector[numberOfXNodes - 1]; 
  }
  if(y < yVector[0]) { 
    y = yVector[0]; 
  } else if(y > yVector[numberOfYNodes - 1]) { 
    y = yVector[numberOfYNodes - 1]; 
  }

  // find bins
  idx = FindBinLocationX(x, idx);
  idy = FindBinLocationY(y, idy);

  // interpolate
  if(useBicubic) {
    return BicubicInterpolation(x, y, idx, idy);
  } else {
    G4double x1 = xVector[idx];
    G4double x2 = xVector[idx+1];
    G4double y1 = yVector[idy];
    G4double y2 = yVector[idy+1];
    G4double v11= GetValue(idx,   idy);
    G4double v12= GetValue(idx+1, idy);
    G4double v21= GetValue(idx,   idy+1);
    G4double v22= GetValue(idx+1, idy+1);
    return ((y2 - y)*(v11*(x2 - x) + v12*(x - x1)) + 
	    ((y - y1)*(v21*(x2 - x) + v22*(x - x1))))/((x2 - x1)*(y2 - y1)); 
  }
}

// --------------------------------------------------------------

G4double 
G4Physics2DVector::BicubicInterpolation(G4double x, G4double y,
					size_t idx, size_t idy) const
{
    // Bicubic interpolation according to 
    // 1. H.M. Antia, "Numerical Methods for Scientists and Engineers",
    //    MGH, 1991. 
    // 2. W.H. Press et al., "Numerical recipes. The Art of Scientific 
    //    Computing", Cambridge University Press, 2007. 
    G4double x1 = xVector[idx];
    G4double x2 = xVector[idx+1];
    G4double y1 = yVector[idy];
    G4double y2 = yVector[idy+1];
    G4double f1 = GetValue(idx,   idy);
    G4double f2 = GetValue(idx+1, idy);
    G4double f3 = GetValue(idx+1, idy+1);
    G4double f4 = GetValue(idx,   idy+1);

    G4double dx = x2 - x1;
    G4double dy = y2 - y1;

    G4double h1 = (x - x1)/dx;
    G4double h2 = (y - y1)/dy;   

    G4double h12 = h1*h1;
    G4double h13 = h12*h1;
    G4double h22 = h2*h2;
    G4double h23 = h22*h2;

    // Three derivatives at each of four points (1-4) defining the 
    // subregion are computed by numerical centered differencing from 
    // the functional values already tabulated on the grid. 

    G4double f1x = DerivativeX(idx, idy, dx);
    G4double f2x = DerivativeX(idx+1, idy, dx);
    G4double f3x = DerivativeX(idx+1, idy+1, dx);
    G4double f4x = DerivativeX(idx, idy+1, dx);

    G4double f1y = DerivativeY(idx, idy, dy);
    G4double f2y = DerivativeY(idx+1, idy, dy);
    G4double f3y = DerivativeY(idx+1, idy+1, dy);
    G4double f4y = DerivativeY(idx, idy+1, dy);

    G4double dxy = dx*dy;
    G4double f1xy = DerivativeXY(idx, idy, dxy);
    G4double f2xy = DerivativeXY(idx+1, idy, dxy);
    G4double f3xy = DerivativeXY(idx+1, idy+1, dxy);
    G4double f4xy = DerivativeXY(idx, idy+1, dxy);

    return 
      f1 + f1y*h2 + (3*(f4-f1) - 2*f1y - f4y)*h22 + (2*(f1 - f4) + f1y + f4y)*h23
      + f1x*h1 + f1xy*h1*h2 +(3*(f4x - f1x) - 2*f1xy - f4xy)*h1*h22
      + (2*(f1x - f4x) + f1xy + f4xy)*h1*h23
      + (3*(f2 - f1) - 2*f1x - f2x)*h12 + (3*f2y - 3*f1y - 2*f1xy - f2xy)*h12*h2
      + (9*(f1 - f2 + f3 - f4) + 6*f1x + 3*f2x - 3*f3x - 6*f4x + 6*f1y - 6*f2y
	 - 3*f3y + 3*f4y + 4*f1xy + 2*f2xy + f3xy + 2*f4xy)*h12*h22
      + (6*(-f1 + f2 - f3 + f4) - 4*f1x - 2*f2x + 2*f3x + 4*f4x - 3*f1y
	 + 3*f2y + 3*f3y - 3*f4y - 2*f1xy - f2xy - f3xy - 2*f4xy)*h12*h23
      + (2*(f1 - f2) + f1x + f2x)*h13 + (2*(f1y - f2y) + f1xy + f2xy)*h13*h2
      + (6*(-f1 + f2 -f3 + f4) + 3*(-f1x - f2x + f3x + f4x) - 4*f1y
	 + 4*f2y + 2*f3y - 2*f4y - 2*f1xy - 2*f2xy - f3xy - f4xy)*h13*h22
      + (4*(f1 - f2 + f3 - f4) + 2*(f1x + f2x - f3x - f4x) 
	 + 2*(f1y - f2y - f3y + f4y) + f1xy + f2xy + f3xy + f4xy)*h13*h23;
}

// --------------------------------------------------------------

void 
G4Physics2DVector::PutVectors(const std::vector<G4double>& vecX,
			      const std::vector<G4double>& vecY)
{
  ClearVectors();
  numberOfXNodes = vecX.size();
  numberOfYNodes = vecY.size();
  PrepareVectors();
  for(size_t i = 0; i<numberOfXNodes; ++i) {
    xVector[i] = vecX[i];
  }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    yVector[j] = vecY[j];
  }
}

// --------------------------------------------------------------

void G4Physics2DVector::Store(std::ofstream& out) const
{
  // binning
  G4int prec = out.precision();
  out << G4int(type) << " " << numberOfXNodes << " " << numberOfYNodes 
      << G4endl; 
  out << std::setprecision(5);

  // contents
  for(size_t i = 0; i<numberOfXNodes-1; ++i) {
    out << xVector[i] << "  ";
  }
  out << xVector[numberOfXNodes-1] << G4endl;
  for(size_t j = 0; j<numberOfYNodes-1; ++j) {
    out << yVector[j] << "  ";
  }
  out << yVector[numberOfYNodes-1] << G4endl;
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    for(size_t i = 0; i<numberOfXNodes-1; ++i) {
      out << GetValue(i, j) << "  ";
    }
    out << GetValue(numberOfXNodes-1,j) << G4endl;
  }
  out.precision(prec);
  out.close();
}

// --------------------------------------------------------------

G4bool G4Physics2DVector::Retrieve(std::ifstream& in)
{
  // initialisation
  ClearVectors();

  // binning
  G4int k;
  in >> k >> numberOfXNodes >> numberOfYNodes;
  if (in.fail() || 0 >= numberOfXNodes || 0 >= numberOfYNodes || 
      numberOfXNodes >= INT_MAX || numberOfYNodes >= INT_MAX) { 
    if( 0 >= numberOfXNodes || numberOfXNodes >= INT_MAX) {
      numberOfXNodes = 0;
    }
    if( 0 >= numberOfYNodes || numberOfYNodes >= INT_MAX) {
      numberOfYNodes = 0;
    }
    return false; 
  }
  PrepareVectors();
  type = G4PhysicsVectorType(k); 

  // contents
  G4double val;
  for(size_t i = 0; i<numberOfXNodes; ++i) {
    in >> xVector[i];
    if (in.fail())  { return false; }
  }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    in >> yVector[j];
    if (in.fail())  { return false; }
  }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    for(size_t i = 0; i<numberOfXNodes; ++i) {
      in >> val;
      if (in.fail())  { return false; }
      PutValue(i, j, val);
    }
  }
  in.close();
  return true;
}

// --------------------------------------------------------------

void 
G4Physics2DVector::ScaleVector(G4double factor)
{
  G4double val;
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    for(size_t i = 0; i<numberOfXNodes; ++i) {
      val = GetValue(i, j)*factor;
      PutValue(i, j, val);
    }
  }
}

// --------------------------------------------------------------

size_t 
G4Physics2DVector::FindBinLocation(G4double z, 
				   const G4PV2DDataVector& v) const
{
  size_t bin;
  size_t binmax = v.size() - 2; 

  if(z <= v[0])           { bin = 0; }
  else if(z >= v[binmax]) { bin = binmax; }
  else {
    bin = std::lower_bound(v.begin(), v.end(), z) - v.begin() - 1;
  }
  return bin;
}

// --------------------------------------------------------------

G4double G4Physics2DVector::FindLinearX(G4double rand, G4double yy, 
					size_t& idy) const
{
  G4double y = yy;

  // no interpolation outside the table
  if(y < yVector[0]) { 
    y = yVector[0]; 
  } else if(y > yVector[numberOfYNodes - 1]) { 
    y = yVector[numberOfYNodes - 1]; 
  }

  // find bins
  idy = FindBinLocationY(y, idy);

  G4double x1 = InterpolateLinearX(*(value[idy]), rand);
  G4double x2 = InterpolateLinearX(*(value[idy+1]), rand);
  G4double res = x1;
  G4double del = yVector[idy+1] - yVector[idy];
  if(del != 0.0) {
    res += (x2 - x1)*(y - yVector[idy])/del;
  }
  return res;
}

// --------------------------------------------------------------

G4double G4Physics2DVector::InterpolateLinearX(G4PV2DDataVector& v, 
					       G4double rand) const
{
  size_t nn = v.size();
  if(1 >= nn) { return 0.0; }
  size_t n1 = 0;
  size_t n2 = nn/2;
  size_t n3 = nn - 1;
  G4double y = rand*v[n3];
  while (n1 + 1 != n3)
    {
      if (y > v[n2])
	{ n1 = n2; }
      else
	{ n3 = n2; }
      n2 = (n3 + n1 + 1)/2;
    }
  G4double res = xVector[n1];
  G4double del = v[n3] - v[n1];
  if(del > 0.0) {
    res += (y - v[n1])*(xVector[n3] - res)/del;
  }
  return res;
}

// --------------------------------------------------------------
