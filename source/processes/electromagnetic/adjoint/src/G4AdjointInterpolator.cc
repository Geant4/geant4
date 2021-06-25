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

#include "G4AdjointInterpolator.hh"

G4ThreadLocal G4AdjointInterpolator* G4AdjointInterpolator::fInstance = nullptr;

///////////////////////////////////////////////////////
G4AdjointInterpolator* G4AdjointInterpolator::GetAdjointInterpolator()
{
  return GetInstance();
}

///////////////////////////////////////////////////////
G4AdjointInterpolator* G4AdjointInterpolator::GetInstance()
{
  if(!fInstance)
  {
    fInstance = new G4AdjointInterpolator;
  }
  return fInstance;
}

///////////////////////////////////////////////////////
G4AdjointInterpolator::G4AdjointInterpolator() {}

///////////////////////////////////////////////////////
G4AdjointInterpolator::~G4AdjointInterpolator() {}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::LinearInterpolation(G4double& x, G4double& x1,
                                                    G4double& x2, G4double& y1,
                                                    G4double& y2)
{
  G4double res = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
  return res;
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::LogarithmicInterpolation(
  G4double& x, G4double& x1, G4double& x2, G4double& y1, G4double& y2)
{
  if(y1 <= 0. || y2 <= 0. || x1 <= 0.)
    return LinearInterpolation(x, x1, x2, y1, y2);
  G4double B   = std::log(y2 / y1) / std::log(x2 / x1);
  G4double A   = y1 / std::pow(x1, B);
  G4double res = A * std::pow(x, B);
  return res;
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::ExponentialInterpolation(
  G4double& x, G4double& x1, G4double& x2, G4double& y1, G4double& y2)
{
  G4double B   = (std::log(y2) - std::log(y1)) / (x2 - x1);
  G4double A   = y1 * std::exp(-B * x1);
  G4double res = A * std::exp(B * x);
  return res;
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::Interpolation(G4double& x, G4double& x1,
                                              G4double& x2, G4double& y1,
                                              G4double& y2,
                                              G4String InterPolMethod)
{
  if(InterPolMethod == "Log")
  {
    return LogarithmicInterpolation(x, x1, x2, y1, y2);
  }
  else if(InterPolMethod == "Lin")
  {
    return LinearInterpolation(x, x1, x2, y1, y2);
  }
  else if(InterPolMethod == "Exp")
  {
    return ExponentialInterpolation(x, x1, x2, y1, y2);
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "The interpolation method that you invoked does not exist!\n";
    G4Exception("G4AdjointInterpolator::Interpolation", "adoint001",
                FatalException, ed);
    return 0.;
  }
}

///////////////////////////////////////////////////////
// only valid if x_vec is monotically increasing
size_t G4AdjointInterpolator::FindPosition(G4double& x,
                                           std::vector<G4double>& x_vec, size_t,
                                           size_t)
{
  // most rapid method could be used probably

  size_t ndim = x_vec.size();
  size_t ind1 = 0;
  size_t ind2 = ndim - 1;

  if(ndim > 1)
  {
    if(x_vec[0] < x_vec[1])
    {  // increasing
      do
      {
        size_t midBin = (ind1 + ind2) / 2;
        if(x < x_vec[midBin])
          ind2 = midBin;
        else
          ind1 = midBin;
      } while(ind2 - ind1 > 1);
    }
    else
    {
      do
      {
        size_t midBin = (ind1 + ind2) / 2;
        if(x < x_vec[midBin])
          ind1 = midBin;
        else
          ind2 = midBin;
      } while(ind2 - ind1 > 1);
    }
  }

  return ind1;
}

///////////////////////////////////////////////////////
// only valid if x_vec is monotically increasing
size_t G4AdjointInterpolator::FindPositionForLogVector(
  G4double& log_x, std::vector<G4double>& log_x_vec)
{
  // most rapid method could be used probably
  return FindPosition(log_x, log_x_vec);
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::Interpolate(G4double& x,
                                            std::vector<G4double>& x_vec,
                                            std::vector<G4double>& y_vec,
                                            G4String InterPolMethod)
{
  size_t i = FindPosition(x, x_vec);
  return Interpolation(x, x_vec[i], x_vec[i + 1], y_vec[i], y_vec[i + 1],
                       InterPolMethod);
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::InterpolateWithIndexVector(
  G4double& x, std::vector<G4double>& x_vec, std::vector<G4double>& y_vec,
  std::vector<size_t>& index_vec, G4double x0,
  G4double dx)  // only linear interpolation possible
{
  size_t ind = 0;
  if(x > x0)
    ind = int((x - x0) / dx);
  if(ind >= index_vec.size() - 1)
    ind = index_vec.size() - 2;
  size_t ind1 = index_vec[ind];
  size_t ind2 = index_vec[ind + 1];
  if(ind1 > ind2)
  {
    size_t ind11 = ind1;
    ind1         = ind2;
    ind2         = ind11;
  }
  ind = FindPosition(x, x_vec, ind1, ind2);
  return Interpolation(x, x_vec[ind], x_vec[ind + 1], y_vec[ind],
                       y_vec[ind + 1], "Lin");
}

///////////////////////////////////////////////////////
G4double G4AdjointInterpolator::InterpolateForLogVector(
  G4double& log_x, std::vector<G4double>& log_x_vec,
  std::vector<G4double>& log_y_vec)
{
  size_t i = FindPositionForLogVector(log_x, log_x_vec);

  G4double log_y = LinearInterpolation(log_x, log_x_vec[i], log_x_vec[i + 1],
                                       log_y_vec[i], log_y_vec[i + 1]);
  return log_y;
}
