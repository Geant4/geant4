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
//  Class:    G4AdjointInterpolator
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Used by  G4AdjointCSManager for interpolation purpose.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointInterpolator_h
#define G4AdjointInterpolator_h 1

#include "globals.hh"

#include <vector>

class G4AdjointInterpolator
{
 public:
  static G4AdjointInterpolator* GetAdjointInterpolator();
  static G4AdjointInterpolator* GetInstance();

  ~G4AdjointInterpolator();

  // Caution: everywhere it is considered that x_vec increases monotically

  G4double LinearInterpolation(G4double& x, G4double& x1, G4double& x2,
                               G4double& y1, G4double& y2);

  G4double LogarithmicInterpolation(G4double& x, G4double& x1, G4double& x2,
                                    G4double& y1, G4double& y2);

  G4double ExponentialInterpolation(G4double& x, G4double& x1, G4double& x2,
                                    G4double& y1, G4double& y2);

  G4double Interpolation(G4double& x, G4double& x1, G4double& x2, G4double& y1,
                         G4double& y2, G4String InterPolMethod = "Log");

  size_t FindPosition(G4double& x, std::vector<G4double>& x_vec,
                      size_t ind_min = 0, size_t ind_max = 0);

  size_t FindPositionForLogVector(G4double& x, std::vector<G4double>& x_vec);

  // xvec should monotically increase
  G4double Interpolate(G4double& x, std::vector<G4double>& x_vec,
                       std::vector<G4double>& y_vec,
                       G4String InterPolMethod = "Log");

  G4double InterpolateWithIndexVector(
    G4double& x, std::vector<G4double>& x_vec, std::vector<G4double>& y_vec,
    std::vector<size_t>& index_vec, G4double x0,
    G4double dx);  // xvec should monotically increase

  G4double InterpolateForLogVector(G4double& x, std::vector<G4double>& x_vec,
                                   std::vector<G4double>& y_vec);

  G4AdjointInterpolator(G4AdjointInterpolator&) = delete;
  G4AdjointInterpolator& operator=(const G4AdjointInterpolator& right) = delete;

 private:
  G4AdjointInterpolator();

  static G4ThreadLocal G4AdjointInterpolator* fInstance;
};
#endif
