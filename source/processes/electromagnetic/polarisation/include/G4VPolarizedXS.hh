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
// File name:     G4VPolarizedXS
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   (virtual) interface class
//
//   provides readable but efficient routines to determine
//   polarization for the final state of a given process
//   empoying the differential cross section

#ifndef G4VPolarizedXS_h
#define G4VPolarizedXS_h 1

#include "G4StokesVector.hh"

class G4VPolarizedXS
{
 public:
  G4VPolarizedXS();
  virtual ~G4VPolarizedXS();

  virtual void Initialize(G4double, G4double, G4double,
                          const G4StokesVector& p0, const G4StokesVector& p1,
                          G4int flag = 0) = 0;

  virtual G4double XSection(const G4StokesVector& pol2,
                            const G4StokesVector& pol3) = 0;

  virtual G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
                                 const G4StokesVector& pol0,
                                 const G4StokesVector& pol1);

  // return expected mean polarisation
  virtual G4StokesVector GetPol2();
  virtual G4StokesVector GetPol3();

  // return basic kinematics properties
  //   minimal gamma value in TotalXSection
  inline G4double GetYmin() { return fYmin; }
  //   minimal energy fraction in TotalXSection
  virtual G4double GetXmin(G4double y);
  //   maximal energy fraction in TotalXSection
  virtual G4double GetXmax(G4double y);

  inline void SetMaterial(G4double A, G4double Z, G4double coul)
  {
    fA    = A;
    fZ    = Z;
    fCoul = coul;
  }

  G4VPolarizedXS& operator=(const G4VPolarizedXS& right) = delete;
  G4VPolarizedXS(const G4VPolarizedXS&)                  = delete;

 protected:
  // define kinematics properties
  inline void SetXmin(G4double xmin) { fXmin = xmin; }
  inline void SetXmax(G4double xmax) { fXmax = xmax; }
  inline void SetYmin(G4double ymin) { fYmin = ymin; }

  // kinematic properties
  G4double fXmin, fXmax, fYmin;
  // material properties
  G4double fA, fZ;
  G4double fCoul;
};

#endif
