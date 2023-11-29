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
// G4StatDouble
//
// Class description:
//
// Class providing simple "one variable statistics" capability.

// Original Author: Giovanni Santin (ESA) - October 2005 in GRAS tool
// Adaptation and comments by: John Apostolakis (CERN) - November 2011
// --------------------------------------------------------------------
#ifndef G4StatDouble_hh
#define G4StatDouble_hh 1

#include "globals.hh"

class G4StatDouble
{
 public:
  G4StatDouble();
  G4StatDouble(G4double);
  virtual ~G4StatDouble() = default;

  G4StatDouble(const G4StatDouble&) = default;

  G4StatDouble& operator=(const G4double& rhs)
  {
    reset();
    fill(rhs);
    return *this;
  }
  G4StatDouble& operator=(const G4StatDouble& rhs) = default;
  G4StatDouble& operator+=(const G4double& rhs)
  {
    fill(rhs);
    return *this;
  }
  G4StatDouble& operator+=(const G4StatDouble& rhs)
  {
    add(&rhs);
    return *this;
  }

  void reset();
  void fill(G4double x, G4double weight = 1.);
  // Add new data point: value "x" with weight
  void scale(G4double);
  // Reset scale

  G4double mean() const;
  G4double rms();
  // The moments

  G4double mean(G4double ext_sum_w) const;
  // Mean scaled to sum of weights
  G4double rms(G4double ext_sum_w, G4int ext_n);
  // RMS  scaled to sum of weights

  void add(const G4StatDouble*);
  // merge 2 statistics

  inline G4int n() const { return m_n; }
  inline G4double sum_w() const { return m_sum_w; }
  inline G4double sum_w2() const { return m_sum_w2; }
  inline G4double sum_wx() const { return m_sum_wx; }
  inline G4double sum_wx2() const { return m_sum_wx2; }

 protected:
  G4double rms(G4double sum_wx, G4double sum_wx2, G4double sum_w, G4int n);

 protected:
  G4double m_sum_wx  = 0.0;
  G4double m_sum_wx2 = 0.0;
  G4int m_n          = 0;
  G4double m_sum_w   = 0.0;
  G4double m_sum_w2  = 0.0;
  G4double m_scale   = 0.0;
};

#endif
