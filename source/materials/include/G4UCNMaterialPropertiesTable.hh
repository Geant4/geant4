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
//---------------------------------------------------------------------------
//
// ClassName:   G4UCNMaterialPropertiesTable
//
// Class description:
//
// A derived class of G4MaterialPropertiesTable in order to save the look-up
// table for the microroughness probability. The derived class has four
// pointers to G4double-arrays:
// (1) integral prob. for reflection
// (2) maximum probability for reflection (needed for accept-reject method)
// (3) integral prob. for transmission
// (4) maximum probability for transmission
//
// 12-05-14, adopted from Stefan Heule (PSI) Thesis by P.Gumplinger
//           http://ucn.web.psi.ch/papers/stefanheule_thesis2008.pdf
//           reported in F. Atchison et al., Eur. Phys. J. A 44, 23–29 (2010)
//           Thanks to Geza Zsigmond

#ifndef G4UCNMATERIALPROPERTIESTABLE_HH
#define G4UCNMATERIALPROPERTIESTABLE_HH 1

#include "G4MaterialPropertiesTable.hh"

class G4UCNMaterialPropertiesTable : public G4MaterialPropertiesTable
{
public:

  G4UCNMaterialPropertiesTable();
  ~G4UCNMaterialPropertiesTable() override;

  // returns the pointer to the mr-reflection table
  G4double* GetMicroRoughnessTable();

  // returns the pointer to the mr-transmission table
  G4double* GetMicroRoughnessTransTable();

  // Assigns double-array to the table-pointers, currently not used
  void LoadMicroRoughnessTables(G4double*, G4double*, G4double*, G4double*);

  // Creates new double arrays and assigns them to the table pointers
  void InitMicroRoughnessTables();

  // Reads the MR-parameters from the corresponding fields and starts
  // the computation of the mr-tables
  void ComputeMicroRoughnessTables();

  // returns the integral prob. value for a theta_i - E pair
  G4double GetMRIntProbability (G4double, G4double);

  // returns the maximum prob. value for a theta_i - E pair
  G4double GetMRMaxProbability (G4double, G4double);

  // sets the maximum prob. value for a theta_i - E pair
  void SetMRMaxProbability (G4double, G4double, G4double);

  // returns the mr-prob.

  // arguments:
  //         1) theta_i
  //         2) Energy
  //         3) V_F
  //         4) theta_o
  //         5) phi_o

  G4double GetMRProbability (G4double, G4double, G4double, G4double, G4double);

  // returns the integral transmission prob. value for a theta_i - E pair
  G4double GetMRIntTransProbability (G4double, G4double);

  // returns the maximum transmission prob. for a theta_i - E pair
  G4double GetMRMaxTransProbability (G4double, G4double);

  // sets the maximum prob. value for a theta_i - E pair
  void SetMRMaxTransProbability (G4double, G4double, G4double);

  // returns the mr-transmission-prob.

  // arguments:
  //         1) theta_i
  //         2) E
  //         3) V_F
  //         4) theta_o
  //         5) phi_o

  G4double GetMRTransProbability (G4double, G4double,
                                  G4double, G4double, G4double);

  // Checks if the validity condition for the microroughness model are
  // satisfied, cf. Steyerl-paper p. 175
  G4bool ConditionsValid (G4double E, G4double VFermi, G4double theta_i);

  // Checks if the validity conditions for the transmission of the
  // microroughness model are satisfied
  G4bool TransConditionsValid (G4double E, G4double VFermi, G4double theta_i);

  // Adds the values for mr-related units to the MaterialPropertiesTable

  // arguments:
  //         1) w
  //         2) b
  //         3) number of angles theta_i in the look-up tables
  //         4) number of energies in the look-up tables
  //         5) minimum value of theta_i
  //         6) maximum value of theta_i
  //         7) minimum value of E
  //         8) maximum value of E
  //         9) number of angles theta_o in the look-up table calculation
  //        10) number of angles phi_o   in the look-up table calculation
  //        11) angular cut

  void SetMicroRoughnessParameters(G4double, G4double, G4int, G4int, G4double,
                                   G4double, G4double, G4double, G4int, G4int,
                                   G4double);

  // returns b
  G4double GetRMS() const;

  // returns w
  G4double GetCorrLen() const;

private:

  // Pointer to the integral reflection probability table
  G4double* theMicroRoughnessTable;

  // Pointer to the maximum reflection probability table
  G4double* maxMicroRoughnessTable;

  // Pointer to the integral transmission probability table
  G4double* theMicroRoughnessTransTable;

  // Pointer to the maximum transmission probability table
  G4double* maxMicroRoughnessTransTable;

  G4double theta_i_min;
  G4double theta_i_max;
  G4double Emin;
  G4double Emax;
  G4int no_theta_i;
  G4int noE;
  G4double theta_i_step;
  G4double E_step;

  // RMS roughness and correlation length
  G4double b, w;
  G4double AngCut;
};

// ==========================================================================
// inline functions
// ==========================================================================

inline G4double G4UCNMaterialPropertiesTable::GetRMS() const {return b;}
inline G4double G4UCNMaterialPropertiesTable::GetCorrLen() const {return w;}

#endif
