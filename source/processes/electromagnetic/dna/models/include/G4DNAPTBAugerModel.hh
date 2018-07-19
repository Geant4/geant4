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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
// $Id: G4DNAPTBAugerModel.cc,v 1.11
// 
//
// -------------------------------------------------------------------
//
// Geant4 Header G4DNAPTBAugerModel
//  
// -------------------------------------------------------------------
//
// Class description:
// Implementation of atomic deexcitation 
//
// -------------------------------------------------------------------

#ifndef G4DNAPTBAugerModel_h
#define G4DNAPTBAugerModel_h 1

#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include <vector>

class G4AtomicTransitionManager;
class G4VhShellCrossSection;
class G4EmCorrections;
class G4Material;

/*!
 * \brief The G4DNAPTBAugerModel class
 * Implement the PTB Auger model
 */
class G4DNAPTBAugerModel
{  
public: 
  
   /*!
   * \brief G4DNAPTBAugerModel
   * Constructor
   * \param modelName
   */
  G4DNAPTBAugerModel(const G4String &modelName);

  /*!
   * \brief ~G4DNAPTBAugerModel
   * Destructor
   */
  virtual ~G4DNAPTBAugerModel();
   

  /*!
   * \brief Initialise
   * Set the verbose value
   */
  virtual void Initialise();

  /*!
   * \brief SetCutForAugerElectrons
   * Set the cut for the auger electrons production
   * \param cut
   */
  void SetCutForAugerElectrons(G4double cut);

  /*!
   * \brief ComputeAugerEffect
   * Main method to be called by the ionisation model.
   * \param fvect
   * \param materialNameIni
   * \param bindingEnergy
   */
  void ComputeAugerEffect(std::vector<G4DynamicParticle *> *fvect, const G4String& materialNameIni, G4double bindingEnergy);

private:
 
  const G4String modelName; ///< name of the auger model

  G4int verboseLevel;
  G4double minElectronEnergy;

  /*!
   * \brief GenerateAugerWithRandomDirection
   * Generates the auger particle
   * \param fvect
   * \param kineticEnergy
   */
  void GenerateAugerWithRandomDirection(std::vector<G4DynamicParticle*>* fvect, G4double kineticEnergy);

  /*!
   * \brief CalculAugerEnergyFor
   * \param atomId
   * \return the auger particle energy
   */
  G4double CalculAugerEnergyFor(G4int atomId);

  /*!
   * \brief DetermineIonisedAtom
   * \param atomId
   * \param materialName
   * \param bindingEnergy
   * \return the id of the chosen ionised atom
   */
  G4int DetermineIonisedAtom(G4int atomId, const G4String &materialName, G4double bindingEnergy);

  // copy constructor and hide assignment operator
  G4DNAPTBAugerModel(G4DNAPTBAugerModel &);  // prevent copy-construction
  G4DNAPTBAugerModel & operator=(const G4DNAPTBAugerModel &right);  // prevent assignement
 
};

#endif




