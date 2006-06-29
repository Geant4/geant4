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
// $Id: G4PenelopeRayleigh.hh,v 1.5 2006-06-29 19:36:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 01 Dec 2002   L.Pandola      1st implementation
// 14 Feb 2003   MG Pia         Removed data member material
// 18 Mar 2004   M.Mendenhall   Introduced SamplingTable 
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic process, Rayleigh effect
// Penelope model
// -------------------------------------------------------------------

#ifndef G4PENELOPERAYLEIGH_HH
#define G4PENELOPERAYLEIGH_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PenelopeIntegrator.hh"
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4Material;
class G4DataVector;

class G4PenelopeRayleigh : public G4VDiscreteProcess {

public:

  typedef std::map<const G4Material *,std::pair<G4DataVector *, G4DataVector *> >
  SamplingTablePair;

  G4PenelopeRayleigh(const G4String& processName ="PenRayleigh");
  
  ~G4PenelopeRayleigh();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  G4double MolecularFormFactor(G4double x); 
  G4double DifferentialCrossSection (G4double e);
 
  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

 

protected:

  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

private:

  // Hide copy constructor and assignment operator as private 
  G4PenelopeRayleigh& operator=(const G4PenelopeRayleigh &right);
  G4PenelopeRayleigh(const G4PenelopeRayleigh& );
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  //void InizialiseSampling(const G4Material* material);
  void InizialiseSampling();

  G4DataVector* samplingFunction_x;
  G4DataVector* samplingFunction_y;

  SamplingTablePair SamplingTables;

  G4double samplingConstant;
  G4double facte; //cross section factor

  G4VEMDataSet* meanFreePathTable;
  
  const G4Material* material;
  const G4int nBins;
  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;
};

#endif

 












