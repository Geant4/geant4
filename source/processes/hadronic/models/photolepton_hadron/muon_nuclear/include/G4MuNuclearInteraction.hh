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
// $Id: G4MuNuclearInteraction.hh,v 1.7 2009-03-04 19:09:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuNuclearInteraction physics process ---------
//                by Laszlo Urban, May 1998
// ************************************************************

#ifndef G4MuNuclearInteraction_h
#define G4MuNuclearInteraction_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4PionZero.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParametrizedHadronicVertex.hh"
#include "G4HadronicProcessType.hh"
 
class G4MuNuclearInteraction : public G4VDiscreteProcess
 
{ 
public:
 
  G4MuNuclearInteraction(const G4String& processName = "muNuclear");
 
  virtual ~G4MuNuclearInteraction();

  void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void PreparePhysicsTable(const G4ParticleDefinition& ParticleType);

  virtual void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

  virtual void PrintInfoDefinition() ;

  virtual G4VParticleChange *PostStepDoIt(const G4Track& track,
					  const G4Step& Step  ) ;

protected:

  virtual G4double GetMeanFreePath(const G4Track& track,
				   G4double previousStepSize,
				   G4ForceCondition* condition ) ;
 
  virtual G4double ComputeMeanFreePath( const G4ParticleDefinition* ParticleType,
					G4double KineticEnergy,
					const G4Material* aMaterial);

  void ComputePartialSumSigma(  const G4ParticleDefinition* ParticleType,
				G4double KineticEnergy,
				const G4Material* aMaterial);

  virtual G4double ComputeMicroscopicCrossSection(
				const G4ParticleDefinition* ParticleType,
				G4double KineticEnergy,
				G4double AtomicNumber,
				G4double AtomicMass);

  virtual G4double ComputeDMicroscopicCrossSection(
                                const G4ParticleDefinition* ParticleType,
				G4double KineticEnergy,
				G4double AtomicNumber,
				G4double AtomicMass,
				G4double epsilon);

private:

  G4MuNuclearInteraction & operator=(const G4MuNuclearInteraction &right);
  G4MuNuclearInteraction(const G4MuNuclearInteraction&);

  G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  void MakeSamplingTables( const G4ParticleDefinition* ParticleType );

private:

  G4PhysicsTable* theMeanFreePathTable;
  G4PhysicsTable* theCrossSectionTable ;        

  G4OrderedTable PartialSumSigma;     

  G4double LowestKineticEnergy;  
  G4double HighestKineticEnergy;   
  G4int TotBin;      

  //cut from R.P. Kokoulin
  const G4double CutFixed ;
  // for the atomic weight conversion
  G4double GramPerMole ;

  const G4MuonMinus* theMuonMinus;
  const G4MuonPlus* theMuonPlus;
  const G4PionZero* thePionZero;

  // tables for sampling ..............
  static G4int nzdat,ntdat,NBIN ;
  static G4double zdat[5],adat[5],tdat[8] ;
  static G4double ya[1001],proba[5][8][1001] ;
     
  G4ParametrizedHadronicVertex theHadronicVertex;
};
  
#endif
