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
// $Id: G4ParticleChangeForTransport.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// ------------------------------------------------------------
//   Implemented for the new scheme                 10 May. 1998  H.Kurahige
//   Added theMaterialChange                        16 FEb. 2000  H.Kurahige
//   Remove thePolarizationChange		    12 Feb. 2001  H.Kurashige
//   Modification for G4TouchableHandle             22 Oct. 2001  R.Chytracek
//   Add MaterialCutsCouple                          8 Oct. 2002  H.Kurashige
//
// Class Description
//  This class is a concrete class for ParticleChange for transportation
//        
#ifndef G4ParticleChangeForTransport_h
#define G4ParticleChangeForTransport_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"

class G4MaterialCutsCouple;
class G4VSensitiveDetector;

class G4ParticleChangeForTransport: public G4ParticleChange
{ 
  public:
    // default constructor
    G4ParticleChangeForTransport();

    // destructor
    virtual ~G4ParticleChangeForTransport();

  protected:
    // hide copy constructor and assignment operator as protected
    G4ParticleChangeForTransport(const G4ParticleChangeForTransport &right);
    G4ParticleChangeForTransport & operator=(const G4ParticleChangeForTransport &right);

  public: // with description
    // ----------------------------------------------------
    // --- the following methods are for updating G4Step -----   
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process    
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForAtRest(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // A physics process gives the final state of the particle 
    // based on information of G4Track (or equivalently the PreStepPoint)
 
    virtual void Initialize(const G4Track&);
    // Initialize all propoerties by using G4Track information
           
    // ----------------------------------------------------
    //--- methods to keep information of the final state--
    //  IMPORTANT NOTE: Although the name of the class and methods are
    //   "Change", what it stores (and returns in get) are the "FINAL" 
    //   values of the Position, Momentum, etc.

    const G4TouchableHandle& GetTouchableHandle() const;
    void  SetTouchableHandle(const G4TouchableHandle& fTouchable);
    //  Get/Set the touchable of the current particle.
    //  Note: Touchable in PostStepPoint will be updated only after PostStepDoIt

    G4Material* GetMaterialInTouchable() const;
    void SetMaterialInTouchable(G4Material* fMaterial);
    //  Get/Propose the material in the touchable of the current particle.

    const G4MaterialCutsCouple* GetMaterialCutsCoupleInTouchable() const;
    void SetMaterialCutsCoupleInTouchable(const G4MaterialCutsCouple* fMaterialCutsCouple);
    //  Get/Set the materialCutsCouple in the touchable of the current particle.

    G4VSensitiveDetector* GetSensitiveDetectorInTouchable() const;
    void SetSensitiveDetectorInTouchable(G4VSensitiveDetector* fSensitiveDetector);
    //  Get/Set the sensitive detector in the touchable of the current particle.

    G4bool GetMomentumChanged() const;
    void SetMomentumChanged(G4bool b);

  public:
    virtual void DumpInfo() const;

  protected:
    G4TouchableHandle theTouchableHandle;
    //  The changed touchable of a given particle.

  public:

    // Prototype implementation of smooth representation of curved trajectories.
    // Auxiliary points are ThreeVectors for now; change to G4AuxiliaryPoints.

    inline void SetPointerToVectorOfAuxiliaryPoints( std::vector<G4ThreeVector>* theNewVectorPointer );
    inline std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const;

  private:
    G4bool     isMomentumChanged;
    //  The flag which is set if momentum is changed in current step
    G4Material* theMaterialChange;
    const G4MaterialCutsCouple* theMaterialCutsCoupleChange;
    G4VSensitiveDetector* theSensitiveDetectorChange;
     // The material (and MaterialCutsCouple) where given track
     // currently locates

  private:
    std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer;
};

#include "G4ParticleChangeForTransport.icc"

#endif






