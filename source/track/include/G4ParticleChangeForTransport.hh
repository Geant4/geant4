//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParticleChangeForTransport.hh,v 1.9 2002-11-01 15:55:30 jacek Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


class G4ParticleChangeForTransport: public G4ParticleChange
{ 
  public:
    // default constructor
    G4ParticleChangeForTransport();

    // destructor
    virtual ~G4ParticleChangeForTransport();

  protected:
    // hide copy constructor and assignment operaor as protected
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

    G4Material* GetMaterialChange() const;
    void SetMaterialChange(G4Material* fMaterial);
    //  Get/Set the material in the touchable of the current particle.

    G4bool GetMomentumChanged() const;
    void SetMomentumChanged(G4bool b);

  public:
    virtual void DumpInfo() const;

  protected:
    G4TouchableHandle theTouchableHandle;
    //  The changed touchable of a given particle.

  private:
    G4bool     isMomentumChanged;
    //  The flag which is set if mometum is changed in this stepi
    G4Material* theMaterialChange;
     //  The material where given track currently locates

  // Prototyping implementation of smooth representation of curved
  // trajectories. (jacek 30/10/2002)
public:
  // Auxiliary points are ThreeVectors for now; change to
  // G4AuxiliaryPoints or some such (jacek 30/10/2002)
  void SetPointerToVectorOfAuxiliaryPoints( G4std::vector<G4ThreeVector>* theNewVectorPointer ) {
    fpVectorOfAuxiliaryPointsPointer = theNewVectorPointer;
  }
  G4std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const {
    return fpVectorOfAuxiliaryPointsPointer;
  }
private:
  // Explicity including the word "Pointer" in the name as I keep
  // forgetting the * (jacek 30/10/2002)
  G4std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer;

};

#include "G4ParticleChangeForTransport.icc"

#endif
