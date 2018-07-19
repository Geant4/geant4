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
// $Id: G4AdjointPrimaryGeneratorAction.hh 98735 2016-08-09 10:54:06Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPrimaryGeneratorAction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	10-01-2007 creation by L. Desorgher 	
//		1-11-2009 Splitting of G4AdjointPrimaryGeneratorAction in two classes  G4AdjointPrimaryGeneratorAction and G4AdjointPrimaryGenerator L.Desorgher
//			   		
//
//-------------------------------------------------------------
//	Documentation:
//		This class represents the PrimaryGeneratorAction that is used during the entire adjoint simulation.
//		It uses the class G4AdjointPrimaryGenerator to generate randomly adjoint primary particles on a user selected 
//		adjoint source (External surface of a volume or Sphere).
//		The spectrum of the primary  adjoint particles is set as 1/E with user defined max and min energy.
//		The weight of the primary is set according to ReverseMC theory as  w=log(Emax/Emin)*E*adjoint_source_area*pi/n, with E the energy of the
//		particle, n the number of adjoint primary particles of same type that will be generated during the simulation.  
//		Different types of adjoint particles are generated event after event in order
//		to cover all the type of primaries and secondaries needed for the simulation. For example if reverse e- ionisation, brem, photo electric effect, and
//		compton are considered both adjoint gamma and adjoint e- will be considered alternatively as adjoint primary. 
//		The user can decide to consider/neglect some type of particle by using the macro
//		commands   /adjoint/ConsiderAsPrimary and /adjoint/NeglectAsPrimary. If  an adjoint primary or its secondary has reached the external surface, 
//		in the next event a fwd primary particle equivalent to the last generated adjoint primary is generated with the same position, energy but opposite direction
//		and the forward tracking phase starts. 
//		
//		
//
#ifndef G4AdjointPrimaryGeneratorAction_h
#define G4AdjointPrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include"G4ThreeVector.hh"
#include <vector>
#include <map>
#include <iterator>

class G4AdjointPosOnPhysVolGenerator;
class G4ParticleGun;
class G4Event;
class G4AdjointPrimaryGenerator;
class G4ParticleDefinition;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class G4AdjointPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public: //constructor, destructor

    G4AdjointPrimaryGeneratorAction();    
   ~G4AdjointPrimaryGeneratorAction();

  public: //public methods
    
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(const G4String& val) { rndmFlag = val;}
    void SetEmin(G4double val);
    void SetEmax(G4double val); 
    void SetEminIon(G4double val);
    void SetEmaxIon(G4double val);
    void SetSphericalAdjointPrimarySource(G4double radius, G4ThreeVector pos);
    void SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(const G4String& volume_name);
    void ConsiderParticleAsPrimary(const G4String& particle_name);
    void NeglectParticleAsPrimary(const G4String& particle_name);
    void SetPrimaryIon(G4ParticleDefinition* adjointIon, G4ParticleDefinition* fwdIon);
    void UpdateListOfPrimaryParticles();
    inline size_t GetNbOfAdjointPrimaryTypes(){return ListOfPrimaryAdjParticles.size();}
    inline std::vector<G4ParticleDefinition*>* GetListOfPrimaryFwdParticles(){
                                             return &ListOfPrimaryFwdParticles;}
    inline const G4String& GetPrimaryIonName(){return ion_name;}
    inline void SetNbPrimaryFwdGammasPerEvent(G4int nb) {nb_fwd_gammas_per_event=nb;}
    inline void SetNbAdjointPrimaryGammasPerEvent(G4int nb) {nb_adj_primary_gammas_per_event=nb;}
    inline void SetNbAdjointPrimaryElectronsPerEvent(G4int nb) {nb_adj_primary_electrons_per_event=nb;}
    inline  G4ParticleDefinition* GetLastGeneratedFwdPrimaryParticle(){return ListOfPrimaryFwdParticles[index_particle];}

  private: //private methods

    G4double ComputeEnergyDistWeight(G4double energy, G4double E1, G4double E2);
  
  private: //attributes
    
    G4String  rndmFlag;	  //flag for a rndm impact point
    
    //The generator of primary vertex except for weight
    G4AdjointPrimaryGenerator* theAdjointPrimaryGenerator;
  
    //Emin and Emax energies of the adjoint source
    //---------------------------------------------
    G4double Emin;
    G4double Emax;
    G4double EminIon;
    G4double EmaxIon;

    //List of type of primary adjoint and forward  particle used in the simulation
    //---------------------------------------------------------------------------
    std::vector<G4ParticleDefinition*> ListOfPrimaryFwdParticles;
    std::vector<G4ParticleDefinition*> ListOfPrimaryAdjParticles;
    std::map<G4String, G4bool> PrimariesConsideredInAdjointSim; //if true considered if false not considered


    size_t index_particle;

    G4ThreeVector  pos,  direction, p; 
   
    G4String type_of_adjoint_source; //Spherical ExtSurfaceOfAVolume
    G4double radius_spherical_source;
    G4ThreeVector center_spherical_source;
    G4int nb_fwd_gammas_per_event;
    G4int nb_adj_primary_gammas_per_event;
    G4int nb_adj_primary_electrons_per_event;
    
    //For simulation with ions
    //--------------------------
    G4ParticleDefinition* fwd_ion;
    G4ParticleDefinition* adj_ion;
    G4String ion_name;
    //disable copy constructor and assignement operator
    G4AdjointPrimaryGeneratorAction(const G4AdjointPrimaryGeneratorAction&);
    G4AdjointPrimaryGeneratorAction& operator=(const G4AdjointPrimaryGeneratorAction&);
};
#endif


