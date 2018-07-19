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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $ID$
/// \file Run.hh
/// \brief Definition of the Run class

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>
#include "G4Molecule.hh"

class DetectorConstruction;
class G4ParticleDefinition;
//class NeuronLoadDataFile;
class G4Molecule;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);   
    void AddPrimaryLET(G4double let);   
    void SetTrackLength(G4double tracklen);  
    void CountProcesses(const G4VProcess* process);
    void ParticleCount(G4String, G4double);
    void ParticleCountNeuron(G4String, G4double); 
    void MoleculeCount(G4String, G4double);
    void MoleculeCountNeuron(G4Molecule* molecule); 
    void AddEdep (G4double edep);
    void AddEflow (G4double eflow);                   
    void ParticleFlux(G4String, G4double); 
                          
    virtual void Merge(const G4Run*);
    void EndOfRun();     

// Edep in all volume
  G4double GetEdepALL(){return fEdepAll;}
  void SetEdepALL(G4double vall){ fEdepAll = vall;}
  void AddEdepALL (G4double vall)
  { 
   fEdepAll += vall; 
   fEdepAll_err += vall*vall;
  }   
// 0. Edep in homogeneous Medium
  G4double GetEdepMedium(){return fEdepMedium;}
  void SetEdepMedium(G4double vall){ fEdepMedium = vall;}
  void AddEdepMedium (G4double vall)
  { 
   fEdepMedium += vall; 
   fEdepMedium_err += vall*vall;
  }  
// 1. Edep in Bounding Slice Volume
  G4double GetEdepSlice(){return fEdepSlice;}
  void SetEdepSlice(G4double vall){ fEdepSlice = vall;}
  void AddEdepSlice (G4double vall)
  { 
   fEdepSlice += vall; 
   fEdepSlice_err += vall*vall;
  }
// 2. Edep in Soma volume
  G4double GetEdepSoma(){return fEdepSoma;}
  void SetEdepSoma(G4double vall){ fEdepSoma = vall;}
  void AddEdepSoma (G4double vall)
  { 
   fEdepSoma += vall; 
   fEdepSoma_err += vall*vall;
  }
  
  void AddSomaCompart(G4int i, G4double x){ fSoma3DEdep[i] +=x;}
  G4double GetSomaCompart(G4int i){ return fSoma3DEdep[i];}  
  //G4double GetSomalength(G4int i){ return fSoma3DEdep[i];}  
  //G4ThreeVector GetVectSoma(G4int i) {return fSomaCompartments[i];}
  
// 3. Edep in Dendrites volume
  G4double GetEdepDend(){return fEdepDend;}
  void SetEdepDend(G4double vall){ fEdepDend = vall;}
  void AddEdepDend (G4double vall)
  { 
   fEdepDend += vall; 
   fEdepDend_err += vall*vall;
  }
  
  void AddDendCompart(G4int i, G4double x){ fDend3DEdep[i] +=x;}
  G4double GetDendCompart(G4int i){ return fDend3DEdep[i];}  
  
// 4. Edep in Axon volume
  G4double GetEdepAxon(){return fEdepAxon;}
  void SetEdepAxon(G4double vall){ fEdepAxon = vall;}
  void AddEdepAxon (G4double vall)
  { 
   fEdepAxon += vall; 
   fEdepAxon_err += vall*vall;
  }
  
  void AddAxonCompart(G4int i, G4double x){ fAxon3DEdep[i] +=x;}
  G4double GetAxonCompart(G4int i){ return fAxon3DEdep[i];}  
  
// 5. Edep in whole Neuron volume
  G4double GetEdepNeuron(){return fEdepNeuron;}
  void SetEdepNeuron(G4double vall){ fEdepNeuron = vall;}
  void AddEdepNeuron (G4double vall)
  { 
   fEdepNeuron += vall; 
   fEdepNeuron_err += vall*vall;
  }  

  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };

    //NeuronLoadDataFile * fNeuronLoadParamz;  
    DetectorConstruction* fDetector;
    G4ParticleDefinition* fParticle;
    G4double              fEkin;
    G4double     fLET, fLET2;
    G4double     ftrackLength;
    G4double   fTrackLen,  fTrackLen2;
    G4int fNParticle;
 
    G4double * fSoma3DEdep;
    G4double * fDend3DEdep ;
    G4double * fAxon3DEdep ;  
    G4double * fNeuron3DEdep;  
    G4double fEdepAll,  fEdepAll_err,fEdepMedium, fEdepMedium_err, fEdepSlice,
    fEdepSlice_err,fEdepSoma,fEdepSoma_err, fEdepDend,  fEdepDend_err,fEdepAxon,
    fEdepAxon_err,fEdepNeuron,  fEdepNeuron_err;
  
    G4int fmolNum, fmolNum2;
    G4double fEnergyFlow,    fEnergyFlow2;  
    std::map<G4String,G4int>   fMoleculeNumber; 
    std::map<G4String,G4int>        fProcCounter;
    std::map<G4String,ParticleData> fParticleDataMap1;                    
    std::map<G4String,ParticleData> fParticleDataMap2;
};

#endif
