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
// $Id: Tst26RunAction.hh,v 1.2 2003-02-01 18:14:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst26RunAction_h
#define Tst26RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"

#include "g4std/vector"

typedef  G4std::vector<G4double> MyVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26DetectorConstruction;
class Tst26PrimaryGeneratorAction;

class G4Run;

#ifndef G4NOHIST
namespace AIDA {
 class ITree;
 class IHistogram1D;
} 
#endif

class Tst26RunAction : public G4UserRunAction
{
  public:
  
    Tst26RunAction(Tst26DetectorConstruction*, Tst26PrimaryGeneratorAction*);
   ~Tst26RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    inline void initializePerEvent();
           void fillPerEvent();
    inline void fillPerTrack(G4double,G4double);
    inline void fillPerStep (G4double,G4int,G4int);
    inline void particleFlux(G4ParticleDefinition*,G4int);
 
  private:
  
    void bookHisto();
    void cleanHisto();
    
  private:
    
    Tst26DetectorConstruction*   Tst26Det;
    Tst26PrimaryGeneratorAction* Tst26Kin;
    
    G4int nLbin;    
    MyVector dEdL;
    MyVector sumELongit;
    MyVector sumE2Longit;
    MyVector sumELongitCumul;
    MyVector sumE2LongitCumul;
    
    G4int nRbin;    
    MyVector dEdR;
    MyVector sumERadial;
    MyVector sumE2Radial;
    MyVector sumERadialCumul;
    MyVector sumE2RadialCumul;
        
    MyVector gammaFlux;
    MyVector electronFlux;
    MyVector positronFlux;
    
    G4double ChargTrLength;
    G4double sumChargTrLength;
    G4double sum2ChargTrLength;
    
    G4double NeutrTrLength;
    G4double sumNeutrTrLength;
    G4double sum2NeutrTrLength;
                  
#ifndef G4NOHIST
    AIDA::ITree* tree;             // the tree should only be deleted at the end
    AIDA::IHistogram1D* histo[12];   // (after writing the histos to file)
#endif    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Tst26RunAction::initializePerEvent()
{
  //initialize arrays of energy deposit per bin     
  for (G4int i=0; i<nLbin; i++)
     { dEdL[i] = 0.; }
     
  for (G4int j=0; j<nRbin; j++)
     { dEdR[j] = 0.; }     
  
  //initialize tracklength 
    ChargTrLength = NeutrTrLength = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Tst26RunAction::fillPerTrack(G4double charge, G4double trkLength)
{
  if (charge != 0.) ChargTrLength += trkLength;
  else              NeutrTrLength += trkLength;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Tst26RunAction::fillPerStep(G4double dEstep, G4int Lbin, G4int Rbin)
{
  dEdL[Lbin] += dEstep; dEdR[Rbin] += dEstep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Tst26RunAction::particleFlux(G4ParticleDefinition* particle, G4int Lplan)
{
       if (particle == G4Gamma::Gamma())          gammaFlux[Lplan]++;
  else if (particle == G4Electron::Electron()) electronFlux[Lplan]++;
  else if (particle == G4Positron::Positron()) positronFlux[Lplan]++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

