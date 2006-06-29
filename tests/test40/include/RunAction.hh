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
// $Id: RunAction.hh,v 1.2 2006-06-29 22:03:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 08.03.01 Hisaya: Adapted MyVector for STL   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4DataVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;

class G4Run;

namespace AIDA {
  class ITree;
  class IHistogram1D;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:

    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    inline void initializePerEvent();
           void fillPerEvent();
    inline void fillPerTrack(G4double,G4double);
    inline void fillPerStep (G4double,G4int,G4int);
    inline void particleFlux(G4ParticleDefinition*,G4int);

  private:

    void bookHisto();
    void saveHisto();

  private:

    DetectorConstruction*   Det;
    PrimaryGeneratorAction* Kin;

    G4int nLbin;
    G4DataVector dEdL;
    G4DataVector sumELongit;
    G4DataVector sumE2Longit;
    G4DataVector sumELongitCumul;
    G4DataVector sumE2LongitCumul;

    G4int nRbin;
    G4DataVector dEdR;
    G4DataVector sumERadial;
    G4DataVector sumE2Radial;
    G4DataVector sumERadialCumul;
    G4DataVector sumE2RadialCumul;

    G4DataVector gammaFlux;
    G4DataVector electronFlux;
    G4DataVector positronFlux;

    G4double ChargTrLength;
    G4double sumChargTrLength;
    G4double sum2ChargTrLength;

    G4double NeutrTrLength;
    G4double sumNeutrTrLength;
    G4double sum2NeutrTrLength;

    AIDA::ITree* tree;             // the tree should only be deleted at the end
    AIDA::IHistogram1D* histo[12];   // (after writing the histos to file)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void RunAction::initializePerEvent()
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
void RunAction::fillPerTrack(G4double charge, G4double trkLength)
{
  if (charge != 0.) ChargTrLength += trkLength;
  else              NeutrTrLength += trkLength;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void RunAction::fillPerStep(G4double dEstep, G4int Lbin, G4int Rbin)
{
  if(Lbin<nLbin && Lbin>=0) dEdL[Lbin] += dEstep; 
  if(Rbin<nRbin && Rbin>=0) dEdR[Rbin] += dEstep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

inline
void RunAction::particleFlux(G4ParticleDefinition* particle, G4int Lplan)
{
  if(Lplan<nLbin && Lplan>=0) { 
         if (particle == G4Gamma::Gamma())          gammaFlux[Lplan]+=1.0;
    else if (particle == G4Electron::Electron()) electronFlux[Lplan]+=1.0;
    else if (particle == G4Positron::Positron()) positronFlux[Lplan]+=1.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

