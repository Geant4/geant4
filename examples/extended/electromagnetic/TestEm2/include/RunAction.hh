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
// $Id: RunAction.hh,v 1.6 2004/09/17 10:51:38 maire Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

typedef  std::vector<G4double> MyVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;
class RunActionMessenger;

class G4Run;

namespace AIDA {
  class IAnalysisFactory;
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
    
     // Acceptance parameters
     void     SetEdepAndRMS(G4ThreeVector);
     
     G4double GetAverageEdep() const    {return edeptrue;};
     G4double GetRMSEdep() const        {return rmstrue;};
     G4double GetLimitEdep() const      {return limittrue;};

     // Histogram name and type
     void SetHistoName(G4String& val)   {histoName = val;};
     void SetHistoType(G4String& val)   {histoType = val;};
     
     const G4String& GetHistoName() const  {return histoName;};
     const G4String& GetHistoType() const  {return histoType;};
     
  private:

    void bookHisto();
    void cleanHisto();

  private:

    DetectorConstruction*   Det;
    PrimaryGeneratorAction* Kin;
    RunActionMessenger*     runMessenger;
    
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

    G4double edeptrue;
    G4double rmstrue;
    G4double limittrue;
    
    G4String histoName;
    G4String histoType;
    
    AIDA::IAnalysisFactory* af;
    AIDA::ITree*            tree;
    AIDA::IHistogram1D*     histo[12];
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
  dEdL[Lbin] += dEstep; dEdR[Rbin] += dEstep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

inline
void RunAction::particleFlux(G4ParticleDefinition* particle, G4int Lplan)
{
       if (particle == G4Gamma::Gamma())          gammaFlux[Lplan]++;
  else if (particle == G4Electron::Electron()) electronFlux[Lplan]++;
  else if (particle == G4Positron::Positron()) positronFlux[Lplan]++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

