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
/// \file electromagnetic/TestEm2/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "g4root.hh"
////#include "g4xml.hh"
////#include "g4hbook.hh"

#include <vector>
typedef  std::vector<G4double> MyVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;
class RunActionMessenger;

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:

    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

           void InitializePerEvent();
           void FillPerEvent();
    inline void FillPerTrack(G4double,G4double);
    inline void FillPerStep (G4double,G4int,G4int);
    inline void AddStep (G4double q);
    
    void SetVerbose(G4int val)  {fVerbose = val;};
    
     // Acceptance parameters
     void     SetEdepAndRMS(G4ThreeVector);
     
     G4double GetAverageEdep() const    {return fEdeptrue;};
     G4double GetRMSEdep() const        {return fRmstrue;};
     G4double GetLimitEdep() const      {return fLimittrue;};

     // Histogram name and type
     void SetHistoName(G4String& val)   {fHistoName[0] = val;};
     
     const G4String& GetHistoName() const  {return fHistoName[1];};
     
  private:

    void BookHisto();
    void CleanHisto();

  private:

    DetectorConstruction*   fDet;
    PrimaryGeneratorAction* fKin;
    RunActionMessenger*     fRunMessenger;
    
    G4int f_nLbin;
    MyVector f_dEdL;
    MyVector fSumELongit;
    MyVector fSumE2Longit;
    MyVector fSumELongitCumul;
    MyVector fSumE2LongitCumul;

    G4int f_nRbin;
    MyVector f_dEdR;
    MyVector fSumERadial;
    MyVector fSumE2Radial;
    MyVector fSumERadialCumul;
    MyVector fSumE2RadialCumul;

    G4double fChargTrLength;
    G4double fSumChargTrLength;
    G4double fSum2ChargTrLength;

    G4double fNeutrTrLength;
    G4double fSumNeutrTrLength;
    G4double fSum2NeutrTrLength;

    G4double fEdeptrue;
    G4double fRmstrue;
    G4double fLimittrue;

    G4double fChargedStep;
    G4double fNeutralStep;    
    
    G4int    fVerbose;
    
    G4String fHistoName[2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void RunAction::FillPerTrack(G4double charge, G4double trkLength)
{
  if (charge != 0.) fChargTrLength += trkLength;
  else              fNeutrTrLength += trkLength;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void RunAction::FillPerStep(G4double dEstep, G4int Lbin, G4int Rbin)
{
  f_dEdL[Lbin] += dEstep; f_dEdR[Rbin] += dEstep;
}

inline void RunAction::AddStep(G4double q)
{
  if (q == 0.0) { fNeutralStep += 1.0; }
  else          { fChargedStep += 1.0; }  
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

