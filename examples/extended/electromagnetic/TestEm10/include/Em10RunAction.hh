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
// $Id: Em10RunAction.hh,v 1.5 2006-06-29 16:38:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10RunAction_h
#define Em10RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10RunMessenger;
class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10RunAction : public G4UserRunAction
{
  public:
    Em10RunAction();
   ~Em10RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void CountEvent();
    void CountParticles(G4double,G4double);
    void AddEP(G4double,G4double);
    void AddEdeps(G4double Eabs); 
    void AddTrackLength(G4double tlabs); 
    void AddnStepsCharged(G4double ns);
    void AddnStepsNeutral(G4double ns);

    void AddTrRef(G4double tr,G4double ref);

    void FillEn(G4double En);
    void FillTh(G4double Th);
    void FillThBack(G4double Th);
    void FillR(G4double R);
    void FillTt(G4double Tt);
    void FillTb(G4double Tt);
    void FillTsec(G4double T);
    void FillGammaSpectrum(G4double E);
    void FillNbOfSteps(G4double nstep);
    void Fillvertexz(G4double z);

    void SethistName(G4String name) ;

    void SetnbinStep(G4int nbin);
    void SetSteplow(G4double Slow);
    void SetStephigh(G4double Shigh);

    void SetnbinEn(G4int nbin);
    void SetEnlow(G4double Elow);
    void SetEnhigh(G4double Enhigh);

    void SetnbinTt(G4int nbin);
    void SetTtlow(G4double Ttlow);
    void SetTthigh(G4double Tthigh);

    void SetnbinTb(G4int nbin);
    void SetTblow(G4double Tblow);
    void SetTbhigh(G4double Tbhigh);

    void SetnbinTsec(G4int nbin);
    void SetTseclow(G4double Tlow);
    void SetTsechigh(G4double Thigh);

    void SetnbinTh(G4int nbin);
    void SetThlow(G4double Thlow);
    void SetThhigh(G4double Thhigh);

    void SetnbinThBack(G4int nbin);
    void SetThlowBack(G4double Thlow);
    void SetThhighBack(G4double Thhigh);

    void SetnbinR(G4int nbin);
    void SetRlow(G4double Rlow);
    void SetRhigh(G4double Rhigh);

    void SetnbinGamma(G4int nbin);
    void SetElowGamma(G4double Elow);
    void SetEhighGamma(G4double Ehigh);

    void Setnbinzvertex(G4int nbin);
    void Setzlow(G4double z);
    void Setzhigh(G4double z);
    
    void  SetRndmFreq(G4int val) {saveRndm = val;}
    G4int GetRndmFreq()          {return saveRndm;}
    
  private:

    void bookHisto();

  private:

    G4String histName ;

    G4double EnergySumAbs,EnergySquareSumAbs;
    G4double tlSumAbs,tlsquareSumAbs;
    G4double nStepSumCharged,nStepSum2Charged ;
    G4double nStepSumNeutral,nStepSum2Neutral ;
    G4double TotNbofEvents;
    G4double SumCharged,Sum2Charged,SumNeutral,Sum2Neutral;
    G4double Selectron,Spositron;

    G4double Transmitted,Reflected ;

    G4double entryStep,underStep,overStep,distStep[200];
    G4double Steplow,Stephigh,dStep;
    G4int    nbinStep;
    G4double entryEn,underEn,overEn,distEn[200];
    G4double Enlow,Enhigh,dEn;
    G4int    nbinEn;
    G4double entryTt,underTt,overTt,distTt[200];
    G4double Ttlow,Tthigh,dTt;
    G4int    nbinTt;
    G4double Ttmean,Tt2mean;
    G4double entryTb,underTb,overTb,distTb[200];
    G4double Tblow,Tbhigh,dTb;
    G4int    nbinTb;
    G4double Tbmean,Tb2mean;
    G4double entryTsec,underTsec,overTsec,distTsec[200];
    G4double Tseclow,Tsechigh,dTsec;
    G4int    nbinTsec;
    G4double entryTh,underTh,overTh,distTh[200];
    G4double Thlow,Thhigh,dTh;
    G4int    nbinTh;
    G4double entryThback,underThback,overThback,distThback[200];
    G4double Thlowback,Thhighback,dThback;
    G4int    nbinThback;
    G4double entryR ,underR ,overR ,distR[200];
    G4double Rlow,Rhigh,dR;
    G4int    nbinR;
    G4double Rmean,R2mean;
    G4double entryGamma,underGamma,overGamma,distGamma[200];
    G4double ElowGamma,EhighGamma,dEGamma;
    G4int nbinGamma ;
    G4double entryvertexz,undervertexz,oververtexz,distvertexz[200];
    G4double zlow,zhigh,dz;
    G4int nbinvertexz;
 
    Em10RunMessenger* runMessenger;
    G4int saveRndm;    
};

#endif

