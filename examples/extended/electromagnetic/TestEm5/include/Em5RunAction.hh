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
// $Id: Em5RunAction.hh,v 1.9 2002-06-05 15:43:42 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em5RunAction_h
#define Em5RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "g4std/iostream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em5RunMessenger;
class G4Run;

#ifndef G4NOHIST
 class HepTupleManager;
 class HepHistogram;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em5RunAction : public G4UserRunAction
{
  public:
    Em5RunAction();
   ~Em5RunAction();

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

    G4double GetdTh() { return dTh;};
    G4double GetdThback() { return dThback;};
#ifndef G4NOHIST
    HepHistogram* GetHisto(G4int id) {return histo[id];}
#endif

  private:

    void bookHisto();

  private:

    G4String histName ;
    
#ifndef G4NOHIST    
    HepTupleManager* hbookManager;
    HepHistogram* histo[10] ;    
#endif
    
    G4double EnergySumAbs,EnergySquareSumAbs;
    G4double tlSumAbs,tlsquareSumAbs;
    G4double nStepSumCharged,nStepSum2Charged ;
    G4double nStepSumNeutral,nStepSum2Neutral ;
    G4double TotNbofEvents;
    G4double SumCharged,Sum2Charged,SumNeutral,Sum2Neutral;
    G4double Selectron,Spositron;

    G4double Transmitted,Reflected ;

    G4double Steplow,Stephigh;
    G4int    nbinStep;
    G4double Enlow,Enhigh;
    G4int    nbinEn;
    G4double Ttlow,Tthigh;
    G4int    nbinTt;
    G4double Tblow,Tbhigh;
    G4int    nbinTb;
    G4double Tseclow,Tsechigh;
    G4int    nbinTsec;
    G4double Thlow,Thhigh,dTh;
    G4int    nbinTh;
    G4double Thlowback,Thhighback,dThback;
    G4int    nbinThback;
    G4double Rlow,Rhigh;
    G4int    nbinR;
    G4double ElowGamma,EhighGamma;
    G4int nbinGamma ;
    G4double zlow,zhigh;
    G4int nbinvertexz;
 
    Em5RunMessenger* runMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

