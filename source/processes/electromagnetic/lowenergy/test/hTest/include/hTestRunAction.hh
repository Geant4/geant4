// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The user has a possibility to define and to fill his histograms in this class.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef hTestRunAction_h
#define hTestRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4IonC12.hh"
#include "G4IonAr40.hh"
#include "G4IonFe56.hh"
#include "G4Electron.hh"
#include "g4std/iostream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestRunMessenger;

#ifndef G4NOHIST
class HepTupleManager;
class HepTuple;
class HepHistogram;
#endif

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestRunAction : public G4UserRunAction
{
public: // Without description

    hTestRunAction();
   ~hTestRunAction();

public: // With description
 
    void BeginOfRunAction(const G4Run*);
  // In this method histogramms are booked

    void   EndOfRunAction(const G4Run*);
  // In this method bookHisto method is called in which histogramms are filled

public: // Without description

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

#ifndef G4NOHIST
    void SethistName(G4String name) ;
#endif

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
#ifndef G4NOHIST
    void bookHisto();
    HepTuple* GetNtuple();
#endif
    void SaveToTuple(G4String parname,G4double val);
    void SaveToTuple(G4String parname,G4double val,G4double defval);
    void SaveEvent();

  private:

#ifndef G4NOHIST
    void FillLowEnergyTest();
#endif

  private:

#ifndef G4NOHIST
    G4String histName ;
    HepTupleManager* hbookManager;
    HepHistogram *histo1, *histo2, *histo3, *histo4, *histo5 ;
    HepHistogram *histo6, *histo7, *histo8, *histo9, *histo10;
    HepHistogram *histo11, *histo12, *histo13, *histo14, *histo15;
    HepHistogram *histo21, *histo22, *histo23, *histo24, *histo25;
    HepHistogram *histo31, *histo32, *histo33, *histo34, *histo35;
    HepHistogram *histo41, *histo42, *histo43, *histo44, *histo45;
    HepHistogram *histo51, *histo52, *histo53, *histo54, *histo55;
    HepHistogram *histo61, *histo62, *histo63, *histo64, *histo65;
    HepHistogram *histo71, *histo72, *histo73, *histo74, *histo75, *histo76;
    HepHistogram *histo81, *histo82, *histo83, *histo84, *histo85, *histo86;
    HepTuple *ntup;
#endif

    G4double LowestEnergy,HighestEnergy;
    G4int TotBin;

    G4ParticleDefinition* theProton ;
    const G4Electron* theElectron ;

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
 
    hTestRunMessenger* runMessenger ;

};

#endif

