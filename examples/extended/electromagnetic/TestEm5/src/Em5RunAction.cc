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
// $Id: Em5RunAction.cc,v 1.14 2002-08-08 13:52:19 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Em5RunAction.hh"
#include "Em5RunMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "Randomize.hh"

#ifndef G4NOHIST
 #include "AIDA/IAnalysisFactory.h"
 #include "AIDA/ITreeFactory.h"
 #include "AIDA/ITree.h"
 #include "AIDA/IHistogramFactory.h"
 #include "AIDA/IHistogram1D.h"
 #include "AIDA/IAxis.h"
 #include "AIDA/IAnnotation.h"
 #include "AIDA/ITupleFactory.h"
 #include "AIDA/ITuple.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5RunAction::Em5RunAction()
  :histName("histfile"),nbinStep(0),nbinEn(0),nbinTt(0),nbinTb(0),
   nbinTsec(0),nbinTh(0),nbinThback(0),nbinR(0),nbinGamma(0),nbinvertexz(0)
{
  runMessenger = new Em5RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5RunAction::~Em5RunAction()
{
  delete runMessenger;
#ifndef G4NOHIST   
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
 
  delete tree;
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::bookHisto()
{
#ifndef G4NOHIST
 // Creating the analysis factory
 IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 tree = tf->create(histName, false, false, "hbook");

 // Creating a histogram factory, whose histograms will be handled by the tree
 IHistogramFactory* hf = af->createHistogramFactory(*tree);

 // creating histograms
 //
  if(nbinStep>0)
  {
    histo[0] = hf->create1D("1","number of steps/event"
                               ,nbinStep,Steplow,Stephigh);
  }
  if(nbinEn>0)
  {
    histo[1] = hf->create1D("2","energy deposit in absorber(in MeV)"
                               ,nbinEn,Enlow,Enhigh);
  }
  if(nbinTh>0)
  {
    histo[2] = hf->create1D("3","angle distribution at exit(deg)"
                               ,nbinTh,Thlow/deg,Thhigh/deg);
  }
  if(nbinR>0)
  {
    histo[3] = hf->create1D("4","lateral distribution at exit(mm)"
                               ,nbinR ,Rlow,Rhigh);
  }
  if(nbinTt>0)
  {
    histo[4] = hf->create1D("5","kinetic energy of the primary at exit(MeV)"
                               ,nbinTt,Ttlow,Tthigh);
  }
  if(nbinThback>0)
  {
    histo[5] = hf->create1D("6",
                           "angle distribution of backscattered primaries(deg)"
                           ,nbinThback,Thlowback/deg,Thhighback/deg);
  }
  if(nbinTb>0)
  {
    histo[6] = hf->create1D("7",
                           "kinetic energy of the backscattered primaries (MeV)"
                           ,nbinTb,Tblow,Tbhigh);
  }
  if(nbinTsec>0)
  {
    histo[7] = hf->create1D("8","kinetic energy of the charged secondaries(MeV)"
                               ,nbinTsec,Tseclow,Tsechigh);
  }
  if(nbinvertexz>0)
  {
    histo[8] = hf->create1D("9","x of secondary charged vertices(mm)"
                               ,nbinvertexz ,zlow,zhigh);
  }
  if(nbinGamma>0)
  {
    histo[9]= hf->create1D("10",
                      "kinetic energy of gammas escaping the absorber(MeV)"
                                //     ,nbinGamma,ElowGamma,EhighGamma);
                          ,nbinGamma,log10(ElowGamma),log10(EhighGamma));
  }
  
 delete hf;
 delete tf;
 delete af;  
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
      
  EnergySumAbs = 0. ;
  EnergySquareSumAbs = 0.;
  tlSumAbs = 0. ;
  tlsquareSumAbs = 0. ;
  nStepSumCharged = 0. ;
  nStepSum2Charged= 0. ;
  nStepSumNeutral = 0. ;
  nStepSum2Neutral= 0. ;
  TotNbofEvents = 0. ;
  SumCharged=0.;
  SumNeutral=0.;
  Sum2Charged=0.;
  Sum2Neutral=0.;
  Selectron=0.;
  Spositron=0.;

  Transmitted=0.;
  Reflected  =0.;
//  plot definitions     --------------------------------
  if(nbinTh>0)
    dTh = (Thhigh-Thlow)/nbinTh ;

  if(nbinThback>0)
    dThback = (Thhighback-Thlowback)/nbinThback ;

  if (aRun->GetRunID() == 0) bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double sAbs,sigAbs,sigstep,sigcharged,signeutral;

  tlSumAbs /= TotNbofEvents ;
  sAbs = tlsquareSumAbs/TotNbofEvents-tlSumAbs*tlSumAbs ;
  if(sAbs>0.)
    sAbs = sqrt(sAbs/TotNbofEvents) ;
  else
    sAbs = 0. ;
  
  EnergySumAbs /= TotNbofEvents ;
  sigAbs = EnergySquareSumAbs/TotNbofEvents-EnergySumAbs*EnergySumAbs;
  if(sigAbs>0.)
    sigAbs = sqrt(sigAbs/TotNbofEvents);
  else
    sigAbs = 0.;

  nStepSumCharged /= TotNbofEvents ;
  sigstep = nStepSum2Charged/TotNbofEvents-nStepSumCharged*nStepSumCharged;
  if(sigstep>0.)
    sigstep = sqrt(sigstep/TotNbofEvents);
  else
    sigstep = 0.;
  G4double sigch=sigstep ;
  
  nStepSumNeutral /= TotNbofEvents ;
  sigstep = nStepSum2Neutral/TotNbofEvents-nStepSumNeutral*nStepSumNeutral;
  if(sigstep>0.)
    sigstep = sqrt(sigstep/TotNbofEvents);
  else
    sigstep = 0.;
  G4double signe=sigstep ;
  
  SumCharged /= TotNbofEvents;
  sigcharged = Sum2Charged/TotNbofEvents-SumCharged*SumCharged; 
  if(sigcharged>0.)
    sigcharged = sqrt(sigcharged/TotNbofEvents);
  else
    sigcharged = 0. ;
 
  SumNeutral /= TotNbofEvents;
  signeutral = Sum2Neutral/TotNbofEvents-SumNeutral*SumNeutral; 
  if(signeutral>0.)
    signeutral = sqrt(signeutral/TotNbofEvents);
  else
    signeutral = 0. ;
 
  Selectron /= TotNbofEvents ;
  Spositron /= TotNbofEvents ;

  Transmitted /=TotNbofEvents ;
  Reflected   /=TotNbofEvents ;
 G4cout << " ================== run summary =====================" << G4endl;
 G4int prec = G4cout.precision(6);
  G4cout << " end of Run TotNbofEvents = " <<  
           TotNbofEvents << G4endl ;
  G4cout << "    mean charged track length   in absorber=" <<
           tlSumAbs/mm      << " +- " << sAbs/mm    <<
          "  mm  " << G4endl; 
  G4cout << G4endl;
  G4cout << "            mean energy deposit in absorber=" <<
           EnergySumAbs/MeV << " +- " << sigAbs/MeV <<
          "  MeV " << G4endl ;
  G4cout << G4endl ;
  G4cout << " mean number of steps in absorber (charged) =" <<
           nStepSumCharged         << " +- " << sigch     <<
          "      " << G4endl ;
  G4cout << " mean number of steps in absorber (neutral) =" <<
           nStepSumNeutral         << " +- " << signe     <<
          "      " << G4endl ;
  G4cout << G4endl ;
  G4cout << "   mean number of charged secondaries = " <<
           SumCharged << " +- " << sigcharged << G4endl;  
  G4cout << G4endl ;
  G4cout << "   mean number of neutral secondaries = " <<
           SumNeutral << " +- " << signeutral << G4endl;  
  G4cout << G4endl ;
  
  G4cout << "   mean number of e-s =" << Selectron << 
            "  and e+s =" << Spositron << G4endl;
  G4cout << G4endl; 
  
  G4cout << "(number) transmission coeff=" << Transmitted <<
            "  reflection coeff=" << Reflected << G4endl;
  G4cout << G4endl; 

 G4cout.precision(prec);
  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
   
  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::CountEvent()
{
  TotNbofEvents += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::AddnStepsCharged(G4double ns)
{
  nStepSumCharged += ns;
  nStepSum2Charged += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::AddnStepsNeutral(G4double ns)
{
  nStepSumNeutral += ns;
  nStepSum2Neutral += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::AddEdeps(G4double Eabs)
{
  EnergySumAbs += Eabs;
  EnergySquareSumAbs += Eabs*Eabs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::AddTrackLength(G4double tlabs)
{
  tlSumAbs += tlabs;
  tlsquareSumAbs += tlabs*tlabs ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::AddTrRef(G4double tr,G4double ref)
{
  Transmitted += tr ;
  Reflected   += ref;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5RunAction::SethistName(G4String name)
{
  histName = name ;
  G4cout << " hist file = " << histName << G4endl;
}

void Em5RunAction::SetnbinStep(G4int nbin)
{
  nbinStep = nbin ;
  if(nbinStep>0)
  G4cout << " Nb of bins in #step plot = " << nbinStep << G4endl ;
}
void Em5RunAction::SetSteplow(G4double low)
{
  Steplow = low ;
  if(nbinStep>0)
  G4cout << " low  in the #step plot = " << Steplow << G4endl ;
}
void Em5RunAction::SetStephigh(G4double high)
{
  Stephigh = high ;
  if(nbinStep>0)
  G4cout << " high in the #step plot = " << Stephigh << G4endl ;
}
void Em5RunAction::SetnbinEn(G4int nbin)
{
  nbinEn = nbin ;
  if(nbinEn>0)
  G4cout << " Nb of bins in Edep plot = " << nbinEn << G4endl ;
}
void Em5RunAction::SetEnlow(G4double Elow)
{
  Enlow = Elow ;
  if(nbinEn>0)
  G4cout << " Elow  in the  Edep plot = " << Enlow << G4endl ;
}
void Em5RunAction::SetEnhigh(G4double Ehigh)
{
  Enhigh = Ehigh ;
  if(nbinEn>0)
  G4cout << " Ehigh in the  Edep plot = " << Enhigh << G4endl ;
}

void Em5RunAction::SetnbinGamma(G4int nbin)
{
  nbinGamma = nbin ;
  if(nbinGamma>0)
  G4cout << " Nb of bins in gamma spectrum plot = " << nbinGamma << G4endl ;
}
void Em5RunAction::SetElowGamma(G4double Elow)
{
  ElowGamma = Elow ;
  if(nbinGamma>0)
  G4cout << " Elow  in the gamma spectrum plot = " << ElowGamma << G4endl ;
}
void Em5RunAction::SetEhighGamma(G4double Ehigh)
{
  EhighGamma = Ehigh ;
  if(nbinGamma>0)
  G4cout << " Ehigh in the gamma spectrum plot = " << EhighGamma << G4endl ;
}

void Em5RunAction::SetnbinTt(G4int nbin)
{
  nbinTt = nbin ;
  if(nbinTt>0)
  G4cout << " Nb of bins in Etransmisssion plot = " << nbinTt << G4endl ;
}
void Em5RunAction::SetTtlow(G4double Elow)
{
  Ttlow = Elow ;
  if(nbinTt>0)
  G4cout << " Elow  in the  Etransmission plot = " << Ttlow << G4endl ;
}
void Em5RunAction::SetTthigh(G4double Ehigh)
{
  Tthigh = Ehigh ;
  if(nbinTt>0)
  G4cout << " Ehigh in the  Etransmission plot = " << Tthigh << G4endl ;
}
void Em5RunAction::SetnbinTb(G4int nbin)
{
  nbinTb = nbin ;
  if(nbinTb>0)
  G4cout << " Nb of bins in Ebackscattered plot = " << nbinTb << G4endl ;
}
void Em5RunAction::SetTblow(G4double Elow)
{
  Tblow = Elow ;
  if(nbinTb>0)
  G4cout << " Elow  in the  Ebackscattered plot = " << Tblow << G4endl ;
}
void Em5RunAction::SetTbhigh(G4double Ehigh)
{
  Tbhigh = Ehigh ;
  if(nbinTb>0)
  G4cout << " Ehigh in the  Ebackscattered plot = " << Tbhigh << G4endl ;
}

void Em5RunAction::SetnbinTsec(G4int nbin)
{
  nbinTsec = nbin ;
  if(nbinTsec>0)
  G4cout << " Nb of bins in Tsecondary  plot = " << nbinTsec << G4endl ;
}
void Em5RunAction::SetTseclow(G4double Elow)
{
  Tseclow = Elow ;
  if(nbinTsec>0)
  G4cout << " Elow  in the  Tsecondary plot = " << Tseclow << G4endl ;
}
void Em5RunAction::SetTsechigh(G4double Ehigh)
{
  Tsechigh = Ehigh ;
  if(nbinTsec>0)
  G4cout << " Ehigh in the  Tsecondary plot = " << Tsechigh << G4endl ;
}
 
void Em5RunAction::SetnbinR(G4int nbin)
{
  nbinR  = nbin ;
  if(nbinR>0)
  G4cout << " Nb of bins in R plot = " << nbinR << G4endl ;
}
void Em5RunAction::SetRlow(G4double rlow)
{
  Rlow = rlow ;
  if(nbinR>0)
  G4cout << " Rlow  in the  R plot = " << Rlow << G4endl ;
}
void Em5RunAction::SetRhigh(G4double rhigh)
{
  Rhigh = rhigh ;
  if(nbinR>0)
  G4cout << " Rhigh in the R plot = " << Rhigh << G4endl ;
}

void Em5RunAction::Setnbinzvertex(G4int nbin)
{
  nbinvertexz  = nbin ;
  if(nbinvertexz>0)
  G4cout << " Nb of bins in Z plot = " << nbinvertexz << G4endl ;
}
void Em5RunAction::Setzlow(G4double z)
{
  zlow = z ;
  if(nbinvertexz>0)
  G4cout << " zlow  in the  Z plot = " << zlow << G4endl ;
}
void Em5RunAction::Setzhigh(G4double z)
{
  zhigh = z ;
  if(nbinvertexz>0)
  G4cout << " zhigh in the Z plot = " << zhigh << G4endl ;
}

void Em5RunAction::SetnbinTh(G4int nbin)
{
  nbinTh = nbin ;
  if(nbinTh>0)
  G4cout << " Nb of bins in Theta plot = " << nbinTh << G4endl ;
}
void Em5RunAction::SetThlow(G4double Tlow)
{
  Thlow = Tlow ;
  if(nbinTh>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlow << G4endl ;
}
void Em5RunAction::SetThhigh(G4double Thigh)
{
  Thhigh = Thigh ;
  if(nbinTh>0)
  G4cout << " Thigh in the Theta plot = " << Thhigh << G4endl ;
}

void Em5RunAction::SetnbinThBack(G4int nbin)
{
  nbinThback = nbin ;
  if(nbinThback>0)
  G4cout << " Nb of bins in Theta plot = " << nbinThback << G4endl ;
}
void Em5RunAction::SetThlowBack(G4double Tlow)
{
  Thlowback = Tlow ;
  if(nbinThback>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlowback << G4endl ;
}
void Em5RunAction::SetThhighBack(G4double Thigh)
{
  Thhighback = Thigh ;
  if(nbinThback>0)
  G4cout << " Thigh in the Theta plot = " << Thhighback << G4endl ;
}
void Em5RunAction::CountParticles(G4double nch,G4double nne)
{
  SumCharged += nch ;
  SumNeutral += nne ;
  Sum2Charged += nch*nch ;
  Sum2Neutral += nne*nne ;
}
void Em5RunAction::AddEP(G4double nele,G4double npos) 
{
  Selectron += nele;
  Spositron += npos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
