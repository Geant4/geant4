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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"
#include "G4PionMinus.hh"
#include "G4EnergyLossTables.hh"
#include "Test17PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17RunAction::Test17RunAction():
   LowestEnergy(0.01*keV),
   HighestEnergy(100.*TeV),
   TotBin(200),
   theProton (G4Proton::Proton()),
   theElectron (G4Electron::Electron()),
   nbinStep(0),nbinEn(0),nbinTt(0),nbinTb(0),
   nbinTsec(0),nbinTh(0),nbinThback(0),nbinR(0),nbinGamma(0),
   nbinvertexz(0),part0(0)
{
  //runMessenger = new Test17RunMessenger(this);

  nbinStep = 100;
  Steplow  = -0.5;
  Stephigh = 95.5;

  nbinEn = 100;
  Enlow  = 0.0;
  Enhigh = 5.0;

  nbinTh = 90;
  Thlow  = 0.0;
  Thhigh = 90.0;

  nbinR  = 100;
  Rlow   = 0.0;
  Rhigh  = 1.0;

  nbinTt = 100;
  Ttlow  = 0.0;
  Tthigh = 2.0;

  nbinThback=90;
  Thlowback =0.0;
  Thhighback=90.0;

  nbinTb=100;
  Tblow =0.0;
  Tbhigh=2.0;

  nbinTsec=100;
  Tseclow =0.0;
  Tsechigh=2.0;

  nbinvertexz=100;
  zlow=0.0;
  zhigh=10.;

  nbinGamma=100;
  ElowGamma=0.0;
  EhighGamma=2.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17RunAction::~Test17RunAction()
{
  //delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::SaveEvent()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::SaveToTuple(G4String parname,G4double val)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::SaveToTuple(G4String parname,G4double val,G4double defval)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //  G4UImanager* UI = G4UImanager::GetUIpointer();
         
  EnergySumAbs = 0. ;
  EnergySquareSumAbs = 0.;
  tlSumAbs = 0. ;
  tlsquareSumAbs = 0. ;
  nStepSumCharged = 0. ;
  nStepSum2Charged= 0. ;
  nStepSumNeutral = 0. ;
  nStepSum2Neutral= 0. ;
  TotNbofEvents = 0. ;
  TotNbofEvents0 = 0. ;
  SumCharged=0.;
  SumNeutral=0.;
  Sum2Charged=0.;
  Sum2Neutral=0.;
  Selectron=0.;
  Spositron=0.;

  Transmitted=0.;
  Reflected  =0.;
//  plot definitions     --------------------------------
  if(nbinStep>0)
  {
    dStep=(Stephigh-Steplow)/nbinStep;
    entryStep=0.;
    underStep=0.;
    overStep=0.;
    for(G4int ist=0; ist<nbinStep; ist++)
    {
      distStep[ist]=0.;
    }
  }      
  if(nbinEn>0)
  {
    dEn = (Enhigh-Enlow)/nbinEn ;
    entryEn=0.;
    underEn=0.;
    overEn=0.;
    for (G4int ien=0; ien<nbinEn; ien++)
    {
      distEn[ien]=0.;
    }
  }

  if(nbinTt>0)
  {
    dTt = (Tthigh-Ttlow)/nbinTt ;
    entryTt=0.;
    underTt=0.;
    overTt=0.;
    for (G4int itt=0; itt<nbinTt; itt++)
    {
      distTt[itt]=0.;
    }
    Ttmean=0.;
    Tt2mean=0.;
  }
  if(nbinTb>0)
  {
    dTb = (Tbhigh-Tblow)/nbinTb ;
    entryTb=0.;
    underTb=0.;
    overTb=0.;
    for (G4int itt=0; itt<nbinTb; itt++)
    {
      distTb[itt]=0.;
    }
    Tbmean=0.;
    Tb2mean=0.;
  }
  if(nbinTsec>0)
  {
    dTsec = (Tsechigh-Tseclow)/nbinTsec ;
    entryTsec=0.;
    underTsec=0.;
    overTsec=0.;
    for (G4int its=0; its<nbinTsec; its++)
    {
      distTsec[its]=0.;
    }
  }
  if(nbinTh>0)
  {
    dTh = (Thhigh-Thlow)/nbinTh ;
    entryTh=0.;
    underTh=0.;
    overTh=0.;
    for (G4int ith=0; ith<nbinTh; ith++)
    {
      distTh[ith]=0.;
    }
  }

  if(nbinThback>0)
  {
    dThback = (Thhighback-Thlowback)/nbinThback ;
    entryThback=0.;
    underThback=0.;
    overThback=0.;
    for (G4int ithback=0; ithback<nbinThback; ithback++)
    {
      distThback[ithback]=0.;
    }
  }


  if(nbinR >0)
  {
    dR  = (Rhigh-Rlow)/nbinR  ;
    entryR =0.;
    underR =0.;
    overR =0.;
    for (G4int ir =0; ir <nbinR ; ir++)
    {
      distR[ir]=0.;
    }
    Rmean=0.;
    R2mean=0.;
  }

  if(nbinGamma>0)
  {
//    dEGamma = log(EhighGamma/ElowGamma)/nbinGamma ;
    dEGamma = (EhighGamma-ElowGamma)/nbinGamma ;
    entryGamma = 0.;
    underGamma=0.;
    overGamma=0.;
    for (G4int ig=0; ig<nbinGamma; ig++)
    {
      distGamma[ig]=0.;
    }
  } 
  if(nbinvertexz>0)
  {
    dz=(zhigh-zlow)/nbinvertexz;
    entryvertexz=0.;
    undervertexz=0.;
    oververtexz=0.;
    for(G4int iz=0; iz<nbinvertexz; iz++)
    {
      distvertexz[iz]=0.;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double sAbs,sigAbs,sigstep,sigcharged,signeutral;

  tlSumAbs /= TotNbofEvents ;
  sAbs = tlsquareSumAbs/TotNbofEvents-tlSumAbs*tlSumAbs ;
  if(sAbs>0.)
    sAbs = sqrt(sAbs/TotNbofEvents) ;
  else
    sAbs = 0. ;
  
  EnergySumAbs /= TotNbofEvents0 ;
  sigAbs = EnergySquareSumAbs/TotNbofEvents0-EnergySumAbs*EnergySumAbs;
  if(sigAbs>0.) sigAbs = sqrt(sigAbs/TotNbofEvents0);
  else          sigAbs = 0.;

  nStepSumCharged /= TotNbofEvents ;
  sigstep = nStepSum2Charged/TotNbofEvents-nStepSumCharged*nStepSumCharged;
  if(sigstep>0.) sigstep = sqrt(sigstep/TotNbofEvents);
  else           sigstep = 0.;
  

  Transmitted /=TotNbofEvents0 ;
  Reflected   /=TotNbofEvents0 ;

  G4bool   icru  = false;
  G4double protL, protR;

  if(abs(kinEnergy0 - 500.0*keV)<0.1*keV && part0 == theProton) { 
    icru  = true;
    protL = 0.009059*mm;
    protR = 0.008869*mm;
  }

  if(abs(kinEnergy0 - MeV)<0.1*keV && part0->GetParticleName() == "pi-") { 
    icru  = true;
    protL = 0.009158*cm*0.9059/0.8869;
    protR = 0.009158*cm;
  }


  G4cout << " ================== run summary =====================" << G4endl;
  G4cout << G4endl;
  G4int prec = G4cout.precision(6);
  G4cout << " end of Run TotNbofEvents = " <<  TotNbofEvents 
         << " for " << part0->GetParticleName() 
         << " with Ekin = " << kinEnergy0/MeV << " MeV" << G4endl ;
  G4cout << "    Track Length in absorber = " <<
          tlSumAbs/micrometer     << " +- " << sAbs/micrometer   <<
          "  microns " << G4endl; 
  G4cout << "    CSDA  Range  in absorber = " <<
          nStepSumCharged/micrometer     << " +- " << sigstep/micrometer   <<
          "  microns " << G4endl; 
  G4cout << G4endl;
  G4cout << "    Energy deposit in absorber = " <<
           EnergySumAbs/MeV << " +- " << sigAbs/MeV <<
          "  MeV " << G4endl ;
  G4cout << G4endl;

  if(icru) {
    G4cout << "### Comparison with the ICRU49 data: " << G4endl;
    G4cout << "    Track Length (G4 - ICRU49) = " 
           << (tlSumAbs - protL)/micrometer
           << " +- " << sAbs/micrometer
           << " microns " << G4endl; 
    G4cout << "    CSDA  Range  (G4 - ICRU49) = " 
           << (nStepSumCharged - protR)/micrometer 
           << " +- " << sigstep/micrometer
           << " microns " << G4endl;
    G4cout << G4endl ;
  }

  G4cout << "(number) transmission coeff=" << Transmitted <<
            "  reflection coeff=" << Reflected << G4endl;
  G4cout << G4endl; 

 G4cout.precision(prec);
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::CountEvent()
{
  TotNbofEvents += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::AddnStepsCharged(G4double ns)
{
  nStepSumCharged += ns;
  nStepSum2Charged += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::AddnStepsNeutral(G4double ns)
{
  nStepSumNeutral += ns;
  nStepSum2Neutral += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::AddEdeps(G4double Eabs)
{
  EnergySumAbs += Eabs;
  EnergySquareSumAbs += Eabs*Eabs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::AddTrackLength(G4double tlabs)
{
  tlSumAbs += tlabs;
  tlsquareSumAbs += tlabs*tlabs ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::AddTrRef(G4double tr,G4double ref)
{
  Transmitted += tr ;
  Reflected   += ref;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillNbOfSteps(G4double)
{
  TotNbofEvents0 += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillTt(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillTb(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillTsec(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Test17RunAction::FillGammaSpectrum(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillThBack(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::FillR(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::Fillvertexz(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Test17RunAction::SetnbinStep(G4int nbin)
{
  nbinStep = nbin ;
  if(nbinStep>0)
  G4cout << " Nb of bins in #step plot = " << nbinStep << G4endl ;
}
void Test17RunAction::SetSteplow(G4double low)
{
  Steplow = low ;
  if(nbinStep>0)
  G4cout << " low  in the #step plot = " << Steplow << G4endl ;
}
void Test17RunAction::SetStephigh(G4double high)
{
  Stephigh = high ;
  if(nbinStep>0)
  G4cout << " high in the #step plot = " << Stephigh << G4endl ;
}
void Test17RunAction::SetnbinEn(G4int nbin)
{
  nbinEn = nbin ;
  if(nbinEn>0)
  G4cout << " Nb of bins in Edep plot = " << nbinEn << G4endl ;
}
void Test17RunAction::SetEnlow(G4double Elow)
{
  Enlow = Elow ;
  if(nbinEn>0)
  G4cout << " Elow  in the  Edep plot = " << Enlow << G4endl ;
}
void Test17RunAction::SetEnhigh(G4double Ehigh)
{
  Enhigh = Ehigh ;
  if(nbinEn>0)
  G4cout << " Ehigh in the  Edep plot = " << Enhigh << G4endl ;
}

void Test17RunAction::SetnbinGamma(G4int nbin)
{
  nbinGamma = nbin ;
  if(nbinGamma>0)
  G4cout << " Nb of bins in gamma spectrum plot = " << nbinGamma << G4endl ;
}
void Test17RunAction::SetElowGamma(G4double Elow)
{
  ElowGamma = Elow ;
  if(nbinGamma>0)
  G4cout << " Elow  in the gamma spectrum plot = " << ElowGamma << G4endl ;
}
void Test17RunAction::SetEhighGamma(G4double Ehigh)
{
  EhighGamma = Ehigh ;
  if(nbinGamma>0)
  G4cout << " Ehigh in the gamma spectrum plot = " << EhighGamma << G4endl ;
}

void Test17RunAction::SetnbinTt(G4int nbin)
{
  nbinTt = nbin ;
  if(nbinTt>0)
  G4cout << " Nb of bins in Etransmisssion plot = " << nbinTt << G4endl ;
}
void Test17RunAction::SetTtlow(G4double Elow)
{
  Ttlow = Elow ;
  if(nbinTt>0)
  G4cout << " Elow  in the  Etransmission plot = " << Ttlow << G4endl ;
}
void Test17RunAction::SetTthigh(G4double Ehigh)
{
  Tthigh = Ehigh ;
  if(nbinTt>0)
  G4cout << " Ehigh in the  Etransmission plot = " << Tthigh << G4endl ;
}
void Test17RunAction::SetnbinTb(G4int nbin)
{
  nbinTb = nbin ;
  if(nbinTb>0)
  G4cout << " Nb of bins in Ebackscattered plot = " << nbinTb << G4endl ;
}
void Test17RunAction::SetTblow(G4double Elow)
{
  Tblow = Elow ;
  if(nbinTb>0)
  G4cout << " Elow  in the  Ebackscattered plot = " << Tblow << G4endl ;
}
void Test17RunAction::SetTbhigh(G4double Ehigh)
{
  Tbhigh = Ehigh ;
  if(nbinTb>0)
  G4cout << " Ehigh in the  Ebackscattered plot = " << Tbhigh << G4endl ;
}

void Test17RunAction::SetnbinTsec(G4int nbin)
{
  nbinTsec = nbin ;
  if(nbinTsec>0)
  G4cout << " Nb of bins in Tsecondary  plot = " << nbinTsec << G4endl ;
}
void Test17RunAction::SetTseclow(G4double Elow)
{
  Tseclow = Elow ;
  if(nbinTsec>0)
  G4cout << " Elow  in the  Tsecondary plot = " << Tseclow << G4endl ;
}
void Test17RunAction::SetTsechigh(G4double Ehigh)
{
  Tsechigh = Ehigh ;
  if(nbinTsec>0)
  G4cout << " Ehigh in the  Tsecondary plot = " << Tsechigh << G4endl ;
}
 
void Test17RunAction::SetnbinR(G4int nbin)
{
  nbinR  = nbin ;
  if(nbinR>0)
  G4cout << " Nb of bins in R plot = " << nbinR << G4endl ;
}
void Test17RunAction::SetRlow(G4double rlow)
{
  Rlow = rlow ;
  if(nbinR>0)
  G4cout << " Rlow  in the  R plot = " << Rlow << G4endl ;
}
void Test17RunAction::SetRhigh(G4double rhigh)
{
  Rhigh = rhigh ;
  if(nbinR>0)
  G4cout << " Rhigh in the R plot = " << Rhigh << G4endl ;
}

void Test17RunAction::Setnbinzvertex(G4int nbin)
{
  nbinvertexz  = nbin ;
  if(nbinvertexz>0)
  G4cout << " Nb of bins in Z plot = " << nbinvertexz << G4endl ;
}
void Test17RunAction::Setzlow(G4double z)
{
  zlow = z ;
  if(nbinvertexz>0)
  G4cout << " zlow  in the  Z plot = " << zlow << G4endl ;
}
void Test17RunAction::Setzhigh(G4double z)
{
  zhigh = z ;
  if(nbinvertexz>0)
  G4cout << " zhigh in the Z plot = " << zhigh << G4endl ;
}

void Test17RunAction::SetnbinTh(G4int nbin)
{
  nbinTh = nbin ;
  if(nbinTh>0)
  G4cout << " Nb of bins in Theta plot = " << nbinTh << G4endl ;
}
void Test17RunAction::SetThlow(G4double Tlow)
{
  Thlow = Tlow ;
  if(nbinTh>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlow << G4endl ;
}
void Test17RunAction::SetThhigh(G4double Thigh)
{
  Thhigh = Thigh ;
  if(nbinTh>0)
  G4cout << " Thigh in the Theta plot = " << Thhigh << G4endl ;
}

void Test17RunAction::SetnbinThBack(G4int nbin)
{
  nbinThback = nbin ;
  if(nbinThback>0)
  G4cout << " Nb of bins in Theta plot = " << nbinThback << G4endl ;
}
void Test17RunAction::SetThlowBack(G4double Tlow)
{
  Thlowback = Tlow ;
  if(nbinThback>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlowback << G4endl ;
}
void Test17RunAction::SetThhighBack(G4double Thigh)
{
  Thhighback = Thigh ;
  if(nbinThback>0)
  G4cout << " Thigh in the Theta plot = " << Thhighback << G4endl ;
}
void Test17RunAction::CountParticles(G4double nch,G4double nne)
{
  SumCharged += nch ;
  SumNeutral += nne ;
  Sum2Charged += nch*nch ;
  Sum2Neutral += nne*nne ;
}
void Test17RunAction::AddEP(G4double nele,G4double npos) 
{
  Selectron += nele;
  Spositron += npos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













