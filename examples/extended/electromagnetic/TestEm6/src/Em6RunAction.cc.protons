// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em6RunAction.hh"
#include "Em6RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include <iomanip.h>

#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>

#include "G4hEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4ionLowEnergyIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6RunAction::Em6RunAction()
  :histName("histfile"),nbinStep(0.),nbinEn(0.),nbinTt(0.),nbinTb(0.),
   nbinTsec(0.),nbinTh(0.),nbinThback(0.),nbinR(0.),nbinGamma(0.),
   nbinvertexz(0.),histo1(0),histo2(0),histo3(0),histo4(0),histo5(0),
   histo6(0),histo7(0),histo8(0),histo9(0),histo10(0),
   theProton (G4Proton::Proton()),
   theElectron ( G4Electron::Electron() ),
   LowestEnergy(0.01*keV),
   HighestEnergy(100.*TeV),
   TotBin(200)
{
  runMessenger = new Em6RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6RunAction::~Em6RunAction()
{
  delete runMessenger;
  if(histo1) delete histo1 ;
  if(histo2) delete histo2 ;
  if(histo3) delete histo3 ;
  if(histo4) delete histo4 ;
  if(histo5) delete histo5 ;
  if(histo6) delete histo6 ;
  if(histo7) delete histo7 ;
  if(histo8) delete histo8 ;
  if(histo9) delete histo9 ;
  if(histo10) delete histo10 ;
  delete histo11 ;
  delete histo12 ;
  delete histo13 ;
  delete histo14 ;
  delete histo15 ;
  delete histo21 ;
  delete histo22 ;
  delete histo23 ;
  delete histo24 ;
  delete histo25 ;
  delete histo31 ;
  delete histo32 ;
  delete histo33 ;
  delete histo34 ;
  delete histo35 ;
  delete histo41 ;
  delete histo42 ;
  delete histo43 ;
  delete histo44 ;
  delete histo45 ;
  delete histo51 ;
  delete histo52 ;
  delete histo53 ;
  delete histo54 ;
  delete histo55 ;
  delete histo61 ;
  delete histo62 ;
  delete histo63 ;
  delete histo64 ;
  delete histo65 ;
  delete histo71 ;
  delete histo72 ;
  delete histo73 ;
  delete histo74 ;
  delete histo75 ;
  delete histo76 ;
  delete histo81 ;
  delete histo82 ;
  delete histo83 ;
  delete histo84 ;
  delete histo85 ;
  delete histo86 ;
  delete hbookManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::bookHisto()
{
  // init hbook
  hbookManager = new HBookFile(histName, 68);
  assert (hbookManager != 0);

  // book histograms
  if(nbinStep>0)
  {
    histo1 = hbookManager->histogram("number of steps/event"
                                   ,nbinStep,Steplow,Stephigh) ;
    assert (histo1 != 0);
  }
  if(nbinEn>0)
  {
    histo2 = hbookManager->histogram("energy deposit in absorber(in MeV)"
                                     ,nbinEn,Enlow,Enhigh) ;
    assert (histo2 != 0);
  }
  if(nbinTh>0)
  {
    histo3 = hbookManager->histogram("angle distribution at exit(deg)"
                                     ,nbinTh,Thlow/deg,Thhigh/deg) ;
    assert (histo3 != 0);
  }
  if(nbinR>0)
  {
    histo4 = hbookManager->histogram("lateral distribution at exit(mm)"
                                     ,nbinR ,Rlow,Rhigh)  ;
    assert (histo4 != 0);
  }
  if(nbinTt>0)
  {
    histo5 = hbookManager->histogram("kinetic energy of the primary at exit(MeV)"
                                     ,nbinTt,Ttlow,Tthigh)  ;
    assert (histo5 != 0);
  }
  if(nbinThback>0)
  {
    histo6 = hbookManager->histogram("angle distribution of backscattered primaries(deg)"
                                     ,nbinThback,Thlowback/deg,Thhighback/deg) ;
    assert (histo6 != 0);
  }
  if(nbinTb>0)
  {
    histo7 = hbookManager->histogram("kinetic energy of the backscattered primaries (MeV)"
                                     ,nbinTb,Tblow,Tbhigh)  ;
    assert (histo7 != 0);
  }
  if(nbinTsec>0)
  {
    histo8 = hbookManager->histogram("kinetic energy of the charged secondaries (MeV)"
                                     ,nbinTsec,Tseclow,Tsechigh)  ;
    assert (histo8 != 0);
  }
  if(nbinvertexz>0)
  {
    histo9 = hbookManager->histogram("z of secondary charged vertices(mm)"
                                     ,nbinvertexz ,zlow,zhigh)  ;
    assert (histo9 != 0);
  }
  if(nbinGamma>0)
  {
    histo10= hbookManager->histogram("kinetic energy of gammas escaping the absorber (MeV)"
                                //     ,nbinGamma,ElowGamma,EhighGamma)  ;
                                ,nbinGamma,log10(ElowGamma),log10(EhighGamma))  ;
    assert (histo10 != 0);
  }

  // Test on G4hLowEnergyIonisation
  histo11 = hbookManager->histogram("proton 40 keV ionisation (keV*cm2/10^15!atoms) Z77p"
                                                  ,92,0.5,92.5) ;
  histo12 = hbookManager->histogram("He 10 keV ionisation (keV*cm2/10^15!atoms) Z77 p"
                                                  ,92,0.5,92.5) ;
  histo13 = hbookManager->histogram("He4 effective charge in Carbon"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo14 = hbookManager->histogram("C12 ionisation in Al (MeV/(mg/cm^2)) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;

  histo15 = hbookManager->histogram("Ar40 ionisation in Al (MeV/(mg/cm^2)) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;

  histo21 = hbookManager->histogram("proton 40 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977p"
                                                  ,92,0.5,92.5) ;
  histo22 = hbookManager->histogram("proton 100 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977p"
                                                  ,92,0.5,92.5) ;
  histo23 = hbookManager->histogram("proton 400 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977p"
                                                  ,92,0.5,92.5) ;
  histo24 = hbookManager->histogram("proton 1 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977p"
                                                  ,92,0.5,92.5) ;
  histo25 = hbookManager->histogram("proton 4 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977p"
                                                  ,92,0.5,92.5) ;
  histo31 = hbookManager->histogram("proton ionisation in Carbon (MeV/mm) Geant4 Ziegler"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo32 = hbookManager->histogram("proton ionisation in Si (MeV/mm) Geant4 Ziegler"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo33 = hbookManager->histogram("proton ionisation in Cu (MeV/mm) Geant4 Ziegler"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo34 = hbookManager->histogram("proton ionisation in CH4 (MeV/mm) Geant4 Ziegler"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo35 = hbookManager->histogram("He effective charge for Au(79)"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;

  histo41 = hbookManager->histogram("proton 40 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU_49p"
                                                  ,92,0.5,92.5) ;
  histo42 = hbookManager->histogram("proton 100 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU_49p"
                                                  ,92,0.5,92.5) ;
  histo43 = hbookManager->histogram("proton 400 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU_49p"
                                                  ,92,0.5,92.5) ;
  histo44 = hbookManager->histogram("proton 1 MeV ionisation (keV*cm2/10^15!atoms) G4 ICRU_49p"
                                                  ,92,0.5,92.5) ;
  histo45 = hbookManager->histogram("proton 4 MeV ionisation (keV*cm2/10^15!atoms) G4 ICRU_49p"
                                                  ,92,0.5,92.5) ;
  histo51 = hbookManager->histogram("proton ionisation in Al (MeV/mm) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo52 = hbookManager->histogram("proton ionisation in Fe (MeV/mm) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo53 = hbookManager->histogram("proton ionisation in Pb (MeV/mm) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo54 = hbookManager->histogram("proton ionisation in H2O (MeV/mm) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;
  histo55 = hbookManager->histogram("proton ionisation in Graphite (MeV/mm) Geant4"
                                   ,TotBin,log10(LowestEnergy),log10(HighestEnergy)) ;


  histo61 = hbookManager->histogram("proton 40 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo62 = hbookManager->histogram("proton 100 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo63 = hbookManager->histogram("proton 400 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo64 = hbookManager->histogram("proton 1 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo65 = hbookManager->histogram("proton 4 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;

  histo71 = hbookManager->histogram("He 40 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo72 = hbookManager->histogram("He 100 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo73 = hbookManager->histogram("He 400 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo74 = hbookManager->histogram("He 1 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo75 = hbookManager->histogram("He 4 MeV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo76 = hbookManager->histogram("He 10 keV ionisation (keV*cm2/10^15!atoms) G4 Ziegler1977He"
                                                  ,92,0.5,92.5) ;
  histo81 = hbookManager->histogram("He 40 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;
  histo82 = hbookManager->histogram("He 100 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;
  histo83 = hbookManager->histogram("He 400 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;
  histo84 = hbookManager->histogram("He 1 MeV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;
  histo85 = hbookManager->histogram("He 4 MeV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;
  histo86 = hbookManager->histogram("He 10 keV ionisation (keV*cm2/10^15!atoms) G4 ICRU49He"
                                                  ,92,0.5,92.5) ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;

  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/clear/view");
    UI->ApplyCommand("/vis/draw/current");
  }
      
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
    dEGamma = log(EhighGamma/ElowGamma)/nbinGamma ;
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

  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::EndOfRunAction(const G4Run* aRun)
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
 G4cout << " ================== run summary =====================" << endl;
 G4int prec = G4cout.precision(6);
  G4cout << " end of Run TotNbofEvents = " <<  
           TotNbofEvents << endl ;
  G4cout << "    mean charged track length   in absorber=" <<
           tlSumAbs/mm      << " +- " << sAbs/mm    <<
          "  mm  " << endl; 
  G4cout << endl;
  G4cout << "            mean energy deposit in absorber=" <<
           EnergySumAbs/MeV << " +- " << sigAbs/MeV <<
          "  MeV " << endl ;
  G4cout << endl ;
  G4cout << " mean number of steps in absorber (charged) =" <<
           nStepSumCharged         << " +- " << sigch     <<
          "      " << endl ;
  G4cout << " mean number of steps in absorber (neutral) =" <<
           nStepSumNeutral         << " +- " << signe     <<
          "      " << endl ;
  G4cout << endl ;
  G4cout << "   mean number of charged secondaries = " <<
           SumCharged << " +- " << sigcharged << endl;  
  G4cout << endl ;
  G4cout << "   mean number of neutral secondaries = " <<
           SumNeutral << " +- " << signeutral << endl;  
  G4cout << endl ;
  
  G4cout << "   mean number of e-s =" << Selectron << 
            "  and e+s =" << Spositron << endl;
  G4cout << endl; 
  
  G4cout << "(number) transmission coeff=" << Transmitted <<
            "  reflection coeff=" << Reflected << endl;
  G4cout << endl; 

  if(nbinStep>0)
  {G4double E , dnorm, norm ;
   G4cout << "   step number/event distribution " << endl ;
   G4cout << "#entries=" << entryStep << "    #underflows=" << underStep <<
             "    #overflows=" << overStep << endl ;
   if( entryStep>0.)
   {
     E = Steplow - dStep ;
     norm = TotNbofEvents ;
     G4cout << " bin nb   nsteplow     entries     normalized " << endl ;
     for(G4int iss=0; iss<nbinStep; iss++)
     {
      E += dStep ;
      dnorm = distStep[iss]/norm;
      G4cout << setw(5) << iss << setw(10) << E << 
                setw(12) << distStep[iss] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }     
  }
  if(nbinEn>0)
  {G4double E , dnorm, norm,fmax,Emp,width ;
   Emp=-999.999 ;
   G4cout << " energy deposit distribution " << endl ;
   G4cout << "#entries=" << entryEn << "    #underflows=" << underEn <<
             "    #overflows=" << overEn << endl ;
   if( entryEn>0.)
   {
     E = Enlow - dEn ;
     norm = TotNbofEvents*dEn ;
     G4cout << " bin nb      Elow      entries     normalized " << endl ;
     fmax = 0. ;
     for(G4int ien=0; ien<nbinEn; ien++)
     {
      E += dEn ;
      if(distEn[ien]>fmax)
      {
        fmax=distEn[ien];
        Emp = E ;
      }
      dnorm = distEn[ien]/norm;
      G4cout << setw(5) << ien << setw(10) << E << 
                setw(12) << distEn[ien] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
     G4int ii ;
     G4double E1,E2 ;
     E1=-1.e6 ;
     E2=+1.e6 ;
     E = Enlow -dEn ;
     ii = -1;
     for(G4int i1=0; i1<nbinEn; i1++)
     {
      E += dEn ;
     if(ii<0)
     {
      if(distEn[i1] >= 0.5*fmax)
      {
        E1=E ;
        ii=i1 ;
      }
     }
     }
     E = Enlow -dEn ;
     for(G4int i2=0; i2<nbinEn; i2++)
     {
      E += dEn ;
      if(distEn[i2] >= 0.5*fmax)
        E2=E ;
     }
     G4cout << " Emp = " << setw(15) << Emp/MeV << "   width="
            << setw(15) << (E2-E1)/MeV <<   "  MeV " << endl;
     G4cout << endl ;
   }     
  }
  if(nbinTt>0)
  {G4double E , dnorm, norm ,sig;
   G4cout << " transmitted energy distribution " << endl ;
   G4cout << "#entries=" << entryTt << "    #underflows=" << underTt <<
             "    #overflows=" << overTt << endl ;
   if( entryTt>0.)
   {
     Ttmean /= entryTt;
     sig=Tt2mean/entryTt-Ttmean*Ttmean ;
     if(sig<=0.)
       sig=0.;
     else
       sig=sqrt(sig/entryTt) ;
     G4cout << " mean energy of transmitted particles=" << Ttmean/keV << 
               " +- " << sig/keV << "  keV." << endl;
     E = Ttlow - dTt ;
     norm = TotNbofEvents*dTt ;
     G4cout << " bin nb      Elow      entries     normalized " << endl ;
     for(G4int itt=0; itt<nbinTt; itt++)
     {
      E += dTt ;
      dnorm = distTt[itt]/norm;
      G4cout << setw(5) << itt << setw(10) << E << 
                setw(12) << distTt[itt] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }     
  }
  if(nbinTb>0)
  {G4double E , dnorm, norm ,sig;
   G4cout << " backscattered energy distribution " << endl ;
   G4cout << "#entries=" << entryTb << "    #underflows=" << underTb <<
             "    #overflows=" << overTb << endl ;
   if( entryTb>0.)
   {
     Tbmean /= entryTb;
     sig=Tb2mean/entryTb-Tbmean*Tbmean ;
     if(sig<=0.)
       sig=0.;
     else
       sig=sqrt(sig/entryTb) ;
     G4cout << " mean energy of backscattered particles=" << Tbmean/keV << 
               " +- " << sig/keV << "  keV." << endl;
     E = Tblow - dTb ;
     norm = TotNbofEvents*dTb ;
     G4cout << " bin nb      Elow      entries     normalized " << endl ;
     for(G4int itt=0; itt<nbinTb; itt++)
     {
      E += dTb ;
      dnorm = distTb[itt]/norm;
      G4cout << setw(5) << itt << setw(10) << E << 
                setw(12) << distTb[itt] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }     
  }
  if(nbinTsec>0)
  {G4double E , dnorm, norm ;
   G4cout << " energy distribution of charged secondaries " << endl ;
   G4cout << "#entries=" << entryTsec << "    #underflows=" << underTsec <<
             "    #overflows=" << overTsec << endl ;
   if( entryTsec>0.)
   {
     E = Tseclow - dTsec ;
     norm = TotNbofEvents*dTsec ;
     G4cout << " bin nb      Elow      entries     normalized " << endl ;
     for(G4int itt=0; itt<nbinTsec; itt++)
     {
      E += dTsec ;
      dnorm = distTsec[itt]/norm;
      G4cout << setw(5) << itt << setw(10) << E << 
                setw(12) << distTsec[itt] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }     
  }

  if(nbinR >0)
  {G4double R , dnorm, norm,sig  ;
   G4cout << "  R  distribution " << endl ;
   G4cout << "#entries=" << entryR  << "    #underflows=" << underR  <<
             "    #overflows=" << overR  << endl ;
   if( entryR >0.)
   {
     Rmean /= entryR;
     sig = R2mean/entryR - Rmean*Rmean;
     if(sig<=0.) sig=0. ;
     else        sig = sqrt(sig/entryR) ;
     G4cout << " mean lateral displacement at exit=" << Rmean/mm << " +- "
            << sig/mm << "  mm." << endl ; 
     R = Rlow - dR  ;
     norm = TotNbofEvents*dR  ;
     G4cout << " bin nb      Rlow      entries     normalized " << endl ;
     for(G4int ier=0; ier<nbinR ; ier++)
     {
      R+= dR  ;
      dnorm = distR[ier]/norm;
      G4cout << setw(5) << ier << setw(10) << R  <<
                setw(12) << distR[ier] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }
  }

  if(nbinTh>0)
  {G4double Th,Thdeg, dnorm, norm,fac0,fnorm,pere,Thpere,Thmean,sum;
   G4cout << "      angle   distribution " << endl ;
   G4cout << "#entries=" << entryTh << "    #underflows=" << underTh <<
             "    #overflows=" << overTh << endl ;
   if( entryTh>0.)
   {
     Th= Thlow - dTh ;
     norm = TotNbofEvents ;
     if(distTh[0] == 0.)
       fac0 = 1. ;
     else
       fac0 = 1./distTh[0] ;
     pere = 1./exp(1.) ;

     G4cout << " bin nb  Thlowdeg      Thlowrad      " <<
               " entries         normalized " << endl ;
     Thpere = 0. ;
     sum = 0. ;
     Thmean = 0. ;
     for(G4int ien=0; ien<nbinTh; ien++)
     {
      Th+= dTh ;
      Thdeg = Th*180./pi ;
      sum += distTh[ien] ;
      Thmean += distTh[ien]*(Th+0.5*dTh) ;
      dnorm = distTh[ien]/norm;
      fnorm = fac0*distTh[ien] ;
      if( fnorm > pere)
        Thpere = Th ; 
      G4cout << setw(5) << ien << setw(10) << Thdeg << "   " <<
                setw(10) << Th << "  " <<   
                setw(12) << distTh[ien] << "  " <<
                setw(12) << dnorm << "  " << setw(12) << fnorm <<endl ;
     }
     Thmean /= sum ;
     G4cout << endl;
     G4cout << " mean = " << Thmean << "  rad  or " << 180.*Thmean/pi <<
               " deg." << endl;
     G4cout << " theta(1/e)=" << Thpere << " - " << Thpere+dTh << " rad   "
            << " or " << 180.*Thpere/pi << " - " << 180.*(Thpere+dTh)/pi 
            << " deg." << endl;
     G4cout << endl;
   }
  }

  if(nbinThback>0)
  {G4double Thb,Thdegb, dnormb, normb,fac0b,fnormb,pereb,Thpereb,Thmeanb,sumb;
   G4cout << " backscattering angle   distribution " << endl ;
   G4cout << "#entries=" << entryThback << "    #underflows=" << underThback <<
             "    #overflows=" << overThback << endl ;
   if( entryThback>0.)
   {
     Thb= Thlowback - dThback ;
     normb = TotNbofEvents ;
     if(distThback[0] == 0.)
       fac0b = 1. ;
     else
       fac0b = 1./distThback[0] ;
     pereb = 1./exp(1.) ;

     G4cout << " bin nb  Thlowdeg      Thlowrad      " <<
               " entries         normalized " << endl ;
     Thpereb = 0. ;
     sumb = 0. ;
     Thmeanb = 0. ;
     for(G4int ien=0; ien<nbinThback; ien++)
     {
      Thb+= dThback ;
      Thdegb = Thb*180./pi ;
      sumb += distThback[ien] ;
      Thmeanb += distThback[ien]*(Thb+0.5*dThback) ;
      dnormb = distThback[ien]/normb;
      fnormb = fac0b*distThback[ien] ;
      if( fnormb > pereb)
        Thpereb = Thb ;
      G4cout << setw(5) << ien << setw(10) << Thdegb << "   " <<
                setw(10) << Thb << "  " <<
                setw(12) << distThback[ien] << "  " <<
                setw(12) << dnormb << "  " << setw(12) << fnormb <<endl ;
     }
     Thmeanb /= sumb ;
     G4cout << endl;
     G4cout << " mean = " << Thmeanb << "  rad  or " << 180.*Thmeanb/pi <<
               " deg." << endl;
     G4cout << " theta(1/e)=" << Thpereb << " - " << Thpereb+dThback << " rad   "
            << " or " << 180.*Thpereb/pi << " - " << 180.*(Thpereb+dThback)/pi
            << " deg." << endl;
     G4cout << endl;
   }
  }

  if(nbinGamma>0)
  {G4double E , fact,dnorm, norm  ;
   G4cout << " gamma energy distribution " << endl ;
   G4cout << "#entries=" << entryGamma << "    #underflows=" << underGamma <<
             "    #overflows=" << overGamma << endl ;
   if( entryGamma>0.)
   {
     fact=exp(dEGamma) ;
     E = ElowGamma/fact  ;
     norm = TotNbofEvents*dEGamma;
     G4cout << " bin nb         Elow      entries       normalized " << endl ;
     for(G4int itt=0; itt<nbinGamma; itt++)
     {
      E *= fact ;
      dnorm = distGamma[itt]/norm;
      G4cout << setw(5) << itt << setw(13) << E << 
                setw(12) << distGamma[itt] <<
                setw(15) << dnorm << endl ;
     }
     G4cout << endl;
   }     
  }

  if(nbinvertexz >0)
  {G4double z , dnorm, norm  ;
   G4cout << " vertex Z  distribution " << endl ;
   G4cout << "#entries=" << entryvertexz  << "    #underflows=" << undervertexz  <<
             "    #overflows=" << oververtexz  << endl ;
   if( entryvertexz >0.)
   {
     z =zlow - dz  ;
     norm = TotNbofEvents*dz  ;
     G4cout << " bin nb      zlow      entries     normalized " << endl ;
     for(G4int iez=0; iez<nbinvertexz ; iez++)
     {
      z+= dz  ;
      if(abs(z)<1.e-12) z=0.;
      dnorm = distvertexz[iez]/norm;
      G4cout << setw(5) << iez << setw(10) << z  <<
                setw(12) << distvertexz[iez] <<
                setw(12) << dnorm << endl ;
     }
     G4cout << endl;
   }
  }
  
 G4cout.precision(prec);
  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
    
   // Write histogram file
   FillLowEnergyTest() ;
   hbookManager->write();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::CountEvent()
{
  TotNbofEvents += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::AddnStepsCharged(G4double ns)
{
  nStepSumCharged += ns;
  nStepSum2Charged += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::AddnStepsNeutral(G4double ns)
{
  nStepSumNeutral += ns;
  nStepSum2Neutral += ns*ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::AddEdeps(G4double Eabs)
{
  EnergySumAbs += Eabs;
  EnergySquareSumAbs += Eabs*Eabs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::AddTrackLength(G4double tlabs)
{
  tlSumAbs += tlabs;
  tlsquareSumAbs += tlabs*tlabs ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::AddTrRef(G4double tr,G4double ref)
{
  Transmitted += tr ;
  Reflected   += ref;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillNbOfSteps(G4double ns)
{
  const G4double eps = 1.e-10 ;
  G4double n,bin ;
  G4int ibin;
 
  if(histo1)
  {
    entryStep += 1. ;
 
    if(ns<Steplow)
      underStep += 1. ;
    else if(ns>=Stephigh)
      overStep  += 1. ;
    else
    {
      n = ns+eps ;
      bin = (n-Steplow)/dStep ;
      ibin= bin ;
      distStep[ibin] += 1. ;
    }
   histo1->accumulate(ns) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillEn(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo2)
  {
    entryEn += 1. ;
 
    if(En<Enlow)
      underEn += 1. ;
    else if(En>=Enhigh)
      overEn  += 1. ;
    else
    {
      bin = (En-Enlow)/dEn ;
      ibin= bin ;
      distEn[ibin] += 1. ;
    }
  histo2->accumulate(En/MeV) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillTt(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo5)
  {
    entryTt += 1. ;
    Ttmean += En ;
    Tt2mean += En*En ;

    if(En<Ttlow)
      underTt += 1. ;
    else if(En>=Tthigh)
      overTt  += 1. ;
    else
    {
      bin = (En-Ttlow)/dTt ;
      ibin= bin ;
      distTt[ibin] += 1. ;
    }
  histo5->accumulate(En/MeV) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillTb(G4double En)
{
  G4double bin ;
  G4int ibin;
  
  if(histo7)
  {
    entryTb += 1. ;
    Tbmean += En ;
    Tb2mean += En*En ;

    if(En<Tblow)
      underTb += 1. ;
    else if(En>=Tbhigh)
      overTb  += 1. ;
    else
    {
      bin = (En-Tblow)/dTb ;
      ibin= bin ;
      distTb[ibin] += 1. ;
    }
  histo7->accumulate(En/MeV) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillTsec(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo8)
  {
    entryTsec += 1. ;

    if(En<Tseclow)
      underTsec += 1. ;
    else if(En>=Tsechigh)
      overTsec  += 1. ;
    else
    {
      bin = (En-Tseclow)/dTsec ;
      ibin= bin ;
      distTsec[ibin] += 1. ;
    }
  histo8->accumulate(En/MeV) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillGammaSpectrum(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo10)
  {
    entryGamma += 1. ;

    if(En<ElowGamma)
      underGamma += 1. ;
    else if(En>=EhighGamma)
      overGamma  += 1. ;
    else
    {
      bin = log(En/ElowGamma)/dEGamma;
      ibin= bin ;
      distGamma[ibin] += 1. ;
    }
  histo10->accumulate(log10(En/MeV)) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillTh(G4double Th)
{
  static const G4double cn=pi/(64800.*dTh) ;
  static const G4double cs=pi/
        (64800.*(cos(Thlow)-cos(Thlow+dTh)));      
  G4double bin,Thbin ,wg;
  G4int ibin;

  if(histo3)
  {
    entryTh += 1. ;

    wg = 0.;

    if(Th<Thlow)
      underTh += 1. ;
    else if(Th>=Thhigh)
      overTh  += 1. ;
    else
    {
      bin = (Th-Thlow)/dTh ;
      ibin= bin ;
      Thbin = Thlow+ibin*dTh ;
      if(Th > 0.001*dTh)
        wg=cn/sin(Th) ;
      else
      {  
        G4double thdeg=Th*180./pi;
        G4cout << "theta < 0.001*dth (from plot excluded) theta="
               << setw(12) << setprecision(4) << thdeg << endl;
        wg=0. ; 
      }
      distTh[ibin] += wg  ;
    }

  histo3->accumulate(Th/deg, wg) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillThBack(G4double Th)
{
  static const G4double cn=pi/(64800.*dThback) ;
  static const G4double cs=pi/
        (64800.*(cos(Thlowback)-cos(Thlowback+dThback)));      
  G4double bin,Thbin,wg ;
  G4int ibin;

  if(histo6)
  {
    entryThback += 1. ;

    if(Th<Thlowback)
      underThback += 1. ;
    else if(Th>=Thhighback)
      overThback  += 1. ;
    else
    {
      bin = (Th-Thlowback)/dThback ;
      ibin= bin ;
      Thbin = Thlowback+ibin*dThback ;
      if(Th > 0.001*dThback)
        wg=cn/sin(Th) ;
      else
      {  
        G4double thdeg=Th*180./pi;
        G4cout << "theta < 0.001*dth (from plot excluded) theta="
               << setw(12) << setprecision(4) << thdeg << endl;
        wg=0. ; 
      }
      distThback[ibin] += wg  ;
    }
  histo6->accumulate(Th/deg, wg) ;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillR(G4double R )
{
  G4double bin ;
  G4int ibin;

  if(histo4)
  {
    entryR  += 1. ;
    Rmean += R ;
    R2mean += R*R ;

    if(R <Rlow)
      underR  += 1. ;
    else if(R >=Rhigh)
      overR   += 1. ;
    else
    {
      bin = (R -Rlow)/dR  ;
      ibin= bin ;
      distR[ibin] += 1. ;
    }
  histo4->accumulate(R/mm) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::Fillvertexz(G4double z )
{
  G4double bin ;
  G4int ibin;
  
  if(histo9)
  {
    entryvertexz  += 1. ;

    if(z <zlow)
      undervertexz  += 1. ;
    else if(z >=zhigh)
      oververtexz   += 1. ;
    else
    {
      bin = (z -zlow)/dz  ;
      ibin= bin ;
      distvertexz[ibin] += 1. ;
    }
  histo9->accumulate(z/mm) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::SethistName(G4String name)
{
  histName = name ;
  G4cout << " hist file = " << histName << endl;
}

void Em6RunAction::SetnbinStep(G4int nbin)
{
  nbinStep = nbin ;
  if(nbinStep>0)
  G4cout << " Nb of bins in #step plot = " << nbinStep << endl ;
}
void Em6RunAction::SetSteplow(G4double low)
{
  Steplow = low ;
  if(nbinStep>0)
  G4cout << " low  in the #step plot = " << Steplow << endl ;
}
void Em6RunAction::SetStephigh(G4double high)
{
  Stephigh = high ;
  if(nbinStep>0)
  G4cout << " high in the #step plot = " << Stephigh << endl ;
}
void Em6RunAction::SetnbinEn(G4int nbin)
{
  nbinEn = nbin ;
  if(nbinEn>0)
  G4cout << " Nb of bins in Edep plot = " << nbinEn << endl ;
}
void Em6RunAction::SetEnlow(G4double Elow)
{
  Enlow = Elow ;
  if(nbinEn>0)
  G4cout << " Elow  in the  Edep plot = " << Enlow << endl ;
}
void Em6RunAction::SetEnhigh(G4double Ehigh)
{
  Enhigh = Ehigh ;
  if(nbinEn>0)
  G4cout << " Ehigh in the  Edep plot = " << Enhigh << endl ;
}

void Em6RunAction::SetnbinGamma(G4int nbin)
{
  nbinGamma = nbin ;
  if(nbinGamma>0)
  G4cout << " Nb of bins in gamma spectrum plot = " << nbinGamma << endl ;
}
void Em6RunAction::SetElowGamma(G4double Elow)
{
  ElowGamma = Elow ;
  if(nbinGamma>0)
  G4cout << " Elow  in the gamma spectrum plot = " << ElowGamma << endl ;
}
void Em6RunAction::SetEhighGamma(G4double Ehigh)
{
  EhighGamma = Ehigh ;
  if(nbinGamma>0)
  G4cout << " Ehigh in the gamma spectrum plot = " << EhighGamma << endl ;
}

void Em6RunAction::SetnbinTt(G4int nbin)
{
  nbinTt = nbin ;
  if(nbinTt>0)
  G4cout << " Nb of bins in Etransmisssion plot = " << nbinTt << endl ;
}
void Em6RunAction::SetTtlow(G4double Elow)
{
  Ttlow = Elow ;
  if(nbinTt>0)
  G4cout << " Elow  in the  Etransmission plot = " << Ttlow << endl ;
}
void Em6RunAction::SetTthigh(G4double Ehigh)
{
  Tthigh = Ehigh ;
  if(nbinTt>0)
  G4cout << " Ehigh in the  Etransmission plot = " << Tthigh << endl ;
}
void Em6RunAction::SetnbinTb(G4int nbin)
{
  nbinTb = nbin ;
  if(nbinTb>0)
  G4cout << " Nb of bins in Ebackscattered plot = " << nbinTb << endl ;
}
void Em6RunAction::SetTblow(G4double Elow)
{
  Tblow = Elow ;
  if(nbinTb>0)
  G4cout << " Elow  in the  Ebackscattered plot = " << Tblow << endl ;
}
void Em6RunAction::SetTbhigh(G4double Ehigh)
{
  Tbhigh = Ehigh ;
  if(nbinTb>0)
  G4cout << " Ehigh in the  Ebackscattered plot = " << Tbhigh << endl ;
}

void Em6RunAction::SetnbinTsec(G4int nbin)
{
  nbinTsec = nbin ;
  if(nbinTsec>0)
  G4cout << " Nb of bins in Tsecondary  plot = " << nbinTsec << endl ;
}
void Em6RunAction::SetTseclow(G4double Elow)
{
  Tseclow = Elow ;
  if(nbinTsec>0)
  G4cout << " Elow  in the  Tsecondary plot = " << Tseclow << endl ;
}
void Em6RunAction::SetTsechigh(G4double Ehigh)
{
  Tsechigh = Ehigh ;
  if(nbinTsec>0)
  G4cout << " Ehigh in the  Tsecondary plot = " << Tsechigh << endl ;
}
 
void Em6RunAction::SetnbinR(G4int nbin)
{
  nbinR  = nbin ;
  if(nbinR>0)
  G4cout << " Nb of bins in R plot = " << nbinR << endl ;
}
void Em6RunAction::SetRlow(G4double rlow)
{
  Rlow = rlow ;
  if(nbinR>0)
  G4cout << " Rlow  in the  R plot = " << Rlow << endl ;
}
void Em6RunAction::SetRhigh(G4double rhigh)
{
  Rhigh = rhigh ;
  if(nbinR>0)
  G4cout << " Rhigh in the R plot = " << Rhigh << endl ;
}

void Em6RunAction::Setnbinzvertex(G4int nbin)
{
  nbinvertexz  = nbin ;
  if(nbinvertexz>0)
  G4cout << " Nb of bins in Z plot = " << nbinvertexz << endl ;
}
void Em6RunAction::Setzlow(G4double z)
{
  zlow = z ;
  if(nbinvertexz>0)
  G4cout << " zlow  in the  Z plot = " << zlow << endl ;
}
void Em6RunAction::Setzhigh(G4double z)
{
  zhigh = z ;
  if(nbinvertexz>0)
  G4cout << " zhigh in the Z plot = " << zhigh << endl ;
}

void Em6RunAction::SetnbinTh(G4int nbin)
{
  nbinTh = nbin ;
  if(nbinTh>0)
  G4cout << " Nb of bins in Theta plot = " << nbinTh << endl ;
}
void Em6RunAction::SetThlow(G4double Tlow)
{
  Thlow = Tlow ;
  if(nbinTh>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlow << endl ;
}
void Em6RunAction::SetThhigh(G4double Thigh)
{
  Thhigh = Thigh ;
  if(nbinTh>0)
  G4cout << " Thigh in the Theta plot = " << Thhigh << endl ;
}

void Em6RunAction::SetnbinThBack(G4int nbin)
{
  nbinThback = nbin ;
  if(nbinThback>0)
  G4cout << " Nb of bins in Theta plot = " << nbinThback << endl ;
}
void Em6RunAction::SetThlowBack(G4double Tlow)
{
  Thlowback = Tlow ;
  if(nbinThback>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlowback << endl ;
}
void Em6RunAction::SetThhighBack(G4double Thigh)
{
  Thhighback = Thigh ;
  if(nbinThback>0)
  G4cout << " Thigh in the Theta plot = " << Thhighback << endl ;
}
void Em6RunAction::CountParticles(G4double nch,G4double nne)
{
  SumCharged += nch ;
  SumNeutral += nne ;
  Sum2Charged += nch*nch ;
  Sum2Neutral += nne*nne ;
}
void Em6RunAction::AddEP(G4double nele,G4double npos) 
{
  Selectron += nele;
  Spositron += npos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6RunAction::FillLowEnergyTest( )
{
  G4int J ;
  G4double* DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
  G4hLowEnergyIonisation* hLEIon = new G4hLowEnergyIonisation() ;
  G4ParticleDefinition* theH = G4Proton::Proton() ;

  //G4double fac = 4.026/1.0073 ;
  //G4double cac = 1.0 ;
  //  theProton = G4Proton::Proton() ;
    G4double fac = 1.0 ;
    G4double cac = 2.0 ;
    theProton = G4Alpha::Alpha() ;

    G4ionLowEnergyIonisation* ionLEIon = new G4ionLowEnergyIonisation() ;
    ionLEIon->SetIonDefinition(theProton) ;

    //  G4ionLowEnergyIonisation* ionLEIonC = new G4ionLowEnergyIonisation() ;
    //  G4ParticleDefinition* theC12 = G4IonC12::IonC12() ;
    //  ionLEIonC->SetIonDefinition(theC12) ;
  G4double chc = 6.0 ;

  //  G4ionLowEnergyIonisation* ionLEIonAr = new G4ionLowEnergyIonisation() ;
  //  G4ParticleDefinition* theAr40 = G4IonAr40::IonAr40() ;
  //  ionLEIonAr->SetIonDefinition(theAr40) ;
  G4double cha = 18.0 ;


  //  G4ionLowEnergyIonisation* ionLEIonFe = new G4ionLowEnergyIonisation() ;
  //  G4ParticleDefinition* theFe56 = G4IonFe56::IonFe56() ;
  //  ionLEIonFe->SetIonDefinition(theFe56) ;
  //  G4double chf = 26.0 ;

  //  G4double mrC12 = (theC12->GetPDGMass()) / (theH->GetPDGMass()) ; 
  //  G4double mrAr40 = (theAr40->GetPDGMass()) / (theH->GetPDGMass()) ; 

  G4double de ;

  // Ziegler1977p
  G4double tau = 40*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPower1977H(J, tau) ;
    histo21->accumulate(double(J),de) ;
  }
  tau = 100*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPower1977H(J, tau) ;
    histo22->accumulate(double(J),de) ;
  }
  tau = 400*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPower1977H(J, tau) ;
    histo23->accumulate(double(J),de) ;
  }
  tau = 1000*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPower1977H(J, tau) ;
    histo24->accumulate(double(J),de) ;
  }
  tau = 4000*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPower1977H(J, tau) ;
    histo25->accumulate(double(J),de) ;
  }

  // Ziegler1977He
  tau = 40*keV*4.026/1.0073 ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) /
                  (hLEIon->GetHeEffChargeSquare(J, tau)) ;
    histo61->accumulate(double(J),de) ;
  }
  tau = 100*keV*4.026/1.0073 ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) /
                  (hLEIon->GetHeEffChargeSquare(J, tau)) ;
    histo62->accumulate(double(J),de) ;
  }
  tau = 400*keV*4.026/1.0073 ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) /
                  (hLEIon->GetHeEffChargeSquare(J, tau)) ;
    histo63->accumulate(double(J),de) ;
  }
  tau = 1000*keV*4.026/1.0073 ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) /
                  (hLEIon->GetHeEffChargeSquare(J, tau)) ;
    histo64->accumulate(double(J),de) ;
  }
  tau = 4000*keV*4.026/1.0073 ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) /
                  (hLEIon->GetHeEffChargeSquare(J, tau)) ;
    histo65->accumulate(double(J),de) ;
  }

  // ICRU_49p
  tau = 40*keV ;
 
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPowerICRU_R49p(J, tau, "Ele") ;
    histo41->accumulate(double(J),de) ;
  }
  tau = 100*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPowerICRU_R49p(J, tau, "Ele") ;
    histo42->accumulate(double(J),de) ;
  }
  tau = 400*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPowerICRU_R49p(J, tau, "Ele") ;
    histo43->accumulate(double(J),de) ;
  }
  tau = 1000*keV ;
  for ( J=1; J<93; J++)
  { 
    de = hLEIon->GetStoppingPowerICRU_R49p(J, tau, "Ele") ;
    histo44->accumulate(double(J),de) ;
  }
  tau = 4000*keV ;
   for ( J=1; J<93; J++)
  { 
   de = hLEIon->GetStoppingPowerICRU_R49p(J, tau, "Ele") ;
   histo45->accumulate(double(J),de) ;
  }
 
  for ( J=1; J<93; J++)
  { 
    static G4double S40[92] = {
    6.22, 6.55, 7.61, 10.4, 12.9, 13.9, 15.9, 14.6, 11.7, 11.1,
    14.2, 20.6, 20.7,  22.3, 18.1, 19.3, 27.5, 30.5, 27.9, 29.8,
    28.3, 26.7, 25.1,  22.5, 19.8, 20.1, 18.0, 20.1, 20.4, 23.4,
    27.5, 29.0, 29.3,  31.0, 31.1, 34.8, 31.5, 35.0, 35.5, 37.3,
    38.3, 35.9, 38.0,  34.4, 33.4, 29.8, 30.9, 32.7, 34.9, 36.0,
    40.9, 39.1, 43.1,  45.8, 40.8, 44.1, 44.9, 42.0, 41.0, 40.0,
    39.0, 38.0, 37.1,  38.2, 35.3, 31.5, 29.9, 29.1, 28.3, 27.5, 
    28.1, 28.9, 27.3,  26.3, 29.9, 29.1, 28.4, 25.8, 28.0, 24.9,
    27.3, 30.7, 34.3,  35.4, 35.6, 35.5, 39.8, 42.9, 43.7, 44.1,
    42.4, 41.8} ;
    histo11->accumulate(double(J),S40[J-1]) ;
  }

  //============================ He ions ===========
  // Ziegler1977He
  tau = 40*keV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo71->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo81->accumulate(double(J),de) ;
  }
  tau = 100*keV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo72->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo82->accumulate(double(J),de) ;
  }
  tau = 400*keV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo73->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo83->accumulate(double(J),de) ;
  }
  tau = 1000*keV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo74->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo84->accumulate(double(J),de) ;
  }
  tau = 4*MeV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo75->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo85->accumulate(double(J),de) ;
  }
  tau = 10*keV ;
  for ( J=1; J<93; J++)
  { 
    de = (hLEIon->GetStoppingPower1977He(J, tau)) ;
    histo76->accumulate(double(J),de) ;
    de = hLEIon->GetStoppingPowerICRU_R49He(J, tau ) ;
    histo86->accumulate(double(J),de) ;
  }
 
  for ( J=1; J<93; J++)
  { 
    static G4double Q40[92] = {
    11.8, 16.7, 21.9,  24.1, 34.7, 36.1, 45.5, 46.1, 45.8, 44.9,
    46.9, 51.7, 54.9,  63.2, 63.8, 67.4, 79.8, 80.8, 78.0, 83.2,
    85.4, 84.2, 90.6,  85.7, 84.1, 86.4, 82.0, 77.7, 74.4, 75.1,
    75.6, 80.7, 84.9,  87.0, 88.6, 103.0, 103.0, 113.0, 116.0, 123.0,
    120.0, 115.0, 118.0, 115.0, 113.0, 112.0, 111.0, 110.0, 116.0, 117.0,
    119.0, 123.0, 122.0, 139.0, 138.0, 143.0, 152.0, 141.0, 139.0, 137.0,
    135.0, 134.0, 130.0, 137.0, 129.0, 128.0, 124.0, 125.0, 120.0, 118.0, 
    119.0, 120.0, 121.0, 122.0, 125.0, 127.0, 128.0, 120.0, 125.0, 128.0,
    136.0, 138.0, 144.0, 148.0, 151.0, 162.0, 166.0, 177.0, 179.0, 183.0, 
    177.0, 176.0 } ;   
    histo12->accumulate(double(J),Q40[J-1]) ;
  }
  //============================ He ions ===========


  //  find out materials table
  static const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();

  G4int numOfMaterials = theMaterialTable->length();
  G4double dedx ;
    // create physics vector and fill it

  G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowestEnergy, HighestEnergy, TotBin);


  //  loop for materials - the parametrisation is defined in the G4PhysicsList

  for ( J=0; J<numOfMaterials; J++)
  {
    G4Material* material = (*theMaterialTable)[J];
    if ("Carbon" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5*(aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	//	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
	G4double cf = ionLEIon->GetHeEffChargeSquare(6, tau) / (cac*cac);
        histo31->accumulate(log10(tau),dedx*cf) ;
        histo13->accumulate(log10(tau),cf*cac*cac) ;
      }
    } else if ("Silicon" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo32->accumulate(log10(tau),dedx*cf) ;
      }
    } else if ("Copper" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo33->accumulate(log10(tau),dedx*cf) ;
      }
    } else if ("CH_4" == material->GetChemicalFormula()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
	tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo34->accumulate(log10(tau),dedx*cf) ;
	cf = hLEIon->GetHeEffChargeSquare(79, tau*4.026/1.0073) ;
        histo35->accumulate(log10(tau),cf) ;
      }
    } else if ("Aluminum" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
	dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo51->accumulate(log10(tau),dedx*cf) ;
	/*
	dedx = G4EnergyLossTables::GetPreciseDEDX(theC12,tau*mrC12,material) ;
	cf = ionLEIonC->GetIonEffChargeSquare(material, tau*mrC12, chc) / (chc*chc);
        histo14->accumulate(log10(tau), dedx*cf / 270.0  ) ;
	dedx = G4EnergyLossTables::GetPreciseDEDX(theAr40,tau*mrAr40,material) ;
	cf = ionLEIonAr->GetIonEffChargeSquare(material, tau*mrAr40, cha) / (cha*cha);
        histo15->accumulate(log10(tau), dedx*cf / 270.0  ) ;
	*/
      }
    } else if ("Iron" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo52->accumulate(log10(tau),dedx*cf) ;
      }
    } else if ("Lead" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	//	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
	G4double cf = hLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac);
        histo53->accumulate(log10(tau),dedx*cf) ;
      }
    } else if ("Water" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
      {
        tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
        dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac) ;
	//        G4cout << "Energy = " << tau << " MeV; dedx = " << dedx << endl;
        histo54->accumulate(log10(tau),dedx*cf) ;
      }
    } else if ("Graphite" == material->GetName()) {
      for (G4int i = 0 ; i < TotBin-1 ; i++)
	{
	  tau = 0.5 * (aVector->GetLowEdgeEnergy(i) + aVector->GetLowEdgeEnergy(i+1)) ;
	  dedx = G4EnergyLossTables::GetPreciseDEDX(theProton,tau/fac,material) ;
	G4double cf = ionLEIon->GetIonEffChargeSquare(material, tau, cac) / (cac*cac) ;
	  histo55->accumulate(log10(tau),dedx*cf) ;
	}
    }
  } 
    delete hLEIon ;
    delete ionLEIon ;
    //    delete ionLEIonC ;
    //    delete ionLEIonAr ;
    //    delete ionLEIonFe ;
}













