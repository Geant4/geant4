// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8RunAction.cc,v 1.3 2000-06-27 10:51:11 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#include "Em8RunAction.hh"
#include "Em8RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include <g4std/iomanip>

#include "Randomize.hh"
#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>

//////////////////////////////////////////////////////////////////////////////

Em8RunAction::Em8RunAction()
  :histName("histfile"),nbinStep(0),nbinEn(0),nbinTt(0),nbinTb(0),
   nbinTsec(0),nbinTh(0),nbinThback(0),nbinR(0),nbinGamma(0),
   nbinvertexz(0),histo1(0),histo2(0),histo3(0),histo4(0),histo5(0),
   histo6(0),histo7(0),histo8(0),histo9(0),histo10(0)
{
  runMessenger = new Em8RunMessenger(this);
  saveRndm = 1;  
}

////////////////////////////////////////////////////////////////////////////

Em8RunAction::~Em8RunAction()
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
  delete hbookManager;

}

////////////////////////////////////////////////////////////////////////////////

void Em8RunAction::bookHisto()
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
    histo2 = hbookManager->histogram("Energy Loss (keV)"
                                     ,nbinEn,Enlow/keV,Enhigh/keV) ;
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
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
  { 
      HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
  }  
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)    UI->ApplyCommand("/vis/scene/notifyHandlers");

      
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

//  plot definitions 
   
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

    for (G4int ien=0; ien<nbinEn; ien++)   distEn[ien]=0.;
  }
  if(nbinTt>0)
  {
    dTt = (Tthigh-Ttlow)/nbinTt ;
    entryTt=0.;
    underTt=0.;
    overTt=0.;

    for (G4int itt=0; itt<nbinTt; itt++)  distTt[itt]=0.;

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

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::EndOfRunAction(const G4Run* aRun)
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

  if(nbinStep>0)
  {G4double E , dnorm, norm ;
   G4cout << "   step number/event distribution " << G4endl ;
   G4cout << "#entries=" << entryStep << "    #underflows=" << underStep <<
             "    #overflows=" << overStep << G4endl ;
   if( entryStep>0.)
   {
     E = Steplow - dStep ;
     norm = TotNbofEvents ;
     G4cout << " bin nb   nsteplow     entries     normalized " << G4endl ;
     for(G4int iss=0; iss<nbinStep; iss++)
     {
      E += dStep ;
      dnorm = distStep[iss]/norm;
      G4cout << G4std::setw(5) << iss << G4std::setw(10) << E << 
                G4std::setw(12) << distStep[iss] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }     
  }
  if(nbinEn > 0)
  {
    G4std::ofstream fileOut("distribution.out", G4std::ios::out ) ;
    fileOut.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4std::ofstream normOut("normDist.out", G4std::ios::out ) ;
    normOut.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4double E , dnorm, norm,fmax,Emp,width ;
    Emp=-999.999 ;
    G4cout << " energy deposit distribution " << G4endl ;
    G4cout << "#entries=" << entryEn << "    #underflows=" << underEn <<
             "    #overflows=" << overEn << G4endl ;
    if( entryEn>0.)
    {
      E = Enlow - dEn ;
      norm = TotNbofEvents*1.0 ;   // *dEn ;
      G4cout << " bin nb      Elow      entries     normalized " << G4endl ;
      fmax = 0. ;

      for(G4int ien=0; ien<nbinEn; ien++)
      {
        E += dEn ;

        if(distEn[ien]>fmax)
        {
          fmax = distEn[ien] ;
          Emp = E ;                // most probable roughly
        }
        dnorm = distEn[ien]/norm;

        G4cout << G4std::setw(5) << ien << G4std::setw(10) << E/keV << 
                  G4std::setw(12) << distEn[ien] <<
                  G4std::setw(12) << dnorm << G4endl ;

        fileOut << E/keV << "\t"<< distEn[ien] << G4endl ;
        normOut << E/keV << "\t"<< dnorm << G4endl ;
      }
      G4cout << G4endl;
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

        if(distEn[i2] >= 0.5*fmax)   E2=E ;
      }
      G4cout << " Emp = " << G4std::setw(15) << Emp/MeV << "   width="
            << G4std::setw(15) << (E2-E1)/MeV <<   "  MeV " << G4endl;
      G4cout << G4endl ;
    }     
  }
  if(nbinTt>0)
  {
     G4double E , dnorm, norm ,sig;
   G4cout << " transmitted energy distribution " << G4endl ;
   G4cout << "#entries=" << entryTt << "    #underflows=" << underTt <<
             "    #overflows=" << overTt << G4endl ;
   if( entryTt>0.)
   {
     Ttmean /= entryTt;
     sig=Tt2mean/entryTt-Ttmean*Ttmean ;
     if(sig<=0.)
       sig=0.;
     else
       sig=sqrt(sig/entryTt) ;
     G4cout << " mean energy of transmitted particles=" << Ttmean/keV << 
               " +- " << sig/keV << "  keV." << G4endl;
     E = Ttlow - dTt ;
     norm = TotNbofEvents*dTt ;
     G4cout << " bin nb      Elow      entries     normalized " << G4endl ;
     for(G4int itt=0; itt<nbinTt; itt++)
     {
      E += dTt ;
      dnorm = distTt[itt]/norm;
      G4cout << G4std::setw(5) << itt << G4std::setw(10) << E << 
                G4std::setw(12) << distTt[itt] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }     
  }
  if(nbinTb>0)
  {
     G4double E , dnorm, norm ,sig;
   G4cout << " backscattered energy distribution " << G4endl ;
   G4cout << "#entries=" << entryTb << "    #underflows=" << underTb <<
             "    #overflows=" << overTb << G4endl ;
   if( entryTb>0.)
   {
     Tbmean /= entryTb;
     sig=Tb2mean/entryTb-Tbmean*Tbmean ;
     if(sig<=0.)
       sig=0.;
     else
       sig=sqrt(sig/entryTb) ;
     G4cout << " mean energy of backscattered particles=" << Tbmean/keV << 
               " +- " << sig/keV << "  keV." << G4endl;
     E = Tblow - dTb ;
     norm = TotNbofEvents*dTb ;
     G4cout << " bin nb      Elow      entries     normalized " << G4endl ;
     for(G4int itt=0; itt<nbinTb; itt++)
     {
      E += dTb ;
      dnorm = distTb[itt]/norm;
      G4cout << G4std::setw(5) << itt << G4std::setw(10) << E << 
                G4std::setw(12) << distTb[itt] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }     
  }
  if(nbinTsec>0)
  {G4double E , dnorm, norm ;
   G4cout << " energy distribution of charged secondaries " << G4endl ;
   G4cout << "#entries=" << entryTsec << "    #underflows=" << underTsec <<
             "    #overflows=" << overTsec << G4endl ;
   if( entryTsec>0.)
   {
     E = Tseclow - dTsec ;
     norm = TotNbofEvents*dTsec ;
     G4cout << " bin nb      Elow      entries     normalized " << G4endl ;
     for(G4int itt=0; itt<nbinTsec; itt++)
     {
      E += dTsec ;
      dnorm = distTsec[itt]/norm;
      G4cout << G4std::setw(5) << itt << G4std::setw(10) << E << 
                G4std::setw(12) << distTsec[itt] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }     
  }

  if(nbinR >0)
  {G4double R , dnorm, norm,sig  ;
   G4cout << "  R  distribution " << G4endl ;
   G4cout << "#entries=" << entryR  << "    #underflows=" << underR  <<
             "    #overflows=" << overR  << G4endl ;
   if( entryR >0.)
   {
     Rmean /= entryR;
     sig = R2mean/entryR - Rmean*Rmean;
     if(sig<=0.) sig=0. ;
     else        sig = sqrt(sig/entryR) ;
     G4cout << " mean lateral displacement at exit=" << Rmean/mm << " +- "
            << sig/mm << "  mm." << G4endl ; 
     R = Rlow - dR  ;
     norm = TotNbofEvents*dR  ;
     G4cout << " bin nb      Rlow      entries     normalized " << G4endl ;
     for(G4int ier=0; ier<nbinR ; ier++)
     {
      R+= dR  ;
      dnorm = distR[ier]/norm;
      G4cout << G4std::setw(5) << ier << G4std::setw(10) << R  <<
                G4std::setw(12) << distR[ier] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }
  }

  if(nbinTh>0)
  {G4double Th,Thdeg, dnorm, norm,fac0,fnorm,pere,Thpere,Thmean,sum;
   G4cout << "      angle   distribution " << G4endl ;
   G4cout << "#entries=" << entryTh << "    #underflows=" << underTh <<
             "    #overflows=" << overTh << G4endl ;
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
               " entries         normalized " << G4endl ;
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
      G4cout << G4std::setw(5) << ien << G4std::setw(10) << Thdeg << "   " <<
                G4std::setw(10) << Th << "  " <<   
                G4std::setw(12) << distTh[ien] << "  " <<
                G4std::setw(12) << dnorm << "  " << G4std::setw(12) << fnorm <<G4endl ;
     }
     Thmean /= sum ;
     G4cout << G4endl;
     G4cout << " mean = " << Thmean << "  rad  or " << 180.*Thmean/pi <<
               " deg." << G4endl;
     G4cout << " theta(1/e)=" << Thpere << " - " << Thpere+dTh << " rad   "
            << " or " << 180.*Thpere/pi << " - " << 180.*(Thpere+dTh)/pi 
            << " deg." << G4endl;
     G4cout << G4endl;
   }
  }

  if(nbinThback>0)
  {G4double Thb,Thdegb, dnormb, normb,fac0b,fnormb,pereb,Thpereb,Thmeanb,sumb;
   G4cout << " backscattering angle   distribution " << G4endl ;
   G4cout << "#entries=" << entryThback << "    #underflows=" << underThback <<
             "    #overflows=" << overThback << G4endl ;
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
               " entries         normalized " << G4endl ;
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
      G4cout << G4std::setw(5) << ien << G4std::setw(10) << Thdegb << "   " <<
                G4std::setw(10) << Thb << "  " <<
                G4std::setw(12) << distThback[ien] << "  " <<
                G4std::setw(12) << dnormb << "  " << G4std::setw(12) << fnormb <<G4endl ;
     }
     Thmeanb /= sumb ;
     G4cout << G4endl;
     G4cout << " mean = " << Thmeanb << "  rad  or " << 180.*Thmeanb/pi <<
               " deg." << G4endl;
     G4cout << " theta(1/e)=" << Thpereb << " - " << Thpereb+dThback << " rad   "
            << " or " << 180.*Thpereb/pi << " - " << 180.*(Thpereb+dThback)/pi
            << " deg." << G4endl;
     G4cout << G4endl;
   }
  }

  if(nbinGamma>0)
  {G4double E , fact,dnorm, norm  ;
   G4cout << " gamma energy distribution " << G4endl ;
   G4cout << "#entries=" << entryGamma << "    #underflows=" << underGamma <<
             "    #overflows=" << overGamma << G4endl ;
   if( entryGamma>0.)
   {
     fact=exp(dEGamma) ;
     E = ElowGamma/fact  ;
     norm = TotNbofEvents*dEGamma;
     G4cout << " bin nb         Elow      entries       normalized " << G4endl ;
     for(G4int itt=0; itt<nbinGamma; itt++)
     {
      E *= fact ;
      dnorm = distGamma[itt]/norm;
      G4cout << G4std::setw(5) << itt << G4std::setw(13) << E << 
                G4std::setw(12) << distGamma[itt] <<
                G4std::setw(15) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }     
  }

  if(nbinvertexz >0)
  {G4double z , dnorm, norm  ;
   G4cout << " vertex Z  distribution " << G4endl ;
   G4cout << "#entries=" << entryvertexz  << "    #underflows=" << undervertexz  <<
             "    #overflows=" << oververtexz  << G4endl ;
   if( entryvertexz >0.)
   {
     z =zlow - dz  ;
     norm = TotNbofEvents*dz  ;
     G4cout << " bin nb      zlow      entries     normalized " << G4endl ;
     for(G4int iez=0; iez<nbinvertexz ; iez++)
     {
      z+= dz  ;
      if(abs(z)<1.e-12) z=0.;
      dnorm = distvertexz[iez]/norm;
      G4cout << G4std::setw(5) << iez << G4std::setw(10) << z  <<
                G4std::setw(12) << distvertexz[iez] <<
                G4std::setw(12) << dnorm << G4endl ;
     }
     G4cout << G4endl;
   }
  }
  
 G4cout.precision(prec);
  
  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
  // Write histogram file

  hbookManager->write();
   
  // save Rndm status

  if (saveRndm == 1)
  { 
    HepRandom::showEngineStatus();
    HepRandom::saveEngineStatus("endOfRun.rndm");
  }     
}

///////////////////////////////////////////////////////////////////////////

void Em8RunAction::CountEvent()
{
  TotNbofEvents += 1. ;
}

/////////////////////////////////////////////////////////////////////////

void Em8RunAction::AddnStepsCharged(G4double ns)
{
  nStepSumCharged += ns;
  nStepSum2Charged += ns*ns;
}

////////////////////////////////////////////////////////////////////////

void Em8RunAction::AddnStepsNeutral(G4double ns)
{
  nStepSumNeutral += ns;
  nStepSum2Neutral += ns*ns;
}

////////////////////////////////////////////////////////////////////////////

void Em8RunAction::AddEdeps(G4double Eabs)
{
  EnergySumAbs += Eabs;
  EnergySquareSumAbs += Eabs*Eabs;
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::AddTrackLength(G4double tlabs)
{
  tlSumAbs += tlabs;
  tlsquareSumAbs += tlabs*tlabs ;
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::AddTrRef(G4double tr,G4double ref)
{
  Transmitted += tr ;
  Reflected   += ref;
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillNbOfSteps(G4double ns)
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
      ibin= (G4int)bin ;
      distStep[ibin] += 1. ;
    }
   histo1->accumulate(ns) ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillEn(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo2)
  {
    entryEn += 1. ;
 
    if(En < Enlow)          underEn += 1. ;
    else if( En >= Enhigh)  overEn  += 1. ;
    else
    {
      bin = (En-Enlow)/dEn ;
      ibin= (G4int)bin ;
      distEn[ibin] += 1. ;
    }
    histo2->accumulate(En/keV) ; // was /MeV
  }
}

////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillTt(G4double En)
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
      ibin= (G4int)bin ;
      distTt[ibin] += 1. ;
    }
  histo5->accumulate(En/MeV) ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillTb(G4double En)
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
      ibin= (G4int)bin ;
      distTb[ibin] += 1. ;
    }
  histo7->accumulate(En/MeV) ;
  }
}

///////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillTsec(G4double En)
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
      ibin= (G4int)bin ;
      distTsec[ibin] += 1. ;
    }
  histo8->accumulate(En/MeV) ;
  }
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillGammaSpectrum(G4double En)
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
      ibin= (G4int)bin ;
      distGamma[ibin] += 1. ;
    }
  histo10->accumulate(log10(En/MeV)) ;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillTh(G4double Th)
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
      ibin= (G4int)bin ;
      Thbin = Thlow+ibin*dTh ;
      if(Th > 0.001*dTh)
        wg=cn/sin(Th) ;
      else
      {  
        G4double thdeg=Th*180./pi;
        G4cout << "theta < 0.001*dth (from plot excluded) theta="
               << G4std::setw(12) << setprecision(4) << thdeg << G4endl;
        wg=0. ; 
      }
      distTh[ibin] += wg  ;
    }

  histo3->accumulate(Th/deg, wg) ;
  }
}

//////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillThBack(G4double Th)
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
      ibin= (G4int)bin ;
      Thbin = Thlowback+ibin*dThback ;
      if(Th > 0.001*dThback)
        wg=cn/sin(Th) ;
      else
      {  
        G4double thdeg=Th*180./pi;
        G4cout << "theta < 0.001*dth (from plot excluded) theta="
               << G4std::setw(12) << setprecision(4) << thdeg << G4endl;
        wg=0. ; 
      }
      distThback[ibin] += wg  ;
    }
  histo6->accumulate(Th/deg, wg) ;
  }

}

//////////////////////////////////////////////////////////////////////////

void Em8RunAction::FillR(G4double R )
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
      ibin= (G4int)bin ;
      distR[ibin] += 1. ;
    }
  histo4->accumulate(R/mm) ;
  }
}

/////////////////////////////////////////////////////////////////////////////

void Em8RunAction::Fillvertexz(G4double z )
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
      ibin = (G4int)bin ;
      distvertexz[ibin] += 1. ;
    }
  histo9->accumulate(z/mm) ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Em8RunAction::SethistName(G4String name)
{
  histName = name ;
  G4cout << " hist file = " << histName << G4endl;
}

void Em8RunAction::SetnbinStep(G4int nbin)
{
  nbinStep = nbin ;
  if(nbinStep>0)
  G4cout << " Nb of bins in #step plot = " << nbinStep << G4endl ;
}
void Em8RunAction::SetSteplow(G4double low)
{
  Steplow = low ;
  if(nbinStep>0)
  G4cout << " low  in the #step plot = " << Steplow << G4endl ;
}
void Em8RunAction::SetStephigh(G4double high)
{
  Stephigh = high ;
  if(nbinStep>0)
  G4cout << " high in the #step plot = " << Stephigh << G4endl ;
}

////////////////////////////////////////////////////////////////////////

void Em8RunAction::SetnbinEn(G4int nbin)
{
  nbinEn = nbin ;

  if(nbinEn > 0) G4cout << " Nb of bins in Edep plot = " << nbinEn << G4endl ;
}
void Em8RunAction::SetEnlow(G4double Elow)
{
  Enlow = Elow ;
  if(nbinEn>0) G4cout << " Elow  in the  Edep plot = " << Enlow << G4endl ;
}
void Em8RunAction::SetEnhigh(G4double Ehigh)
{
  Enhigh = Ehigh ;

  if(nbinEn>0) G4cout << " Ehigh in the  Edep plot = " << Enhigh << G4endl ;
}

/////////////////////////////////////////////////////////////////////////

void Em8RunAction::SetnbinGamma(G4int nbin)
{
  nbinGamma = nbin ;
  if(nbinGamma>0)
  G4cout << " Nb of bins in gamma spectrum plot = " << nbinGamma << G4endl ;
}
void Em8RunAction::SetElowGamma(G4double Elow)
{
  ElowGamma = Elow ;
  if(nbinGamma>0)
  G4cout << " Elow  in the gamma spectrum plot = " << ElowGamma << G4endl ;
}
void Em8RunAction::SetEhighGamma(G4double Ehigh)
{
  EhighGamma = Ehigh ;
  if(nbinGamma>0)
  G4cout << " Ehigh in the gamma spectrum plot = " << EhighGamma << G4endl ;
}

void Em8RunAction::SetnbinTt(G4int nbin)
{
  nbinTt = nbin ;
  if(nbinTt>0)
  G4cout << " Nb of bins in Etransmisssion plot = " << nbinTt << G4endl ;
}
void Em8RunAction::SetTtlow(G4double Elow)
{
  Ttlow = Elow ;
  if(nbinTt>0)
  G4cout << " Elow  in the  Etransmission plot = " << Ttlow << G4endl ;
}
void Em8RunAction::SetTthigh(G4double Ehigh)
{
  Tthigh = Ehigh ;
  if(nbinTt>0)
  G4cout << " Ehigh in the  Etransmission plot = " << Tthigh << G4endl ;
}
void Em8RunAction::SetnbinTb(G4int nbin)
{
  nbinTb = nbin ;
  if(nbinTb>0)
  G4cout << " Nb of bins in Ebackscattered plot = " << nbinTb << G4endl ;
}
void Em8RunAction::SetTblow(G4double Elow)
{
  Tblow = Elow ;
  if(nbinTb>0)
  G4cout << " Elow  in the  Ebackscattered plot = " << Tblow << G4endl ;
}
void Em8RunAction::SetTbhigh(G4double Ehigh)
{
  Tbhigh = Ehigh ;
  if(nbinTb>0)
  G4cout << " Ehigh in the  Ebackscattered plot = " << Tbhigh << G4endl ;
}

void Em8RunAction::SetnbinTsec(G4int nbin)
{
  nbinTsec = nbin ;
  if(nbinTsec>0)
  G4cout << " Nb of bins in Tsecondary  plot = " << nbinTsec << G4endl ;
}
void Em8RunAction::SetTseclow(G4double Elow)
{
  Tseclow = Elow ;
  if(nbinTsec>0)
  G4cout << " Elow  in the  Tsecondary plot = " << Tseclow << G4endl ;
}
void Em8RunAction::SetTsechigh(G4double Ehigh)
{
  Tsechigh = Ehigh ;
  if(nbinTsec>0)
  G4cout << " Ehigh in the  Tsecondary plot = " << Tsechigh << G4endl ;
}
 
void Em8RunAction::SetnbinR(G4int nbin)
{
  nbinR  = nbin ;
  if(nbinR>0)
  G4cout << " Nb of bins in R plot = " << nbinR << G4endl ;
}
void Em8RunAction::SetRlow(G4double rlow)
{
  Rlow = rlow ;
  if(nbinR>0)
  G4cout << " Rlow  in the  R plot = " << Rlow << G4endl ;
}
void Em8RunAction::SetRhigh(G4double rhigh)
{
  Rhigh = rhigh ;
  if(nbinR>0)
  G4cout << " Rhigh in the R plot = " << Rhigh << G4endl ;
}

void Em8RunAction::Setnbinzvertex(G4int nbin)
{
  nbinvertexz  = nbin ;
  if(nbinvertexz>0)
  G4cout << " Nb of bins in Z plot = " << nbinvertexz << G4endl ;
}
void Em8RunAction::Setzlow(G4double z)
{
  zlow = z ;
  if(nbinvertexz>0)
  G4cout << " zlow  in the  Z plot = " << zlow << G4endl ;
}
void Em8RunAction::Setzhigh(G4double z)
{
  zhigh = z ;
  if(nbinvertexz>0)
  G4cout << " zhigh in the Z plot = " << zhigh << G4endl ;
}

void Em8RunAction::SetnbinTh(G4int nbin)
{
  nbinTh = nbin ;
  if(nbinTh>0)
  G4cout << " Nb of bins in Theta plot = " << nbinTh << G4endl ;
}
void Em8RunAction::SetThlow(G4double Tlow)
{
  Thlow = Tlow ;
  if(nbinTh>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlow << G4endl ;
}
void Em8RunAction::SetThhigh(G4double Thigh)
{
  Thhigh = Thigh ;
  if(nbinTh>0)
  G4cout << " Thigh in the Theta plot = " << Thhigh << G4endl ;
}

void Em8RunAction::SetnbinThBack(G4int nbin)
{
  nbinThback = nbin ;
  if(nbinThback>0)
  G4cout << " Nb of bins in Theta plot = " << nbinThback << G4endl ;
}
void Em8RunAction::SetThlowBack(G4double Tlow)
{
  Thlowback = Tlow ;
  if(nbinThback>0)
  G4cout << " Tlow  in the  Theta plot = " << Thlowback << G4endl ;
}
void Em8RunAction::SetThhighBack(G4double Thigh)
{
  Thhighback = Thigh ;
  if(nbinThback>0)
  G4cout << " Thigh in the Theta plot = " << Thhighback << G4endl ;
}
void Em8RunAction::CountParticles(G4double nch,G4double nne)
{
  SumCharged += nch ;
  SumNeutral += nne ;
  Sum2Charged += nch*nch ;
  Sum2Neutral += nne*nne ;
}
void Em8RunAction::AddEP(G4double nele,G4double npos) 
{
  Selectron += nele;
  Spositron += npos;
}

//
//
////////////////////////////////////////////////////////////////////////




