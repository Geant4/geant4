// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17RunAction.cc,v 1.1 1999-11-30 18:01:56 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst17RunAction.hh"
#include "Tst17RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <iomanip.h>
#include "G4VVisManager.hh"

#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17RunAction::Tst17RunAction()
  :histName("histfile"),nbinStep(0),nbinEn(0),nbinTt(0),nbinTb(0),
   nbinTsec(0),nbinTh(0),nbinThback(0),nbinR(0),nbinGamma(0),
   nbinvertexz(0),histo1(0),histo2(0),histo3(0),histo4(0)
{
  runMessenger = new Tst17RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17RunAction::~Tst17RunAction()
{
  delete runMessenger;

  if(histo1) delete histo1 ;
  if(histo2) delete histo2 ;
  if(histo3) delete histo3 ;
  if(histo4) delete histo4 ;
  delete hbookManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Tst17RunAction::bookHisto()
{
  // init hbook
  hbookManager = new HBookFile(histName, 68);
  assert (hbookManager != 0);

  // book histograms
  if(nbinEn>0)
  {
    histo1 = hbookManager->histogram("energy deposit in absorber(in MeV)"
                                     ,nbinEn,Enlow,Enhigh) ;
    assert (histo1 != 0);
  }

  if(nbinTt>0)
  {
    histo2 = hbookManager->histogram("kinetic energy of the primary at exit(MeV)"
                                     ,nbinTt,Ttlow,Tthigh)  ;
    assert (histo2 != 0);
  }

  if(nbinTsec>0)
  {
    histo3 = hbookManager->histogram("kinetic energy of the charged secondaries (MeV)"
                                     ,nbinTsec,Tseclow,Tsechigh)  ;
    assert (histo3 != 0);
  }

  if(nbinGamma>0)
  {
    histo4 = hbookManager->histogram("kinetic energy of gammas escaping the absorber (MeV)"
                                ,nbinGamma,log10(ElowGamma),log10(EhighGamma))  ;
    assert (histo4 != 0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;

  G4UImanager* UI = G4UImanager::GetUIpointer();
   
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

  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  if(pVVisManager){

    UI->ApplyCommand("/vis/clear/view");
    UI->ApplyCommand("/vis/draw/current");
  } 

  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17RunAction::EndOfRunAction(const G4Run* aRun)
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
 G4cout.precision(6);
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

  if(nbinEn>0){

    G4double E , dnorm, norm,fmax,Emp,width ;
    Emp=-999.999 ;
    G4cout <<" energy deposit distribution " << endl ;
    G4cout <<"#entries=" << entryEn << "    #underflows=" << underEn 
	   <<"    #overflows=" << overEn << endl ;

   if( entryEn>0.){

     E = Enlow - dEn ;
     norm = TotNbofEvents*dEn ;
     G4cout << " bin nb      Elow      entries     normalized " << endl ;
     fmax = 0. ;
     for(G4int ien=0; ien<nbinEn; ien++){

      E += dEn ;
      if(distEn[ien]>fmax){

        fmax=distEn[ien];
        Emp = E ;
      }
      dnorm = distEn[ien]/norm;
      G4cout << setw(5) << ien << setw(10) << E 
	     << setw(12) << distEn[ien] 
	     << setw(12) << dnorm << endl ;
     }
     
     G4cout << endl;
     G4int ii ;
     G4double E1,E2 ;
     E1=-1.e6 ;
     E2=+1.e6 ;
     E = Enlow -dEn ;
     ii = -1;
     for(G4int i1=0; i1<nbinEn; i1++){
       
       E += dEn ;
       if(ii<0){
	 
	 if(distEn[i1] >= 0.5*fmax){
	   
	   E1=E ;
	   ii=i1 ;
	 }
       }
     }

     E = Enlow -dEn ;
     for(G4int i2=0; i2<nbinEn; i2++){

       E += dEn ;
       if(distEn[i2] >= 0.5*fmax)
	 E2=E ;
     }

     G4cout << " Emp = " << setw(15) << Emp/MeV << "   width="
            << setw(15) << (E2-E1)/MeV <<   "  MeV " << endl;
     G4cout << endl ;
   }     
  }

  if(nbinTt>0){

    G4double E , dnorm, norm ,sig;
    G4cout << " transmitted energy distribution " << endl ;
    G4cout << "#entries=" << entryTt << "    #underflows=" << underTt 
	   << "    #overflows=" << overTt << endl ;

    if( entryTt>0.){

      Ttmean /= entryTt;
      sig=Tt2mean/entryTt-Ttmean*Ttmean ;
      if(sig<=0.)
	sig=0.;
      else
	sig=sqrt(sig/entryTt) ;
      G4cout << " mean energy of transmitted particles=" << Ttmean/keV 
	     << " +- " << sig/keV << "  keV." << endl;

      E = Ttlow - dTt ;
      norm = TotNbofEvents*dTt ;

      G4cout << " bin nb      Elow      entries     normalized " << endl ;

     for(G4int itt=0; itt<nbinTt; itt++){
       
       E += dTt ;
       dnorm = distTt[itt]/norm;
       G4cout << setw(5) << itt << setw(10) << E 
	      << setw(12) << distTt[itt] 
	      << setw(12) << dnorm << endl ;
     }
     
     G4cout << endl;
    }     
  }

  if(nbinTb>0){

    G4double E , dnorm, norm ,sig;
    G4cout << " backscattered energy distribution " << endl ;
    G4cout << "#entries=" << entryTb << "    #underflows=" << underTb 
	   << "    #overflows=" << overTb << endl ;

   if( entryTb>0.){

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

   // Write histogram file
  hbookManager->write();

  if (G4VVisManager::GetConcreteInstance()){

    G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17RunAction::CountEvent()
{
  TotNbofEvents += 1. ;
}

void Tst17RunAction::AddnStepsCharged(G4double ns)
{
  nStepSumCharged += ns;
  nStepSum2Charged += ns*ns;
}
void Tst17RunAction::AddnStepsNeutral(G4double ns)
{
  nStepSumNeutral += ns;
  nStepSum2Neutral += ns*ns;
}
void Tst17RunAction::AddEdeps(G4double Eabs)
{
  EnergySumAbs += Eabs;
  EnergySquareSumAbs += Eabs*Eabs;
}
void Tst17RunAction::AddTrackLength(G4double tlabs)
{
  tlSumAbs += tlabs;
  tlsquareSumAbs += tlabs*tlabs ;

}
void Tst17RunAction::AddTrRef(G4double tr,G4double ref)
{
  Transmitted += tr ;
  Reflected   += ref;
}

void Tst17RunAction::FillEn(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo1)
  {
    entryEn += 1. ;
 
    if(En<Enlow)
      underEn += 1. ;
    else if(En>=Enhigh)
      overEn  += 1. ;
    else
    {
      bin = (En-Enlow)/dEn ;
      ibin = (G4int) bin ;
      distEn[ibin] += 1. ;
    }
    histo1->accumulate(En/MeV) ;
  }
}

void Tst17RunAction::FillTt(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo2)
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
      ibin = (G4int) bin ;
      distTt[ibin] += 1. ;
    }
    histo2->accumulate(En/MeV) ;
  }
}

void Tst17RunAction::FillTsec(G4double En)
{
  G4double bin ;
  G4int ibin;
  if(histo3)
  {
    entryTsec += 1. ;

    if(En<Tseclow)
      underTsec += 1. ;
    else if(En>=Tsechigh)
      overTsec  += 1. ;
    else
    {
      bin = (En-Tseclow)/dTsec ;
      ibin = (G4int) bin ;
      distTsec[ibin] += 1. ;
    }
    histo3->accumulate(En/MeV);
  }
}

void Tst17RunAction::FillGammaSpectrum(G4double En)
{
  G4double bin ;
  G4int ibin;

  if(histo4)
  {
    entryGamma += 1. ;

    if(En<ElowGamma)
      underGamma += 1. ;
    else if(En>=EhighGamma)
      overGamma  += 1. ;
    else
    {
      bin = log(En/ElowGamma)/dEGamma;
      ibin = (G4int) bin ;
      distGamma[ibin] += 1. ;
    }
    histo4->accumulate(log10(En/MeV)) ;
  }
}

void Tst17RunAction::SethistName(G4String name)
{
  histName = name ;
  G4cout << " hist file = " << histName << endl;
}

void Tst17RunAction::SetnbinEn(G4int nbin)
{
  nbinEn = nbin ;
  if(nbinEn>0)
  G4cout << " Nb of bins in Edep plot = " << nbinEn << endl ;
}
void Tst17RunAction::SetEnlow(G4double Elow)
{
  Enlow = Elow ;
  if(nbinEn>0)
  G4cout << " Elow  in the  Edep plot = " << Enlow << endl ;
}
void Tst17RunAction::SetEnhigh(G4double Ehigh)
{
  Enhigh = Ehigh ;
  if(nbinEn>0)
  G4cout << " Ehigh in the  Edep plot = " << Enhigh << endl ;
}

void Tst17RunAction::SetnbinGamma(G4int nbin)
{
  nbinGamma = nbin ;
  if(nbinGamma>0)
  G4cout << " Nb of bins in gamma spectrum plot = " << nbinGamma << endl ;
}
void Tst17RunAction::SetElowGamma(G4double Elow)
{
  ElowGamma = Elow ;
  if(nbinGamma>0)
  G4cout << " Elow  in the gamma spectrum plot = " << ElowGamma << endl ;
}
void Tst17RunAction::SetEhighGamma(G4double Ehigh)
{
  EhighGamma = Ehigh ;
  if(nbinGamma>0)
  G4cout << " Ehigh in the gamma spectrum plot = " << EhighGamma << endl ;
}

void Tst17RunAction::SetnbinTt(G4int nbin)
{
  nbinTt = nbin ;
  if(nbinTt>0)
  G4cout << " Nb of bins in Etransmisssion plot = " << nbinTt << endl ;
}
void Tst17RunAction::SetTtlow(G4double Elow)
{
  Ttlow = Elow ;
  if(nbinTt>0)
  G4cout << " Elow  in the  Etransmission plot = " << Ttlow << endl ;
}
void Tst17RunAction::SetTthigh(G4double Ehigh)
{
  Tthigh = Ehigh ;
  if(nbinTt>0)
  G4cout << " Ehigh in the  Etransmission plot = " << Tthigh << endl ;
}

void Tst17RunAction::SetnbinTsec(G4int nbin)
{
  nbinTsec = nbin ;
  if(nbinTsec>0)
  G4cout << " Nb of bins in Tsecondary  plot = " << nbinTsec << endl ;
}

void Tst17RunAction::SetTseclow(G4double Elow)
{
  Tseclow = Elow ;
  if(nbinTsec>0)
  G4cout << " Elow  in the  Tsecondary plot = " << Tseclow << endl ;
}
void Tst17RunAction::SetTsechigh(G4double Ehigh)
{
  Tsechigh = Ehigh ;
  if(nbinTsec>0)
  G4cout << " Ehigh in the  Tsecondary plot = " << Tsechigh << endl ;
}
 
void Tst17RunAction::CountParticles(G4double nch,G4double nne)
{
  SumCharged += nch ;
  SumNeutral += nne ;
  Sum2Charged += nch*nch ;
  Sum2Neutral += nne*nne ;
}
void Tst17RunAction::AddEP(G4double nele,G4double npos) 
{
  Selectron += nele;
  Spositron += npos;
}


