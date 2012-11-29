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
// $Id$
//
//
// G4 Physics class: G4QFreeScattering for quasi-free scattering
// Created: M.V. Kossov, CERN/ITEP(Moscow), 26-OCT-11
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 29-Oct-11
// 
// ----------------------------------------------------------------------
// Short description: Provides quasi-free scattering on nuclear nucleons.
// ----------------------------------------------------------------------

//#define debug
//#define pdebug
//#define ppdebug
//#define nandebug

#include "G4QFreeScattering.hh"
#include "G4SystemOfUnits.hh"

// initialisation of statics
std::vector<std::pair<G4double,G4double>*> G4QFreeScattering::vX; // ETPointers to LogTable

G4QFreeScattering::G4QFreeScattering()
{
#ifdef pdebug
  G4cout<<"***^^^*** G4QFreeScattering singletone is created ***^^^***"<<G4endl;
#endif
}

G4QFreeScattering::~G4QFreeScattering()
{
  std::vector<std::pair<G4double,G4double>*>::iterator pos2;
  for(pos2=vX.begin(); pos2<vX.end(); pos2++) delete [] * pos2;
  vX.clear();
}

// Returns Pointer to the G4VQCrossSection class
G4QFreeScattering* G4QFreeScattering::GetPointer()
{
  static G4QFreeScattering theQFS;   // *** Static body of the Quasi-Free Scattering ***
  return &theQFS;
}

// Calculatio pair(hN_el,hN_tot) (mb): p in GeV/c, index(PDG,F) (see FetchElTot)
std::pair<G4double,G4double> G4QFreeScattering::CalcElTot(G4double p, G4int I)
{
  // ---------> Each parameter set can have not more than nPoints=128 parameters
  static const G4double lmi=3.5;       // min of (lnP-lmi)^2 parabola
  static const G4double pbe=.0557;     // elastic (lnP-lmi)^2 parabola coefficient
  static const G4double pbt=.3;        // total (lnP-lmi)^2 parabola coefficient
  static const G4double pmi=.1;        // Below that fast LE calculation is made
  static const G4double pma=1000.;     // Above that fast HE calculation is made
  G4double El=0.;                      // prototype of the elastic hN cross-section
  G4double To=0.;                      // prototype of the total hN cross-section
  if(p<=0.)
  {
    //G4cout<<"-Warning-G4QFreeScattering::CalcElTot: p="<<p<<" is 0 or negative"<<G4endl;
    return std::make_pair(El,To);
  }
  if     (!I)                          // pp/nn
  {
#ifdef debug
    G4cout<<"G4QFreeScatter::CalcElTot:I=0, p="<<p<<", pmi="<<pmi<<", pma="<<pma<<G4endl;
#endif
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.00012+p2*.2);
      To=El;
#ifdef debug
      G4cout<<"G4QFreeScatter::CalcElTot:I=0i, El="<<El<<", To="<<To<<", p2="<<p2<<G4endl;
#endif
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
#ifdef debug
      G4cout<<"G4QFreeScat::CalcElTot:I=0a, El="<<El<<", To="<<To<<", lp2="<<lp2<<G4endl;
#endif
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.00012+p2*.2);
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(pbe*lp2+6.72+32.6/p)/(1.+rp2/p);
      To=LE+(pbt*lp2+38.2+52.7*rp2)/(1.+2.72*rp2*rp2);
#ifdef debug
      G4cout<<"G4QFreeScat::CalcElTot:0,E="<<El<<",T="<<To<<",s="<<p2<<",l="<<lp2<<G4endl;
#endif
    }
  }
  else if(I==1)                        // np/pn
  {
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.00012+p2*(.051+.1*p2));
      To=El;
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.00012+p2*(.051+.1*p2));
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(pbe*lp2+6.72+30./p)/(1.+.49*rp2/p);
      To=LE+(pbt*lp2+38.2)/(1.+.54*rp2*rp2);
    }
  }
  else if(I==2)                        // pimp/pipn
  {
    G4double lp=std::log(p);
    if(p<pmi)
    {
      G4double lr=lp+1.27;
      El=1.53/(lr*lr+.0676);
      To=El*3;
    }
    else if(p>pma)
    {
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      El=pbe*ld2+2.4+7./sp;
      To=pbt*ld2+22.3+12./sp;
    }
    else
    {
      G4double lr=lp+1.27;                    // p1
      G4double LE=1.53/(lr*lr+.0676);         // p2, p3        
      G4double ld=lp-lmi;                     // p4 (lmi=3.5)
      G4double ld2=ld*ld;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      G4double lm=lp+.36;                     // p5
      G4double md=lm*lm+.04;                  // p6
      G4double lh=lp-.017;                    // p7
      G4double hd=lh*lh+.0025;                // p8
      El=LE+(pbe*ld2+2.4+7./sp)/(1.+.7/p4)+.6/md+.05/hd;//p9(pbe=.0557),p10,p11,p12,p13,p14
      To=LE*3+(pbt*ld2+22.3+12./sp)/(1.+.4/p4)+1./md+.06/hd;
    }
  }
  else if(I==3)                        // pipp/pimn
  {
    G4double lp=std::log(p);
    if(p<pmi)
    {
      G4double lr=lp+1.27;
      G4double lr2=lr*lr;
      El=13./(lr2+lr2*lr2+.0676);
      To=El;
    }
    else if(p>pma)
    {
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      El=pbe*ld2+2.4+6./sp;
      To=pbt*ld2+22.3+5./sp;
    }
    else
    {
      G4double lr=lp+1.27;                   // p1
      G4double lr2=lr*lr;
      G4double LE=13./(lr2+lr2*lr2+.0676);   // p2, p3
      G4double ld=lp-lmi;                    // p4 (lmi=3.5)
      G4double ld2=ld*ld;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      G4double lm=lp-.32;                    // p5
      G4double md=lm*lm+.0576;               // p6
      El=LE+(pbe*ld2+2.4+6./sp)/(1.+3./p4)+.7/md; // p7(pbe=.0557), p8, p9, p10, p11
      To=LE+(pbt*ld2+22.3+5./sp)/(1.+1./p4)+.8/md;
    }
  }
  else if(I==4)                        // Kmp/Kmn/K0p/K0n
  {

    if(p<pmi)
    {
      G4double psp=p*std::sqrt(p);
      El=5.2/psp;
      To=14./psp;
    }
    else if(p>pma)
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      El=pbe*ld2+2.23;
      To=pbt*ld2+19.5;
    }
    else
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      G4double psp=p*sp;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double lm=p-.39;
      G4double md=lm*lm+.000156;
      G4double lh=p-1.;
      G4double hd=lh*lh+.0156;
      El=5.2/psp+(pbe*ld2+2.23)/(1.-.7/sp+.075/p4)+.004/md+.15/hd;
      To=14./psp+(pbt*ld2+19.5)/(1.-.21/sp+.52/p4)+.006/md+.30/hd;
    }
  }
  else if(I==5)                        // Kpp/Kpn/aKp/aKn
  {
    if(p<pmi)
    {
      G4double lr=p-.38;
      G4double lm=p-1.;
      G4double md=lm*lm+.372;   
      El=.7/(lr*lr+.0676)+2./md;
      To=El+.6/md;
    }
    else if(p>pma)
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      El=pbe*ld2+2.23;
      To=pbt*ld2+19.5;
    }
    else
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      G4double lr=p-.38;
      G4double LE=.7/(lr*lr+.0676);
      G4double sp=std::sqrt(p);
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double lm=p-1.;
      G4double md=lm*lm+.372;
      El=LE+(pbe*ld2+2.23)/(1.-.7/sp+.1/p4)+2./md;
      To=LE+(pbt*ld2+19.5)/(1.+.46/sp+1.6/p4)+2.6/md;
    }
  }
  else if(I==6)                        // hyperon-N
  {
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.002+p2*(.12+p2));
      To=El;
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double sp=std::sqrt(p);
      El=(pbe*lp2+6.72)/(1.+2./sp);
      To=(pbt*lp2+38.2+900./sp)/(1.+27./sp);
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.002+p2*(.12+p2));
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      El=LE+(pbe*lp2+6.72+99./p2)/(1.+2./sp+2./p4);
      To=LE+(pbt*lp2+38.2+900./sp)/(1.+27./sp+3./p4);
    }
  }
  else if(I==7)                        // antibaryon-N
  {
    if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
    }
    else
    {
      G4double ye=std::pow(p,1.25);
      G4double yt=std::pow(p,.35);
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=80./(ye+1.)+pbe*lp2+6.72;
      To=(80./yt+.3)/yt+pbt*lp2+38.2;
    }
  }
  else
  {
    G4cout<<"*Error*G4QFreeScattering::CalcElTot:ind="<<I<<" is not defined (0-7)"<<G4endl;
    G4Exception("G4QFreeScattering::CalcElTot:","23",FatalException,"CHIPScrash");
  }
  if(El>To) El=To;
  return std::make_pair(El,To);
} // End of CalcElTot

// For hadron PDG with momentum Mom (GeV/c) on N (p/n) calculate <sig_el,sig_tot> pair (mb)
// *** Only for external use. For simulation it is better to use GetElTot(...) [AMDB] ***
std::pair<G4double,G4double> G4QFreeScattering::GetElTotXS(G4double p, G4int PDG, G4bool F)
{
  G4int ind=0;                                 // Prototype of the reaction index
  G4bool kfl=true;                             // Flag of K0/aK0 oscillation
  G4bool kf=false;
  if(PDG==90001000) PDG=2212;
  if(PDG==90000001) PDG=2112;
  if(PDG==91000000) PDG=3122;
  if(PDG==130 || PDG==310)
  {
    kf=true;
    if(G4UniformRand()>.5) kfl=false;
  }
  if      ( (PDG == 2212 && F) || (PDG == 2112 && !F) ) ind=0; // pp/nn
  else if ( (PDG == 2112 && F) || (PDG == 2212 && !F) ) ind=1; // np/pn
  else if ( (PDG == -211 && F) || (PDG == 211 && !F) ) ind=2; // pimp/pipn
  else if ( (PDG == 211 && F) || (PDG == -211 && !F) ) ind=3; // pipp/pimn
  else if ( PDG == -321 || PDG == -311 || (kf && !kfl) ) ind=4; // KmN/K0N
  else if ( PDG == 321 || PDG == 311 || (kf && kfl) ) ind=5; // KpN/aK0N
  else if ( PDG >  3000 && PDG <  3335) ind=6; // @@ for all hyperons - take Lambda
  else if ( PDG > -3335 && PDG < -2000) ind=7; // @@ for all anti-baryons (anti-p/anti-n)
  else {
    // VI changed Fatal to Warning, because this method cannot do better
    //    source of problem in classes which make this call
    ind = 0;
    G4cout<<"*Error*G4QFreeScattering::CalcElTotXS: PDG="<<PDG
          <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
    G4Exception("G4QFreeScattering::CalcElTotXS:","22",JustWarning,"CHIPS_crash");
  }
  return CalcElTot(p,ind); // This is a slow direct method, better use FetchElTot
}

// Calculatio pair(hN_el,hN_tot)(mb): p in GeV/c, F=true -> N=proton, F=false -> N=neutron
std::pair<G4double,G4double> G4QFreeScattering::FetchElTot(G4double p, G4int PDG, G4bool F)
{
  static const G4int    nlp=300;         // Number of steps in the S(lnp) logTable(5% step)
  static const G4int    mlp=nlp+1;       // Number of elements in the S(lnp) logTable
  static const G4double lpi=-5.;         // The min ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
  static const G4double lpa=10.;         // The max ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
  static const G4double mi=std::exp(lpi);// The min p of logTabEl(~ 6.7 MeV/c)
  static const G4double ma=std::exp(lpa);// The max p of logTabEl(~ 22. TeV)
  static const G4double dl=(lpa-lpi)/nlp;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  //static const G4double toler=.001;      // Relative Tolarence defining "theSameMomentum"
  static G4double lastP=0.;              // The last momentum for which XS was calculated
  static G4int    lastH=0;               // The last projPDG for which XS was calculated
  static G4bool   lastF=true;            // The last nucleon for which XS was calculated
  static std::pair<G4double,G4double> lastR(0.,0.); // The last result
  // Local Associative Data Base:
  static std::vector<G4int>     vI;      // Vector of index for which XS was calculated
  static std::vector<G4double>  vM;      // Vector of rel max ln(p) initialized in LogTable
  static std::vector<G4int>     vK;      // Vector of topBin number initialized in LogTable
  // Last values of the Associative Data Base:
  static G4int     lastI=0;              // The Last index for which XS was calculated
  static G4double  lastM=0.;             // The Last rel max ln(p) initialized in LogTable
  static G4int     lastK=0;              // The Last topBin number initialized in LogTable
  static std::pair<G4double,G4double>* lastX=0; // The Last ETPointers to LogTable in heap
  // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
  G4int nDB=vI.size();                   // A number of hadrons already initialized in AMDB
  if(PDG==90001000) PDG=2212;
  if(PDG==90000001) PDG=2112;
  if(PDG==91000000) PDG=3122;
#ifdef pdebug
  G4cout<<"G4QFreeScatter::FetchElTot:p="<<p<<",PDG="<<PDG<<",F="<<F<<",nDB="<<nDB<<G4endl;
#endif
  if(nDB && lastH==PDG && lastF==F && p>0. && p==lastP) return lastR;// VI don't use toler.
  //  if(nDB && lastH==PDG && lastF==F && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
  lastH=PDG;
  lastF=F;
  G4int ind=-1;                          // Prototipe of the index of the PDG/F combination
  // i=0: pp(nn), i=1: np(pn), i=2: pimp(pipn), i=3: pipp(pimn), i=4: Kmp(Kmn,K0n,K0p),
  // i=5: Kpp(Kpn,aK0n,aK0p), i=6: Hp(Hn), i=7: app(apn,ann,anp) 
  G4bool kfl=true;                             // Flag of K0/aK0 oscillation
  G4bool kf=false;
  if(PDG==130 || PDG==310)
  {
    kf=true;
    if(G4UniformRand()>.5) kfl=false;
  }
  if      ( (PDG == 2212 && F) || (PDG == 2112 && !F) ) ind=0; // pp/nn
  else if ( (PDG == 2112 && F) || (PDG == 2212 && !F) ) ind=1; // np/pn
  else if ( (PDG == -211 && F) || (PDG == 211 && !F) ) ind=2; // pimp/pipn
  else if ( (PDG == 211 && F) || (PDG == -211 && !F) ) ind=3; // pipp/pimn
  else if ( PDG == -321 || PDG == -311 || (kf && !kfl) ) ind=4; // KmN/K0N
  else if ( PDG == 321 || PDG == 311 || (kf && kfl) ) ind=5; // KpN/aK0N
  else if ( PDG >  3000 && PDG <  3335) ind=6; // @@ for all hyperons - take Lambda
  else if ( PDG > -3335 && PDG < -2000) ind=7; // @@ for all anti-baryons (anti-p/anti-n)
  else
  {
    // VI changed Fatal to Warning, because this method cannot do better
    //    source of problem in classes which make this call
    ind = 0;
    G4cout<<"*Error*G4QFreeScattering::FetchElTot: PDG="<<PDG
          <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
    G4Exception("G4QFreeScattering::FetchElTot:","22",JustWarning,"CHIPS problem");
  }
  if(nDB && lastI==ind && p>0. && p==lastP) return lastR;  // VI do not use toler
  //  if(nDB && lastI==ind && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
  if(p<=mi || p>=ma) return CalcElTot(p,ind);   // @@ Slow calculation ! (Warning?)
  G4bool found=false;
  G4int i=-1;
  if(nDB) for (i=0; i<nDB; ++i) if(ind==vI[i])  // Sirch for this index in AMDB
  {
    found=true;                                 // The index is found
    break;
  }
  G4double lp=std::log(p);
#ifdef pdebug
  G4cout<<"G4QFreeScat::FetchElTot: I="<<ind<<", i="<<i<<",fd="<<found<<",lp="<<lp<<G4endl;
#endif
  if(!nDB || !found)                            // Create new line in the AMDB
  {
#ifdef pdebug
    G4cout<<"G4QFreeScattering::FetchElTot: NewX, ind="<<ind<<", nDB="<<nDB<<G4endl;
#endif
    lastX = new std::pair<G4double,G4double>[mlp]; // Create logarithmic Table for ElTot
    lastI = ind;                                // Remember the initialized inex
    lastK = static_cast<int>((lp-lpi)/dl)+1;    // MaxBin to be initialized in LogTaB
    if(lastK>nlp)
    {
      lastK=nlp;
      lastM=lpa-lpi;
    }
    else lastM = lastK*dl;               // Calculate max initialized ln(p)-lpi for LogTab
    G4double pv=mi;
    for(G4int j=0; j<=lastK; ++j)        // Calculate LogTab values
    {
      lastX[j]=CalcElTot(pv,ind);
#ifdef pdebug
      G4cout<<"G4QFreeScat::FetchElTot:I,j="<<j<<",pv="<<pv<<",E="<<lastX[j].first<<",T="
            <<lastX[j].second<<G4endl;
#endif
      if(j!=lastK) pv*=edl;
    }
    i++;                                 // Make a new record to AMDB and position on it
    vI.push_back(lastI);
    vM.push_back(lastM);
    vK.push_back(lastK);
    vX.push_back(lastX);
  }
  else                                   // The A value was found in AMDB
  {
    lastI=vI[i];
    lastM=vM[i];
    lastK=vK[i];
    lastX=vX[i];
    G4int nextK=lastK+1;
    G4double lpM=lastM+lpi;
#ifdef pdebug
    G4cout<<"G4QFreeScatt::FetchElTo:M="<<lpM<<",l="<<lp<<",K="<<lastK<<",n="<<nlp<<G4endl;
#endif
    if(lp>lpM && lastK<nlp)              // LogTab must be updated
    {
      lastK = static_cast<int>((lp-lpi)/dl)+1; // MaxBin to be initialized in LogTab
#ifdef pdebug
      G4cout<<"G4QFreeScat::FetET: K="<<lastK<<",lp="<<lp<<",li="<<lpi<<",dl="<<dl<<G4endl;
#endif
      if(lastK>nlp)
      {
        lastK=nlp;
        lastM=lpa-lpi;
      }
      else lastM = lastK*dl;             // Calculate max initialized ln(p)-lpi for LogTab
      G4double pv=std::exp(lpM);         // momentum of the last calculated beam
      for(G4int j=nextK; j<=lastK; ++j)  // Calculate LogTab values
      {
        pv*=edl;
        lastX[j]=CalcElTot(pv,ind);
#ifdef pdebug
        G4cout<<"G4QFreeScat::FetchElTot: U:j="<<j<<",p="<<pv<<",E="<<lastX[j].first<<",T="
              <<lastX[j].second<<G4endl;
#endif
      }
    } // End of LogTab update
    if(lastK >= nextK)                   // The AMDB was apdated
    {
      vM[i]=lastM;
      vK[i]=lastK;
    }
  }
  // Now one can use tabeles to calculate the value
  G4double dlp=lp-lpi;                       // Shifted log(p) value
  G4int n=static_cast<int>(dlp/dl);          // Low edge number of the bin
  G4double d=dlp-n*dl;                       // Log shift
  G4double e=lastX[n].first;                 // E-Base
  lastR.first=e+d*(lastX[n+1].first-e)/dl;   // E-Result
  if(lastR.first<0.)  lastR.first = 0.;
  G4double t=lastX[n].second;                // T-Base
  lastR.second=t+d*(lastX[n+1].second-t)/dl; // T-Result
  if(lastR.second<0.) lastR.second= 0.;
#ifdef pdebug
  G4cout<<"=O=>G4QFreeScat::FetchElTot:1st="<<lastR.first<<", 2nd="<<lastR.second<<G4endl;
#endif
  if(lastR.first>lastR.second) lastR.first = lastR.second;
  return lastR;
} // End of FetchElTot

// (Mean Elastic and Mean Total over (Z,N)) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QFreeScattering::GetElTotMean(G4double pIU, G4int hPDG,
                                                             G4int Z,      G4int N)
{
  G4double pGeV=pIU/gigaelectronvolt;
  if(hPDG==90001000) hPDG=2212;
  if(hPDG==90000001) hPDG=2112;
  if(hPDG==91000000) hPDG=3122;
#ifdef pdebug
  G4cout<<"->G4QFreeSc::GetElTotMean: P="<<pIU<<",pPDG="<<hPDG<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
  if(Z<1 && N<1)
  {
    G4cout<<"-Warning-G4QFreeScat::GetElTotMean: Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
    return std::make_pair(0.,0.);
  }
  std::pair<G4double,G4double> hp=FetchElTot(pGeV, hPDG, true);
  std::pair<G4double,G4double> hn=FetchElTot(pGeV, hPDG, false);
#ifdef pdebug
  G4cout<<"-OUT->G4QFreeScat::GetElTotMean: hp("<<hp.first<<","<<hp.second<<"), hn("
        <<hn.first<<","<<hn.second<<")"<<G4endl;
#endif
  G4double A=(Z+N)/millibarn;                // To make the result in independent units(IU)
  return std::make_pair((Z*hp.first+N*hn.first)/A,(Z*hp.second+N*hn.second)/A);
} // End of GetElTotMean

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
std::pair<G4LorentzVector,G4LorentzVector> G4QFreeScattering::Scatter(G4int NPDG,
                                     G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  p4M/=megaelectronvolt;   // Convert 4-momenta in MeV (CHIPS units)
  N4M/=megaelectronvolt;
  G4LorentzVector tot4M=N4M+p4M;
  if(pPDG==90001000) pPDG=2212;
  if(pPDG==90000001) pPDG=2112;
  if(pPDG==91000000) pPDG=3122;
#ifdef ppdebug
  G4cout<<"->G4QFR::Scat:p4M="<<p4M<<",N4M="<<N4M<<",t4M="<<tot4M<<",NPDG="<<NPDG<<G4endl;
#endif
  G4double mT=mNeut;
#ifdef ppdebug
  G4int Z=0;
  G4int N=1;
#endif  
  if(NPDG==2212 || NPDG==90001000)
  {
    mT=mProt;
#ifdef ppdebug
    Z=1;
    N=0;
#endif  
  }
  else if(NPDG!=2112 && NPDG!=90000001)
  {
    G4cout<<"Error:G4QFreeScattering::Scatter:NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
    G4Exception("G4QFreeScattering::Scatter:","21",FatalException,"CHIPScomplain");
    //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Use this if not exception
  }
  G4double mT2=mT*mT;    // a squared mass of the free scattered nuclead cluster (FSNC)
  G4double mP2=p4M.m2(); // a projectile squared mass
  G4double E=(tot4M.m2()-mT2-mP2)/(mT+mT); // a projectile energy in the CMS of FSNC
#ifdef pdebug
  G4cout<<"G4QFS::Scat:qM="<<mT<<",qM2="<<mT2<<",pM2="<<mP2<<",totM2="<<tot4M.m2()<<G4endl;
#endif
  G4double E2=E*E;
  if( E < 0. || E2 < mP2)
  {
#ifdef ppdebug
    G4cout<<"-Warning-G4QFS::Scat:*Negative Energy*E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
  G4double pP2=E2-mP2; // Squared Momentum in pseudo laboratory system (final particles)
  G4double P=std::sqrt(pP2); // Momentum in pseudo laboratory system (final particles) ?
#ifdef ppdebug
  G4cout<<"G4QFreeS::Scatter: Before XS, P="<<P<<",Z="<<Z<<",N="<<N<<",PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QFreeScat::Scatter: pPDG="<<pPDG<<G4endl;
  // @@ check a possibility to separate p, n, or alpha (!)
  G4double dmT=mT+mT;
  G4double s_value=dmT*E2+mP2+mT2;         // Mondelstam s (?)
  G4double maxt=dmT*dmT*pP2/s_value;       // max possible |t|
  G4double mint=0.;                        // Prototype of functional rand -t (MeV^2)
  if(P < 14.) mint=maxt*G4UniformRand();   // S-wave
  else                                     // Calculate slopes (no cash !)
  {
    G4double p4=pP2*pP2/1000000.;          // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4double theB1=(7.2+4.32/p4/(p4+12./P))/(1.+2.5/p4); // p4 in GeV, P in MeV
    mint=-std::log(G4UniformRand())/theB1; // t-chan only
  }
#ifdef ppdebug
  G4cout<<"G4QFS::Scat:PDG="<<pPDG<<",P="<<P<<",-t="<<mint<<"<"<<maxt<<", Z="<<Z<<",N="<<N
        <<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);
  else  G4cout<<"*Warning*G4QFreeScattering::Scatter: -t="<<mint<<G4endl;
#endif
  G4double cost=1.-(mint+mint)/maxt;        // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QFS::Scat:-t="<<mint<<"<"<<maxt<<", cost="<<cost<<", Z="<<Z<<",N="<<N<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost > 1.) cost = 1.;
    else if(cost <-1.) cost =-1.;
    else
    {
      G4cerr<<"G4QFreeScatter::Scat:*NAN* cost="<<cost<<",-t="<<mint<<",tm="<<maxt<<G4endl;
      return std::make_pair(G4LorentzVector(0.,0.,0.,0.), p4M); // Do Nothing Action
    }
  }
  G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(p4M, reco4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QFS::Scat:t="<<tot4M<<tot4M.m()<<",mT="<<mT<<",mP="<<std::sqrt(mP2)<<G4endl;
    //G4Exception("G4QFS::Scat:","009",FatalException,"Decay of ElasticComp");
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
#ifdef ppdebug
  G4cout<<"G4QFS::Scat:p4M="<<p4M<<"+r4M="<<reco4M<<",dr="<<dir4M<<",t4M="<<tot4M<<G4endl;
#endif
  return std::make_pair( reco4M*megaelectronvolt, p4M*megaelectronvolt ); // The Result
} // End of Scatter

G4QHadronVector* G4QFreeScattering::InElF(G4int NPDG, G4LorentzVector N4M,
                                          G4int pPDG, G4LorentzVector p4M)
{
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass(); // Pi+ (Pi-: -211)
  static const G4double mK   = G4QPDGCode(321).GetMass(); // K+  (K- : -321)
  static const G4double mK0  = G4QPDGCode(311).GetMass(); // aK0 (K0 : -311)
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  //static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mXiM = G4QPDGCode(3312).GetMass();
  static const G4double mXiZ = G4QPDGCode(3322).GetMass();
  //static const G4double mOmM = G4QPDGCode(3334).GetMass();
  static const G4double third =1./3.;
  static const G4double twothd =2./3.;

  p4M/=megaelectronvolt;   // Convert 4-momenta in MeV (CHIPS units)
  N4M/=megaelectronvolt;
  G4LorentzVector c4M=N4M+p4M;
  if(pPDG==90001000) pPDG=2212;
  if(pPDG==90000001) pPDG=2112;
  if(pPDG==91000000) pPDG=3122;
#ifdef ppdebug
  G4cout<<"->G4QFR::InElF: p4M="<<p4M<<",N4M="<<N4M<<", c4M="<<c4M<<",NPDG="<<NPDG<<G4endl;
#endif
  G4double mT=mNeut;
  G4int Z=0;
  //G4int N=1;
  if(NPDG==2212 || NPDG==90001000)
  {
    NPDG=2212;
    mT=mProt;
    Z=1;
    //N=0;
  }
  else if(NPDG!=2112 && NPDG!=90000001)
  {
    G4cout<<"Error:G4QFreeScattering::InElF:NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
    G4Exception("G4QFreeScattering::InElF:","21",FatalException,"CHIPS complain");
  }
  else NPDG=2112;
  G4double mC2=c4M.m2();
  G4double mC=0.;
  if( mC2 < -0.0000001 )
  {
#ifdef ppdebug
    G4cout<<"-Warning-G4QFS::InElF: Negative compoundMass ="<<mC2<<", c4M="<<c4M<<G4endl;
#endif
    return 0; // Do Nothing Action
  }
  else if ( mC2 > 0.0000001) mC=std::sqrt(mC2); 
#ifdef ppdebug
  G4cout<<"G4QFS::InElF: mC ="<<mC<<G4endl;
#endif
  G4double mP=0.;
  G4double mP2=p4M.m2(); // a projectile squared mass
  if( mP2 < -0.0000001 )
  {
#ifdef ppdebug
    G4cout<<"-Warning-G4QFS::InElF: Negative projectileMass ="<<mC2<<", c4M="<<c4M<<G4endl;
#endif
    return 0; // Do Nothing Action
  }
  else if ( mP2 > 0.0000001) mP=std::sqrt(mP2); 
#ifdef pdebug
  G4cout<<"G4QFS::InElF:mT("<<mT<<")+mP("<<mP<<")+mPi0="<<mP+mT+mPi0<<"<? mC="<<mC<<G4endl;
#endif
  if(pPDG > 3334 || pPDG < -321) // Annihilation/Inelastic of anti-barions not implemented
  {
    G4cout<<"-Warning-G4QFreeScat::InElF: pPDG="<<pPDG<<G4endl;
    return 0; // Do Nothing Action
  }
  if(pPDG==130 || pPDG==310)
  {
    if(G4UniformRand()>.5) pPDG = 311; //  K0
    else                   pPDG =-311; // aK0
  }
  G4int mPDG=111;                     // Additional newtral pion by default
  G4double mM=mPi0;                   // Default Pi0 mass
  G4double r=G4UniformRand();         // A random number to split pi0/pi
  if      (pPDG == 2212) // p-N
  {
    if(Z) // p-p
    {
      if(r<twothd) // -> n + Pi+ + p
      {
        pPDG=2112; // n
        mP=mNeut;
        mPDG= 211; // Pi+
        mM=mPi;
      }
    }
    else  // p-n
    {
      if(r<third) // -> n + Pi+ + n
      {
        pPDG=2112; // n
        mP=mNeut;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<twothd) // -> p + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 2112) // n-N
  {
    if(Z) // n-p
    {
      if(r<third) // -> n + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<twothd) // -> p + Pi- + p
      {
        pPDG=2212; // p
        mP=mProt;
        mPDG=-211; // Pi-
        mM=mPi;
      }
    }
    else  // n-n
    {
      if(r<twothd) // -> p + Pi- + n
      {
        pPDG=2212; // p
        mP=mProt;
        mPDG= -211; // Pi-
        mM=mPi;
      }
    }
  }
  else if (pPDG ==  211) // pip-N
  {
    if(Z) // pip-p
    {
      if(r<0.5) // -> Pi+ + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // pip-n
    {
      if(r<third) // -> Pi0 + Pi0 + p
      {
        pPDG= 111; // Pi0
        mP=mPi0;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<twothd) // -> Pi+ + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == -211) // pim-N
  {
    if(Z) // pim-p
    {
      if(r<third) // -> Pi0 + Pi0 + n
      {
        pPDG= 111; // Pi0
        mP=mPi0;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<twothd) // -> Pi- + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // pim-n
    {
      if(r<0.5) // -> Pi- + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == -321) // Km-N
  {
    if(Z) // Km-p
    {
      if(r<0.25) // K0 + Pi0 + n
      {
        pPDG=-311;  // K0
        mP=mK0;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.5) // -> K- + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.75) // -> K0 + Pi- + p
      {
        pPDG=-311;  // K0
        mP=mK0;
        mPDG=-211; // Pi-
        mM=mPi;
      }
    }
    else  // Km-n
    {
      if(r<third) // -> K- + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<twothd) // -> K0 + Pi- + n
      {
        pPDG=-311;  // K0
        mP=mK0;
        mPDG=-211; // Pi-
        mM=mPi;
      }
    }
  }
  else if (pPDG == -311) // K0-N
  {
    if(Z) // K0-p
    {
      if(r<third) // K- + Pi+ + p
      {
        pPDG=-321;  // K-
        mP=mK;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<twothd) // -> K0 + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // K0-n
    {
      if(r<0.25) // -> K- + Pi+ + n
      {
        pPDG=-321; // K-
        mP=mK;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.5) // -> K- + Pi0 + p
      {
        pPDG=-321; // K-
        mP=mK;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<0.75) // -> K0 + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG ==  321) // Kp-N
  {
    if(Z) // Kp-p
    {
      if(r<third) // -> K+ + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<twothd) // -> aK0 + Pi+ + p
      {
        pPDG= 311; // aK0
        mP=mK0;
        mPDG= 211; // Pi+
        mM=mPi;
      }
    }
    else  // Kp-n
    {
      if(r<0.25) // -> aK0 + Pi+ + n
      {
        pPDG= 311; // aK0
        mP=mK0;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.5) // -> Kp + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<0.75) // -> aK0 + Pi0 + p
      {
        pPDG= 311; // aK0
        mP=mK0;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG ==  311) // aK0-N
  {
    if(Z) // aK0-p
    {
      if(r<0.25) // -> K+ + Pi- + p
      {
        pPDG= 321; // K+
        mP=mK;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<0.5) // -> aK0 + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.75) // -> K+ + Pi0 + n
      {
        pPDG= 321; // K+
        mP=mK;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // aK0-n
    {
      if(r<third) // -> aK0 + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<twothd) // -> K+ + Pi- + n
      {
        pPDG= 321; // K+
        mP=mK;
        mPDG=-211; // Pi-
        mM=mPi;
      }
    }
  }
  else if (pPDG == 3122 || pPDG== 3212) // Lambda/Sigma0-N
  {
    if(pPDG == 3212)
    {
      pPDG=3122;
      mP=mLamb;
    }
    if(Z) // L/S0-p
    {
      if(r<0.2) // -> SigP + Pi0 + n
      {
        pPDG=3222; // SigP
        mP=mSigP;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.4) // -> SigP + Pi- + p
      {
        pPDG=3222; // SigP
        mP=mSigP;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<0.6) // -> SigM + Pi+ + p
      {
        pPDG=3112; // SigM
        mP=mSigM;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.8) // -> Lamb + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // L/S0-n
    {
      if(r<0.2) // -> SigM + Pi0 + p
      {
        pPDG=3112; // SigM
        mP=mSigM;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<0.4) // -> SigP + Pi- + n
      {
        pPDG=3222; // SigP
        mP=mSigP;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<0.6) // -> SigM + Pi+ + n
      {
        pPDG=3112; // SigM
        mP=mSigM;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.8) // -> Lamb + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 3112) // Sigma- -N
  {
    if(Z) // SigM-p
    {
      if(r<0.25) // -> Lamb + Pi0 + n
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.5) // -> Lamb + Pi- + p
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<0.75) // -> SigM + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // SigM-n
    {
      if(r<third) // -> Lamb + Pi- + n
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<twothd) // -> SigM + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 3222) // Sigma+ -N
  {
    if(Z) // SigP-p
    {
      if(r<third) // -> Lamb + Pi+ + p
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<twothd) // -> SigP + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // SigP-n
    {
      if(r<0.25) // -> Lamb + Pi0 + p
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<0.5) // -> Lamb + Pi+ + n
      {
        pPDG=3122; // Lamb
        mP=mLamb;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.75) // -> SigP + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 3312) // Xi- -N
  {
    if(Z) // XiM-p
    {
      if(r<0.25) // -> Xi0 + Pi0 + n
      {
        pPDG=3322; // Xi0
        mP=mXiZ;
        NPDG=2112; // n
        mT=mNeut;
      }
      else if(r<0.5) // -> Xi0 + Pi- + p
      {
        pPDG=3322; // Xi0
        mP=mXiZ;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<0.75) // -> XiM + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // XiM-n
    {
      if(r<third) // -> Xi0 + Pi- + n
      {
        pPDG=3322; // Xi0
        mP=mXiZ;
        mPDG=-211; // Pi-
        mM=mPi;
      }
      else if(r<twothd) // -> XiM + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 3322) // Xi0 -N
  {
    if(Z) // Xi0-p
    {
      if(r<third) // -> Xi- + Pi+ + p
      {
        pPDG=3312; // Xi-
        mP=mXiM;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<twothd) // -> Xi0 + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      }
    }
    else  // Xi0-n
    {
      if(r<0.25) // -> Xi- + Pi0 + p
      {
        pPDG=3312; // Xi-
        mP=mXiM;
        NPDG=2212; // p
        mT=mProt;
      }
      else if(r<0.5) // -> Xi- + Pi+ + n
      {
        pPDG=3312; // Xi-
        mP=mXiM;
        mPDG= 211; // Pi+
        mM=mPi;
      }
      else if(r<0.75) // -> Xi0 + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else if (pPDG == 3334) // Om- -N
  {
    if(Z) // OmM-p
    {
      if(r<0.5)    // -> Om- + Pi+ + n
      {
        mPDG= 211; // Pi+
        mM=mPi;
        NPDG=2112; // n
        mT=mNeut;
      } 
    }
    else  // OmM-n
    {
      if(r<0.5)    // -> Om- + Pi- + p
      {
        mPDG=-211; // Pi-
        mM=mPi;
        NPDG=2212; // p
        mT=mProt;
      }
    }
  }
  else
  {
    G4cout<<"*Error*G4QFreeScattering::InElF: PDG="<<pPDG
          <<", while it is defined only for p,n,hyperons(not Omega),pi,K/antiK"<<G4endl;
    G4Exception("G4QFreeScattering::InElF:","22",FatalException,"CHIPS_crash");
  }
  if     (mC-mP-mM-mT <-0.000001) return 0;
  G4QHadronVector* TripQH = new G4QHadronVector; // Proto of the Result
  G4LorentzVector m4M(0.,0.,0.,mM);
  if(mC-mP-mM-mT < 0.000001) // Equal share
  {
    p4M=(mP/mC)*c4M;
    m4M=(mM/mC)*c4M;
    N4M=(mT/mC)*c4M;
  }
  else
  {
    p4M=G4LorentzVector(0.,0.,0.,mP);
    N4M=G4LorentzVector(0.,0.,0.,mT);
    if(!G4QHadron(c4M).DecayIn3(p4M,m4M,N4M))
    {
      G4ExceptionDescription ed;
      ed << "DecayIn3, TotM=" << mC /* << " <? " << mT + mP + mM */ << G4endl;
      G4Exception("G4QFreeScattering::InElF()","HAD_CHPS_0027", FatalException, ed);
    }
    G4QHadron* h1 = new G4QHadron(pPDG, p4M);
    TripQH->push_back(h1); // (delete equivalent, responsibility of users)
#ifdef debug
    G4cout << "G4QFreeScat::InElF: H1=" << pPDG << p4M << G4endl;
#endif
    G4QHadron* h2 = new G4QHadron(mPDG, m4M);
    TripQH->push_back(h2); // (delete equivalent, responsibility of users)
#ifdef debug
    G4cout << "G4QFreeScat::InElF: H2=" << mPDG << m4M << G4endl;
#endif
    G4QHadron* h3 = new G4QHadron(NPDG, N4M);
    TripQH->push_back(h3); // (delete equivalent, responsibility of users)
#ifdef debug
    G4cout << "G4QFreeScat::InElF: H3=" << NPDG << N4M << G4endl;
#endif
  }
  return TripQH; // The Result
} // End of Scatter
