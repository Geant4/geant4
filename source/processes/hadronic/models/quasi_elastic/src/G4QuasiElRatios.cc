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
// $Id: G4QuasiElRatios.cc 90573 2015-06-04 09:38:49Z gcosmo $
//
//
// G4 Physics class: G4QuasiElRatios for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
// 
// ----------------------------------------------------------------------
// This class has been extracted from the CHIPS model. 
// All the dependencies on CHIPS classes have been removed.
// Short description: Provides percentage of quasi-free and quasi-elastic
// reactions in the inelastic reactions.
// ----------------------------------------------------------------------


#include "G4QuasiElRatios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"

namespace  {
    const G4int    nps=150;        // Number of steps in the R(s) LinTable
    const G4int    mps=nps+1;      // Number of elements in the R(s) LinTable
    const G4double sma=150.;       // The first LinTabEl(s=0)=1., s>sma -> logTab
    const G4double ds=sma/nps;     // Step of the linear Table
    const G4int    nls=100;        // Number of steps in the R(lns) logTable
    const G4int    mls=nls+1;      // Number of elements in the R(lns) logTable
    const G4double lsi=5.;         // The min ln(s) logTabEl(s=148.4 < sma=150.)
    const G4double lsa=9.;         // The max ln(s) logTabEl(s=148.4 - 8103. mb)
    const G4double mi=std::exp(lsi);// The min s of logTabEl(~ 148.4 mb)
    const G4double min_s=std::exp(lsa);// The max s of logTabEl(~ 8103. mb)
    const G4double dls=(lsa-lsi)/nls;// Step of the logarithmic Table
    const G4double edls=std::exp(dls);// Multiplication step of the logarithmic Table
    const G4double toler=.01;      // The tolarence mb defining the same cross-section
    const G4double C=1.246;
    const G4double lmi=3.5;       // min of (lnP-lmi)^2 parabola
    const G4double pbe=.0557;     // elastic (lnP-lmi)^2 parabola coefficient
    const G4double pbt=.3;        // total (lnP-lmi)^2 parabola coefficient
    const G4double pmi=.1;        // Below that fast LE calculation is made
    const G4double pma=1000.;     // Above that fast HE calculation is made
    const G4int    nlp=300;         // Number of steps in the S(lnp) logTable(5% step)
    const G4int    mlp=nlp+1;       // Number of elements in the S(lnp) logTable
    const G4double lpi=-5.;         // The min ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
    const G4double lpa=10.;         // The max ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
    const G4double mip=std::exp(lpi);// The min p of logTabEl(~ 6.7 MeV/c)
    const G4double map=std::exp(lpa);// The max p of logTabEl(~ 22. TeV)
    const G4double dlp=(lpa-lpi)/nlp;// Step of the logarithmic Table
    const G4double edlp=std::exp(dlp);// Multiplication step of the logarithmic Table
}


G4QuasiElRatios::G4QuasiElRatios()
{
    vT = new std::vector<G4double*>;
    vL = new std::vector<G4double*>;
    vX = new std::vector<std::pair<G4double,G4double>*>;
    
    lastSRatio=0.;             // The last sigma value for which R was calculated
    lastRRatio=0.;             // The last ratio R which was calculated
    lastARatio=0;             // theLast of calculated A
    lastHRatio=0.;            // theLast of max s initialized in the LinTable
    lastNRatio=0;             // theLast of topBin number initialized in LinTable
    lastMRatio=0.;            // theLast of rel max ln(s) initialized in LogTable
    lastKRatio=0;             // theLast of topBin number initialized in LogTable
    lastTRatio=0;             // theLast of pointer to LinTable in the C++ heap
    lastLRatio=0;             // theLast of pointer to LogTable in the C++ heap
    lastPtot=0.;              // The last momentum for which XS was calculated
    lastHtot=0;               // The last projPDG for which XS was calculated
    lastFtot=true;            // The last nucleon for which XS was calculated
    lastItot=0;              // The Last index for which XS was calculated
    lastMtot=0.;             // The Last rel max ln(p) initialized in LogTable
    lastKtot=0;             // The Last topBin number initialized in LogTable
    lastXtot = new std::pair<G4double,G4double>; // The Last ETPointers to LogTable in heap
    
    PCSmanager=(G4ChipsProtonElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name());
    
    NCSmanager=(G4ChipsNeutronElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name());
}

G4QuasiElRatios::~G4QuasiElRatios()
{
    std::vector<G4double*>::iterator pos;
    for(pos=vT->begin(); pos<vT->end(); pos++)
    { delete [] *pos; }
    vT->clear();
    for(pos=vL->begin(); pos<vL->end(); pos++)
    { delete [] *pos; }
    vL->clear();
    
    std::vector<std::pair<G4double,G4double>*>::iterator pos2;
    for(pos2=vX->begin(); pos2<vX->end(); pos2++)
    { delete [] *pos2; }
    vX->clear();
}

// Calculation of pair(QuasiFree/Inelastic,QuasiElastic/QuasiFree)
std::pair<G4double,G4double> G4QuasiElRatios::GetRatios(G4double pIU, G4int pPDG,
                                                         G4int tgZ,    G4int tgN)
{
    G4double R=0.;
    G4double QF2In=1.;                        // Prototype of QuasiFree/Inel ratio for hN_tot
    G4int tgA=tgZ+tgN;
    if(tgA<2) return std::make_pair(QF2In,R); // No quasi-elastic on the only nucleon
    std::pair<G4double,G4double> ElTot=GetElTot(pIU, pPDG, tgZ, tgN); // mean hN El&Tot(IU)
    //if( ( (pPDG>999 && pIU<227.) || pIU<27.) && tgA>1) R=1.; // @@ TMP to accelerate @lowE
    if(pPDG>999 && pIU<227. && tgZ+tgN>1) R=1.;                // To accelerate @lowE
    else if(ElTot.second>0.)
    {
        R=ElTot.first/ElTot.second;             // El/Total ratio (does not depend on units
        QF2In=GetQF2IN_Ratio(ElTot.second/millibarn, tgZ+tgN);   // QuasiFree/Inelastic ratio
    }
    return std::make_pair(QF2In,R);
}

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section (mb) and A
G4double G4QuasiElRatios::GetQF2IN_Ratio(G4double m_s, G4int A)
{
    // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
    if(m_s<toler || A<2) return 1.;
    if(m_s>min_s) return 0.;
    
    if(A>238)
    {
        G4cout<<"-Warning-G4QuasiElRatio::GetQF2IN_Ratio:A="<<A<<">238, return zero"<<G4endl;
        return 0.;
    }
    G4int nDB=vARatio.size();                  // A number of nuclei already initialized in AMDB
    if(nDB && lastARatio==A && m_s==lastSRatio) return lastRRatio;  // VI do not use tolerance
    G4bool found=false;
    G4int i=-1;
    if(nDB) for (i=0; i<nDB; i++) if(A==vARatio[i]) // Search for this A in AMDB
    {
        found=true;                         // The A value is found
        break;
    }
    if(!nDB || !found)                    // Create new line in the AMDB
    {
        lastARatio = A;
        lastTRatio = new G4double[mps];          // Create the linear Table
        lastNRatio = static_cast<int>(m_s/ds)+1;   // MaxBin to be initialized
        if(lastNRatio>nps)
        {
            lastNRatio=nps;
            lastHRatio=sma;
        }
        else lastHRatio = lastNRatio*ds;              // Calculate max initialized s for LinTab
        G4double sv=0;
        lastTRatio[0]=1.;
        for(G4int j=1; j<=lastNRatio; j++)       // Calculate LogTab values
        {
            sv+=ds;
            lastTRatio[j]=CalcQF2IN_Ratio(sv,A);
        }
        lastLRatio=new G4double[mls];            // Create the logarithmic Table
        if(m_s>sma)                           // Initialize the logarithmic Table
        {
	  G4double ls=std::log(m_s);
            lastKRatio = static_cast<int>((ls-lsi)/dls)+1; // MaxBin to be initialized in LogTaB
            if(lastKRatio>nls)
            {
                lastKRatio=nls;
                lastMRatio=lsa-lsi;
            }
            else lastMRatio = lastKRatio*dls;            // Calculate max initialized ln(s)-lsi for LogTab
            sv=mi;
            for(G4int j=0; j<=lastKRatio; j++)     // Calculate LogTab values
            {
                lastLRatio[j]=CalcQF2IN_Ratio(sv,A);
                if(j!=lastKRatio) sv*=edls;
            }
        }
        else                                // LogTab is not initialized
        {
            lastKRatio = 0;
            lastMRatio = 0.;
        }
        i++;                                // Make a new record to AMDB and position on it
        vARatio.push_back(lastARatio);
        vHRatio.push_back(lastHRatio);
        vNRatio.push_back(lastNRatio);
        vMRatio.push_back(lastMRatio);
        vKRatio.push_back(lastKRatio);
        vT->push_back(lastTRatio);
        vL->push_back(lastLRatio);
    }
    else                                  // The A value was found in AMDB
    {
        lastARatio=vARatio[i];
        lastHRatio=vHRatio[i];
        lastNRatio=vNRatio[i];
        lastMRatio=vMRatio[i];
        lastKRatio=vKRatio[i];
        lastTRatio=(*vT)[i];
        lastLRatio=(*vL)[i];
        if(m_s>lastHRatio)                          // At least LinTab must be updated
        {
            G4int nextN=lastNRatio+1;               // The next bin to be initialized
            if(lastNRatio<nps)
            {
	      G4double sv=lastHRatio; // bug fix by WP

                lastNRatio = static_cast<int>(m_s/ds)+1;// MaxBin to be initialized
                if(lastNRatio>nps)
                {
                    lastNRatio=nps;
                    lastHRatio=sma;
                }
                else lastHRatio = lastNRatio*ds;           // Calculate max initialized s for LinTab

                for(G4int j=nextN; j<=lastNRatio; j++)// Calculate LogTab values
                {
                    sv+=ds;
                    lastTRatio[j]=CalcQF2IN_Ratio(sv,A);
                }
            } // End of LinTab update
            if(lastNRatio>=nextN)
            {
                vHRatio[i]=lastHRatio;
                vNRatio[i]=lastNRatio;
            }
            G4int nextK=lastKRatio+1;
            if(!lastKRatio) nextK=0;
            if(m_s>sma && lastKRatio<nls)             // LogTab must be updated
            {
                G4double sv=std::exp(lastMRatio+lsi); // Define starting poit (lastM will be changed)
                G4double ls=std::log(m_s);
                lastKRatio = static_cast<int>((ls-lsi)/dls)+1; // MaxBin to be initialized in LogTaB
                if(lastKRatio>nls)
                {
                    lastKRatio=nls;
                    lastMRatio=lsa-lsi;
                }
                else lastMRatio = lastKRatio*dls;           // Calculate max initialized ln(s)-lsi for LogTab
                for(G4int j=nextK; j<=lastKRatio; j++)// Calculate LogTab values
                {
                    sv*=edls;
                    lastLRatio[j]=CalcQF2IN_Ratio(sv,A);
                }
            } // End of LogTab update
            if(lastKRatio>=nextK)
            {
                vMRatio[i]=lastMRatio;
                vKRatio[i]=lastKRatio;
            }
        }
    }
    // Now one can use tabeles to calculate the value
    if(m_s<sma)                             // Use linear table
    {
        G4int n=static_cast<int>(m_s/ds);     // Low edge number of the bin
        G4double d=m_s-n*ds;                  // Linear shift
        G4double v=lastTRatio[n];                // Base
        lastRRatio=v+d*(lastTRatio[n+1]-v)/ds;        // Result
    }
    else                                  // Use log table
    {
        G4double ls=std::log(m_s)-lsi;        // ln(s)-l_min
        G4int n=static_cast<int>(ls/dls);    // Low edge number of the bin
        G4double d=ls-n*dls;                 // Log shift
        G4double v=lastLRatio[n];                // Base
        lastRRatio=v+d*(lastLRatio[n+1]-v)/dls;        // Result
    }
    if(lastRRatio<0.) lastRRatio=0.;
    if(lastRRatio>1.) lastRRatio=1.;
    return lastRRatio;
} // End of CalcQF2IN_Ratio

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section and A
G4double G4QuasiElRatios::CalcQF2IN_Ratio(G4double m_s, G4int A)
{
    G4double s2=m_s*m_s;
    G4double s4=s2*s2;
    G4double ss=std::sqrt(std::sqrt(m_s));
    G4double P=7.48e-5*s2/(1.+8.77e12/s4/s4/s2);
    G4double E=.2644+.016/(1.+std::exp((29.54-m_s)/2.49));
    G4double F=ss*.1526*std::exp(-s2*ss*.0000859);
    return C*std::exp(-E*std::pow(G4double(A-1.),F))/std::pow(G4double(A),P);
} // End of CalcQF2IN_Ratio

// Calculatio pair(hN_el,hN_tot) (mb): p in GeV/c, index(PDG,F) (see FetchElTot)
std::pair<G4double,G4double> G4QuasiElRatios::CalcElTot(G4double p, G4int I)
{
    // ---------> Each parameter set can have not more than nPoints=128 parameters

    G4double El=0.;                      // prototype of the elastic hN cross-section
    G4double To=0.;                      // prototype of the total hN cross-section
    if(p<=0.)
    {
        G4cout<<"-Warning-G4QuasiElRatios::CalcElTot: p="<<p<<" is zero or negative"<<G4endl;
        return std::make_pair(El,To);
    }
    if     (!I)                          // pp/nn
    {
        if(p<pmi)
        {
            G4double p2=p*p;
            El=1./(.00012+p2*.2);
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
            G4double LE=1./(.00012+p2*.2);
            G4double lp=std::log(p)-lmi;
            G4double lp2=lp*lp;
            G4double rp2=1./p2;
            El=LE+(pbe*lp2+6.72+32.6/p)/(1.+rp2/p);
            To=LE+(pbt*lp2+38.2+52.7*rp2)/(1.+2.72*rp2*rp2);
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
        G4cout<<"*Error*G4QuasiElRatios::CalcElTot:ind="<<I<<" is not defined (0-7)"<<G4endl;
        G4Exception("G4QuasiElRatios::CalcElTot:","23",FatalException,"QEcrash");
    }
    if(El>To) El=To;
    return std::make_pair(El,To);
} // End of CalcElTot

// For hadron PDG with momentum Mom (GeV/c) on N (p/n) calculate <sig_el,sig_tot> pair (mb)
std::pair<G4double,G4double> G4QuasiElRatios::GetElTotXS(G4double p, G4int PDG, G4bool F)
{
    G4int ind=0;                                 // Prototype of the reaction index
    G4bool kfl=true;                             // Flag of K0/aK0 oscillation
    G4bool kf=false;
    if(PDG==130||PDG==310)
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
        G4cout<<"*Error*G4QuasiElRatios::CalcElTotXS: PDG="<<PDG
        <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
        G4Exception("G4QuasiElRatio::CalcElTotXS:","22",FatalException,"QEcrash");
    }
    return CalcElTot(p,ind);
}

// Calculatio pair(hN_el,hN_tot)(mb): p in GeV/c, F=true -> N=proton, F=false -> N=neutron
std::pair<G4double,G4double> G4QuasiElRatios::FetchElTot(G4double p, G4int PDG, G4bool F)
{
    // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
    G4int nDB=vItot.size();                   // A number of hadrons already initialized in AMDB
    if(nDB && lastHtot==PDG && lastFtot==F && p>0. && p==lastPtot) return lastRtot;// VI don't use toler.
    //  if(nDB && lastH==PDG && lastF==F && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
    lastHtot=PDG;
    lastFtot=F;
    G4int ind=-1;                          // Prototipe of the index of the PDG/F combination
    // i=0: pp(nn), i=1: np(pn), i=2: pimp(pipn), i=3: pipp(pimn), i=4: Kmp(Kmn,K0n,K0p),
    // i=5: Kpp(Kpn,aK0n,aK0p), i=6: Hp(Hn), i=7: app(apn,ann,anp) 
    G4bool kfl=true;                             // Flag of K0/aK0 oscillation
    G4bool kf=false;
    if(PDG==130||PDG==310)
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
        G4cout<<"*Error*G4QuasiElRatios::FetchElTot: PDG="<<PDG
        <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
        G4Exception("G4QuasiELRatio::FetchElTot:","22",FatalException,"QECrash");
    }
    if(nDB && lastItot==ind && p>0. && p==lastPtot) return lastRtot;  // VI do not use toler
    //  if(nDB && lastI==ind && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
    if(p<=mip || p>=map) return CalcElTot(p,ind);   // @@ Slow calculation ! (Warning?)
    G4bool found=false;
    G4int i=-1;
    if(nDB) for (i=0; i<nDB; i++) if(ind==vItot[i])  // Sirch for this index in AMDB
    {
        found=true;                                 // The index is found
        break;
    }
    G4double lp=std::log(p);
    if(!nDB || !found)                            // Create new line in the AMDB
    {
        lastXtot = new std::pair<G4double,G4double>[mlp]; // Create logarithmic Table for ElTot
        lastItot = ind;                                // Remember the initialized inex
        lastKtot = static_cast<int>((lp-lpi)/dlp)+1;    // MaxBin to be initialized in LogTaB
        if(lastKtot>nlp)
        {
            lastKtot=nlp;
            lastMtot=lpa-lpi;
        }
        else lastMtot = lastKtot*dlp;               // Calculate max initialized ln(p)-lpi for LogTab
        G4double pv=mip;
        for(G4int j=0; j<=lastKtot; j++)        // Calculate LogTab values
        {
            lastXtot[j]=CalcElTot(pv,ind);
            if(j!=lastKtot) pv*=edlp;
        }
        i++;                                 // Make a new record to AMDB and position on it
        vItot.push_back(lastItot);
        vMtot.push_back(lastMtot);
        vKtot.push_back(lastKtot);
        vX->push_back(lastXtot);
    }
    else                                   // The A value was found in AMDB
    {
        lastItot=vItot[i];
        lastMtot=vMtot[i];
        lastKtot=vKtot[i];
        lastXtot=(*vX)[i];
        G4int nextK=lastKtot+1;
        G4double lpM=lastMtot+lpi;
        if(lp>lpM && lastKtot<nlp)              // LogTab must be updated
        {
            lastKtot = static_cast<int>((lp-lpi)/dlp)+1; // MaxBin to be initialized in LogTab
            if(lastKtot>nlp)
            {
                lastKtot=nlp;
                lastMtot=lpa-lpi;
            }
            else lastMtot = lastKtot*dlp;           // Calculate max initialized ln(p)-lpi for LogTab
            G4double pv=std::exp(lpM);       // momentum of the last calculated beam
            for(G4int j=nextK; j<=lastKtot; j++)// Calculate LogTab values
            {
                pv*=edlp;
                lastXtot[j]=CalcElTot(pv,ind);
            }
        } // End of LogTab update
        if(lastKtot>=nextK)                   // The AMDB was apdated
        {
            vMtot[i]=lastMtot;
            vKtot[i]=lastKtot;
        }
    }
    // Now one can use tabeles to calculate the value
    G4double dlpp=lp-lpi;                       // Shifted log(p) value
    G4int n=static_cast<int>(dlpp/dlp);          // Low edge number of the bin
    G4double d=dlpp-n*dlp;                       // Log shift
    G4double e=lastXtot[n].first;                 // E-Base
    lastRtot.first=e+d*(lastXtot[n+1].first-e)/dlp;   // E-Result
    if(lastRtot.first<0.)  lastRtot.first = 0.;
    G4double t=lastXtot[n].second;                // T-Base
    lastRtot.second=t+d*(lastXtot[n+1].second-t)/dlp; // T-Result
    if(lastRtot.second<0.) lastRtot.second= 0.;
    if(lastRtot.first>lastRtot.second) lastRtot.first = lastRtot.second;
    return lastRtot;
} // End of FetchElTot

// (Mean Elastic and Mean Total) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiElRatios::GetElTot(G4double pIU, G4int hPDG,
                                                        G4int Z,       G4int N)
{
    G4double pGeV=pIU/gigaelectronvolt;
    if(Z<1 && N<1)
    {
        G4cout<<"-Warning-G4QuasiElRatio::GetElTot:Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
        return std::make_pair(0.,0.);
    }
    std::pair<G4double,G4double> hp=FetchElTot(pGeV, hPDG, true);
    std::pair<G4double,G4double> hn=FetchElTot(pGeV, hPDG, false);
    G4double A=(Z+N)/millibarn;                // To make the result in independent units(IU)
    return std::make_pair((Z*hp.first+N*hn.first)/A,(Z*hp.second+N*hn.second)/A);
} // End of GetElTot

// (Mean Elastic and Mean Total) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiElRatios::GetChExFactor(G4double pIU, G4int hPDG,
                                                             G4int Z, G4int N)
{
    G4double pGeV=pIU/gigaelectronvolt;
    G4double resP=0.;
    G4double resN=0.;
    if(Z<1 && N<1)
    {
        G4cout<<"-Warning-G4QuasiElRatio::GetChExF:Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
        return std::make_pair(resP,resN);
    }
    G4double A=Z+N;
    G4double pf=0.;                              // Possibility to interact with a proton
    G4double nf=0.;                              // Possibility to interact with a neutron
    if   (hPDG==-211||hPDG==-321||hPDG==3112||hPDG==3212||hPDG==3312) pf=Z/(A+N);
    else if(hPDG==211||hPDG==321||hPDG==3222||hPDG==3212||hPDG==3322) nf=N/(A+Z);
    else if(hPDG==-311||hPDG==311||hPDG==130||hPDG==310)
    {
        G4double dA=A+A;
        pf=Z/(dA+N+N);
        nf=N/(dA+Z+Z);
    }
    G4double mult=1.;  // Factor of increasing multiplicity ( ? @@)
    if(pGeV>.5)
    {
        mult=1./(1.+std::log(pGeV+pGeV))/pGeV;
        if(mult>1.) mult=1.;
    }
    if(pf)
    {
        std::pair<G4double,G4double> hp=FetchElTot(pGeV, hPDG, true);
        resP=pf*(hp.second/hp.first-1.)*mult;
    }
    if(nf)
    {
        std::pair<G4double,G4double> hn=FetchElTot(pGeV, hPDG, false);
        resN=nf*(hn.second/hn.first-1.)*mult;
    }
    return std::make_pair(resP,resN);
} // End of GetChExFactor

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
std::pair<G4LorentzVector,G4LorentzVector> G4QuasiElRatios::Scatter(G4int NPDG,
                                                                     G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
    static const G4double mNeut= G4Neutron::Neutron()->GetPDGMass();
    static const G4double mProt= G4Proton::Proton()->GetPDGMass();
    static const G4double mDeut= G4Deuteron::Deuteron()->GetPDGMass();
    static const G4double mTrit= G4Triton::Triton()->GetPDGMass();
    static const G4double mHel3= G4He3::He3()->GetPDGMass();
    static const G4double mAlph= G4Alpha::Alpha()->GetPDGMass();
    
    G4LorentzVector pr4M=p4M/megaelectronvolt;   // Convert 4-momenta in MeV (keep p4M)
    N4M/=megaelectronvolt;
    G4LorentzVector tot4M=N4M+p4M;
    G4double mT=mNeut;
    G4int Z=0;
    G4int N=1;
    if(NPDG==2212||NPDG==90001000)
    {
        mT=mProt;
        Z=1;
        N=0;
    }
    else if(NPDG==90001001)
    {
        mT=mDeut;
        Z=1;
        N=1;
    }
    else if(NPDG==90002001)
    {
        mT=mHel3;
        Z=2;
        N=1;
    }
    else if(NPDG==90001002)
    {
        mT=mTrit;
        Z=1;
        N=2;
    }
    else if(NPDG==90002002)
    {
        mT=mAlph;
        Z=2;
        N=2;
    }
    else if(NPDG!=2112&&NPDG!=90000001)
    {
        G4cout<<"Error:G4QuasiElRatios::Scatter:NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
        G4Exception("G4QuasiElRatios::Scatter:","21",FatalException,"QEcomplain");
        //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
    }
    G4double mT2=mT*mT;
    G4double mP2=pr4M.m2();
    G4double E=(tot4M.m2()-mT2-mP2)/(mT+mT);
    G4double E2=E*E;
    if(E<0. || E2<mP2)
    {
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
    G4double P=std::sqrt(E2-mP2);                   // Momentum in pseudo laboratory system
    // @@ Temporary NN t-dependence for all hadrons
    if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QE::Scatter: pPDG="<<pPDG<<G4endl;
    G4int PDG=2212;                                                // *TMP* instead of pPDG
    if(pPDG==2112||pPDG==-211||pPDG==-321) PDG=2112;               // *TMP* instead of pPDG
    if(!Z && N==1)                 // Change for Quasi-Elastic on neutron
    {
        Z=1;
        N=0;
        if     (PDG==2212) PDG=2112;
        else if(PDG==2112) PDG=2212;
    }
    G4double xSec=0.;                        // Prototype of Recalculated Cross Section *TMP*
    if(PDG==2212) xSec=PCSmanager->GetChipsCrossSection(P, Z, N, PDG); // P CrossSect *TMP*
    else          xSec=NCSmanager->GetChipsCrossSection(P, Z, N, PDG); // N CrossSect *TMP*
    // @@ check a possibility to separate p, n, or alpha (!)
    if(xSec <= 0.)                                    // The cross-section iz 0 -> Do Nothing
    {
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); //Do Nothing Action
    }
    G4double mint=0.;                        // Prototype of functional rand -t (MeV^2) *TMP*
    if(PDG==2212) mint=PCSmanager->GetExchangeT(Z,N,PDG);// P functional rand -t(MeV^2) *TMP*
    else          mint=NCSmanager->GetExchangeT(Z,N,PDG);// N functional rand -t(MeV^2) *TMP*
    G4double maxt=0.;                                    // Prototype of max possible -t
    if(PDG==2212) maxt=PCSmanager->GetHMaxT();           // max possible -t
    else          maxt=NCSmanager->GetHMaxT();           // max possible -t
    G4double cost=1.-(mint+mint)/maxt; // cos(theta) in CMS
    if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
    {
        if     (cost>1.)  cost=1.;
        else if(cost<-1.) cost=-1.;
        else
        {
            G4double tm=0.;
            if(PDG==2212) tm=PCSmanager->GetHMaxT();
            else          tm=NCSmanager->GetHMaxT();
            G4cerr<<"G4QuasiFreeRatio::Scat:*NAN* cost="<<cost<<",-t="<<mint<<",tm="<<tm<<G4endl;
            return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
        }
    }
    G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
    G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
    if(!RelDecayIn2(tot4M, pr4M, reco4M, dir4M, cost, cost))
    {
        G4cerr<<"G4QFR::Scat:t="<<tot4M<<tot4M.m()<<",mT="<<mT<<",mP="<<std::sqrt(mP2)<<G4endl;
        //G4Exception("G4QFR::Scat:","009",FatalException,"Decay of ElasticComp");
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
    return std::make_pair(reco4M*megaelectronvolt,pr4M*megaelectronvolt); // Result
} // End of Scatter

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
// User should himself change the charge (PDG) (e.g. pn->np, pi+n->pi0p, pi-p->pi0n etc.)
std::pair<G4LorentzVector,G4LorentzVector> G4QuasiElRatios::ChExer(G4int NPDG,
                                                                    G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
    static const G4double mNeut= G4Neutron::Neutron()->GetPDGMass();
    static const G4double mProt= G4Proton::Proton()->GetPDGMass();
    G4LorentzVector pr4M=p4M/megaelectronvolt;          // Convert 4-momenta in MeV(keep p4M)
    N4M/=megaelectronvolt;
    G4LorentzVector tot4M=N4M+p4M;
    G4int Z=0;
    G4int N=1;
    G4int sPDG=0;                                        // PDG code of the scattered hadron
    G4double mS=0.;                                      // proto of mass of scattered hadron
    G4double mT=mProt;                                   // mass of the recoil nucleon
    if(NPDG==2212)
    {
        mT=mNeut;
        Z=1;
        N=0;
        if(pPDG==-211) sPDG=111;                           // pi+    -> pi0
        else if(pPDG==-321)
        {
            sPDG=310;                                        // K+     -> K0S
            if(G4UniformRand()>.5) sPDG=130;                 // K+     -> K0L
        }
        else if(pPDG==-311||pPDG==311||pPDG==130||pPDG==310) sPDG=321;  // K0     -> K+ (?)
        else if(pPDG==3112) sPDG=3212;                     // Sigma- -> Sigma0
        else if(pPDG==3212) sPDG=3222;                     // Sigma0 -> Sigma+
        else if(pPDG==3312) sPDG=3322;                     // Xi-    -> Xi0
    }
    else if(NPDG==2112) // Default
    {
        if(pPDG==211)  sPDG=111;                           // pi+    -> pi0
        else if(pPDG==321)
        {
            sPDG=310;                                        // K+     -> K0S
            if(G4UniformRand()>.5) sPDG=130;                 // K+     -> K0L
        }
        else if(pPDG==-311||pPDG==311||pPDG==130||pPDG==310) sPDG=-321; // K0     -> K- (?)
        else if(pPDG==3222) sPDG=3212;                     // Sigma+ -> Sigma0
        else if(pPDG==3212) sPDG=3112;                     // Sigma0 -> Sigma-
        else if(pPDG==3322) sPDG=3312;                     // Xi0    -> Xi-
    }
    else
    {
        G4cout<<"Error:G4QuasiElRatios::ChExer: NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
        G4Exception("G4QuasiElRatios::ChExer:","21",FatalException,"QE complain");
        //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
    }
    if(sPDG) mS=mNeut;
    else
    {
        G4cout<<"Error:G4QuasiElRatios::ChExer: BAD pPDG="<<pPDG<<", NPDG="<<NPDG<<G4endl;
        G4Exception("G4QuasiElRatios::ChExer:","21",FatalException,"QE complain");
        //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
    }
    G4double mT2=mT*mT;
    G4double mS2=mS*mS;
    G4double E=(tot4M.m2()-mT2-mS2)/(mT+mT);
    G4double E2=E*E;
    if(E<0. || E2<mS2)
    {
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
    G4double P=std::sqrt(E2-mS2);                   // Momentum in pseudo laboratory system
    // @@ Temporary NN t-dependence for all hadrons
    G4int PDG=2212;                                                // *TMP* instead of pPDG
    if(pPDG==2112||pPDG==-211||pPDG==-321) PDG=2112;               // *TMP* instead of pPDG
    if(!Z && N==1)                 // Change for Quasi-Elastic on neutron
    {
        Z=1;
        N=0;
        if     (PDG==2212) PDG=2112;
        else if(PDG==2112) PDG=2212;
    }
    G4double xSec=0.;                        // Prototype of Recalculated Cross Section *TMP*
    if(PDG==2212) xSec=PCSmanager->GetChipsCrossSection(P, Z, N, PDG); // P CrossSect *TMP*
    else          xSec=NCSmanager->GetChipsCrossSection(P, Z, N, PDG); // N CrossSect *TMP*
    // @@ check a possibility to separate p, n, or alpha (!)
    if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
    {
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); //Do Nothing Action
    }
    G4double mint=0.;                        // Prototype of functional rand -t (MeV^2) *TMP*
    if(PDG==2212) mint=PCSmanager->GetExchangeT(Z,N,PDG);// P functional rand -t(MeV^2) *TMP*
    else          mint=NCSmanager->GetExchangeT(Z,N,PDG);// N functional rand -t(MeV^2) *TMP*
    G4double maxt=0.;                                    // Prototype of max possible -t
    if(PDG==2212) maxt=PCSmanager->GetHMaxT();           // max possible -t
    else          maxt=NCSmanager->GetHMaxT();           // max possible -t
    G4double cost=1.-mint/maxt;                          // cos(theta) in CMS
    if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
    {
        if     (cost>1.)  cost=1.;
        else if(cost<-1.) cost=-1.;
        else
        {
            G4cerr<<"G4QuasiFreeRatio::ChExer:*NAN* c="<<cost<<",t="<<mint<<",tm="<<maxt<<G4endl;
            return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
        }
    }
    G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
    pr4M=G4LorentzVector(0.,0.,0.,mS);                        // 4mom of the scattered hadron
    G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
    if(!RelDecayIn2(tot4M, pr4M, reco4M, dir4M, cost, cost))
    {
        G4cerr<<"G4QFR::ChEx:t="<<tot4M<<tot4M.m()<<",mT="<<mT<<",mP="<<mS<<G4endl;
        //G4Exception("G4QFR::ChExer:","009",FatalException,"Decay of ElasticComp");
        return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
    return std::make_pair(reco4M*megaelectronvolt,pr4M*megaelectronvolt); // Result
} // End of ChExer

// Calculate ChEx/El ratio (p is in independent units, (Z,N) is target, pPDG is projectile)
G4double G4QuasiElRatios::ChExElCoef(G4double p, G4int Z, G4int N, G4int pPDG) 
{
    
    p/=MeV;                                // Converted from independent units
    G4double A=Z+N;
    if(A<1.5) return 0.;
    G4double Cex=0.;
    if     (pPDG==2212) Cex=N/(A+Z);
    else if(pPDG==2112) Cex=Z/(A+N);
    else G4cout<<"*Warning*G4CohChrgExchange::ChExElCoef: wrong PDG="<<pPDG<<G4endl;
    Cex*=Cex;                         // Coherent processes squares the amplitude
    // @@ This is true only for nucleons: other projectiles must be treated differently
    G4double sp=std::sqrt(p);
    G4double p2=p*p;            
    G4double p4=p2*p2;
    G4double dl1=std::log(p)-5.;
    G4double T=(6.75+.14*dl1*dl1+13./p)/(1.+.14/p4)+.6/(p4+.00013);
    G4double U=(6.25+8.33e-5/p4/p)*(p*sp+.34)/p2/p; 
    G4double R=U/T;
    return Cex*R*R;
}

// Decay of Hadron In2Particles f&s, f is in respect to the direction of HadronMomentumDir
G4bool G4QuasiElRatios::RelDecayIn2(G4LorentzVector& theMomentum, G4LorentzVector& f4Mom, G4LorentzVector& s4Mom,
                                     G4LorentzVector& dir, G4double maxCost, G4double minCost)
{
    G4double fM2 = f4Mom.m2();
    G4double fM  = std::sqrt(fM2);              // Mass of the 1st Hadron
    G4double sM2 = s4Mom.m2();
    G4double sM  = std::sqrt(sM2);              // Mass of the 2nd Hadron
    G4double iM2 = theMomentum.m2();
    G4double iM  = std::sqrt(iM2);              // Mass of the decaying hadron
    G4double vP  = theMomentum.rho();      // Momentum of the decaying hadron
    G4double dE  = theMomentum.e();        // Energy of the decaying hadron
    if(dE<vP)
    {
        G4cerr<<"***G4QHad::RelDecIn2: Tachionic 4-mom="<<theMomentum<<", E-p="<<dE-vP<<G4endl;
        G4double accuracy=.000001*vP;
        G4double emodif=std::fabs(dE-vP);
        //if(emodif<accuracy)
        //{
        G4cerr<<"G4QHadron::RelDecIn2: *Boost* E-p shift is corrected to "<<emodif<<G4endl;
        theMomentum.setE(vP+emodif+.01*accuracy);
        //}
    }
    G4ThreeVector ltb = theMomentum.boostVector();// Boost vector for backward Lorentz Trans.
    G4ThreeVector ltf = -ltb;              // Boost vector for forward Lorentz Trans.
    G4LorentzVector cdir = dir;            // A copy to make a transformation to CMS
    cdir.boost(ltf);                       // Direction transpormed to CMS of the Momentum
    G4ThreeVector vdir = cdir.vect();      // 3-Vector of the direction-particle
    G4ThreeVector vx(0.,0.,1.);            // Ort in the direction of the reference particle
    G4ThreeVector vy(0.,1.,0.);            // First ort orthogonal to the direction
    G4ThreeVector vz(1.,0.,0.);            // Second ort orthoganal to the direction
    if(vdir.mag2() > 0.)                   // the refference particle isn't at rest in CMS
    {
        vx = vdir.unit();                    // Ort in the direction of the reference particle
        G4ThreeVector vv= vx.orthogonal();   // Not normed orthogonal vector (!)
        vy = vv.unit();                      // First ort orthogonal to the direction
        vz = vx.cross(vy);                   // Second ort orthoganal to the direction
    }
    if(maxCost> 1.) maxCost= 1.;
    if(minCost<-1.) minCost=-1.;
    if(maxCost<-1.) maxCost=-1.;
    if(minCost> 1.) minCost= 1.;
    if(minCost> maxCost) minCost=maxCost;
    if(std::fabs(iM-fM-sM)<.00000001)
    {
        G4double fR=fM/iM;
        G4double sR=sM/iM;
        f4Mom=fR*theMomentum;
        s4Mom=sR*theMomentum;
        return true;
    }
    else if (iM+.001<fM+sM || iM==0.)
    {//@@ Later on make a quark content check for the decay
        G4cerr<<"***G4QH::RelDecIn2: fM="<<fM<<"+sM="<<sM<<">iM="<<iM<<",d="<<iM-fM-sM<<G4endl;
        return false;
    }
    G4double d2 = iM2-fM2-sM2;
    G4double p2 = (d2*d2/4.-fM2*sM2)/iM2;    // Decay momentum(^2) in CMS of Quasmon
    if(p2<0.)
    {
        p2=0.;
    }
    G4double p  = std::sqrt(p2);
    G4double ct = maxCost;
    if(maxCost>minCost)
    {
        G4double dcost=maxCost-minCost;
        ct = minCost+dcost*G4UniformRand();
    }
    G4double phi= twopi*G4UniformRand();  // @@ Change 360.*deg to M_TWOPI (?)
    G4double ps=0.;
    if(std::fabs(ct)<1.) ps = p * std::sqrt(1.-ct*ct);
    else
    {
        if(ct>1.) ct=1.;
        if(ct<-1.) ct=-1.;
    }
    G4ThreeVector pVect=(ps*std::sin(phi))*vz+(ps*std::cos(phi))*vy+p*ct*vx;
    
    f4Mom.setVect(pVect);
    f4Mom.setE(std::sqrt(fM2+p2));
    s4Mom.setVect((-1)*pVect);
    s4Mom.setE(std::sqrt(sM2+p2));
    
    if(f4Mom.e()+.001<f4Mom.rho())G4cerr<<"*G4QH::RDIn2:*Boost* f4M="<<f4Mom<<",e-p="
        <<f4Mom.e()-f4Mom.rho()<<G4endl;
    f4Mom.boost(ltb);                        // Lor.Trans. of 1st hadron back to LS
    if(s4Mom.e()+.001<s4Mom.rho())G4cerr<<"*G4QH::RDIn2:*Boost* s4M="<<s4Mom<<",e-p="
        <<s4Mom.e()-s4Mom.rho()<<G4endl;
    s4Mom.boost(ltb);                        // Lor.Trans. of 2nd hadron back to LS
    return true;
} // End of "RelDecayIn2"






