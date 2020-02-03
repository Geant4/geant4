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
/// \file electromagnetic/TestEm17/src/MuCrossSections.cc
/// \brief Implementation of the MuCrossSections class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MuCrossSections.hh"

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;
 
MuCrossSections::MuCrossSections()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuCrossSections::~MuCrossSections()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MuCrossSections::CR_Macroscopic(const G4String& process, 
                                         const G4Material* material,
                                         G4double tkin, G4double ep)
                                          
// return the macroscopic cross section (1/L) in GEANT4 internal units
{
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* NbOfAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();

  G4double SIGMA = 0.;
  G4int nelm = material->GetNumberOfElements(); 
  for (G4int i=0; i < nelm; ++i)
  {
    const G4Element* element = (*theElementVector)[i];
    SIGMA += NbOfAtomsPerVolume[i] * CR_PerAtom(process, element, tkin, ep);
  }
  return SIGMA;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MuCrossSections::CR_PerAtom(const G4String& process, 
                                     const G4Element* element,
                                     G4double tkin, G4double ep)
{
 G4double z = element->GetZ();
 G4double a = element->GetA();

 G4double sigma = 0.;
 if (process == "muBrems")
   sigma = CRB_Mephi(z,a/(g/mole),tkin/GeV,ep/GeV)*(cm2/(g*GeV))*a/Avogadro;
   
 else if (process == "muIoni")
   sigma = CRK_Mephi(z,a/(g/mole),tkin/GeV,ep/GeV)*(cm2/(g*GeV))*a/Avogadro;
   
 //else if (process == "muNucl")
 else if (process == "muonNuclear")
   sigma = CRN_Mephi(z,a/(g/mole),tkin/GeV,ep/GeV)*(cm2/(g*GeV))*a/Avogadro;
   
 else if (process == "muPairProd")
   sigma = CRP_Mephi(z,a/(g/mole),tkin/GeV,ep/GeV)*(cm2/(g*GeV))*a/Avogadro;
   
 return sigma;        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MuCrossSections::CRB_Mephi(double z,double a,double tkin,double ep)

//***********************************************************************
//***        crb_g4_1.inc        in comparison with crb_.inc, following
//***        changes are introduced (September 24th, 1998):
//***                1) function crb_g4 (Z,A,Tkin,EP), Tkin is kinetic energy
//***                2) special case of hidrogen (Z<1.5; b,b1,Dn_star)
//***                Numerical comparison: 5 decimal digits coincide.
//***
//***        Cross section for bremsstrahlung by fast muon
//***        By R.P.Kokoulin, September 1998
//***        Formulae from Kelner,Kokoulin,Petrukhin 1995, Preprint MEPhI
//***        (7,18,19,20,21,25,26); Dn (18) is modified to incorporate
//***        Bugaev's inelatic nuclear correction (28) for Z > 1.
//***********************************************************************
{
//    G4double Z,A,Tkin,EP;
    G4double crb_g4;
    G4double e,v,delta,rab0,z_13,dn,b,b1,dn_star,rab1,fn,epmax1,fe,rab2;
//    
    G4double ame=0.51099907e-3; // GeV
    G4double lamu=0.105658389;        // GeV
    G4double re=2.81794092e-13; // cm
    G4double avno=6.022137e23;
    G4double alpha=1./137.036;
    G4double rmass=lamu/ame; // "207"
    G4double coeff=16./3.*alpha*avno*(re/rmass)*(re/rmass); // cm^2
    G4double sqrte=1.64872127; // sqrt(2.71828...)
    G4double btf=183.;
    G4double btf1=1429.;
    G4double bh=202.4;
    G4double bh1=446.;
//***
        if(ep >= tkin)
        {
          crb_g4=0.;
          return crb_g4;
        }
        e=tkin+lamu;
        v=ep/e;
        delta=lamu*lamu*v/(2.*(e-ep)); // qmin
        rab0=delta*sqrte;
        z_13=pow(z,-0.3333333); //
//***                nuclear size and excitation, screening parameters
        dn=1.54*pow(a,0.27);
        if(z <= 1.5) // special case for hydrogen
        {
          b=bh;
          b1=bh1;
          dn_star=dn;
        }
        else
        {
          b=btf;
          b1=btf1;
          dn_star=pow(dn,(1.-1./z)); // with Bugaev's correction
        }
//***                nucleus contribution logarithm
        rab1=b*z_13;
        fn=log(rab1/(dn_star*(ame+rab0*rab1))*(lamu+delta*(dn_star*sqrte-2.)));
        if(fn < 0.) fn=0.;
//***                electron contribution logarithm
        epmax1=e/(1.+lamu*rmass/(2.*e));
        if(ep >= epmax1)
        {
          fe=0.;
          goto label10;
        }
        rab2=b1*z_13*z_13;
        fe=log(rab2*lamu/((1.+delta*rmass/(ame*sqrte))*(ame+rab0*rab2)));
        if(fe < 0.) fe=0.;
//***
label10:
        crb_g4=coeff*(1.-v*(1.-0.75*v))*z*(z*fn+fe)/(a*ep);
        return crb_g4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MuCrossSections::CRK_Mephi(G4double z,G4double a,G4double tkin,G4double ep)
 
//***********************************************************************
//***        Cross section for knock-on electron production by fast muons
//***        (including bremsstrahlung e-diagrams and rad. correction).
//***        Units: cm^2/(g*GeV); Tkin, ep - GeV.
//***        By R.P.Kokoulin, October 1998
//***        Formulae from Kelner,Kokoulin,Petrukhin, Phys.Atom.Nuclei, 1997
//***        (a bit simplified Kelner's version of Eq.30 - with 2 logarithms).
//***
{
//    G4double Z,A,Tkin,EP;
    G4double crk_g4;
    G4double v,sigma0,a1,a3;
//
    G4double ame=0.51099907e-3; // GeV
    G4double lamu=0.105658389; // GeV
    G4double re=2.81794092e-13; // cm
    G4double avno=6.022137e23;
    G4double alpha=1./137.036;
    G4double lpi=3.141592654;
    G4double bmu=lamu*lamu/(2.*ame);
    G4double coeff0=avno*2.*lpi*ame*re*re;
    G4double coeff1=alpha/(2.*lpi);
//***

    G4double e=tkin+lamu;
    if(e < 0.) {
      G4cout << "CRK: " << tkin << "  " << ep << "  " << e << G4endl;
    }
    G4double epmax=e/(1.+bmu/e);
    if(ep >= epmax)
        {
          crk_g4=0.;
          return crk_g4;
        }
    v=ep/e;
    sigma0=coeff0*(z/a)*(1.-ep/epmax+0.5*v*v)/(ep*ep);
    a1=std::log(1.+2.*ep/ame);
    a3=std::log(4.*e*(e-ep)/(lamu*lamu));
    crk_g4=sigma0*(1.+coeff1*a1*(a3-a1));
    return crk_g4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MuCrossSections::CRN_Mephi(G4double /* z */,G4double a,G4double tkin,G4double ep)

//***********************************************************************
//***        Differential cross section for photonuclear muon interaction.
//***        Formulae from Borog & Petrukhin, 1975
//***        Real photon cross section: Caldwell e.a., 1979
//***        Nuclear shadowing: Brodsky e.a., 1972
//***        Units: cm^2 / g GeV.
//***        CRN_G4_1.inc        January 31st, 1998        R.P.Kokoulin
//***********************************************************************
{
//    G4double Z,A,Tkin,EP;
    G4double crn_g4;
    G4double e,aeff,sigph,v,v1,v2,amu2,up,down;
//***
    G4double lamu=0.105658389; // GeV
    G4double avno=6.022137e23;
    G4double amp=0.9382723; // GeV
    G4double lpi=3.14159265;
    G4double alpha=1./137.036;
//***
    G4double epmin_phn=0.20; // GeV
    G4double alam2=0.400000; // GeV**2
    G4double alam =0.632456; // sqrt(alam2)
    G4double coeffn=alpha/lpi*avno*1e-30; // cm^2/microbarn
//***
        e=tkin+lamu;
        crn_g4=0.;
        if(ep >= e-0.5*amp) return crn_g4;
        if(ep <= epmin_phn) return crn_g4;
        aeff=0.22*a+0.78*pow(a,0.89); // shadowing
        sigph=49.2+11.1*log(ep)+151.8/sqrt(ep); // microbarn
        v=ep/e;
        v1=1.-v;
        v2=v*v;
        amu2=lamu*lamu;
        up=e*e*v1/amu2*(1.+amu2*v2/(alam2*v1));
        down=1.+ep/alam*(1.+alam/(2.*amp)+ep/alam);
        crn_g4=coeffn*aeff/a*sigph/ep*(-v1+(v1+0.5*v2*(1.+2.*amu2/alam2))*log(up/down));
        if(crn_g4 < 0.) crn_g4=0.;
        return crn_g4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MuCrossSections::CRP_Mephi(G4double z,G4double a,G4double tkin,G4double ep)

//**********************************************************************
//***        crp_g4_1.inc        in comparison with crp_m.inc, following
//***        changes are introduced (January 16th, 1998):
//***                1) Avno/A, cm^2/gram GeV
//***                2) zeta_loss(E,Z)         from Kelner 1997, Eqs.(53-54)
//***                3) function crp_g4 (Z,A,Tkin,EP), Tkin is kinetic energy
//***                4) bbb=183        (Thomas-Fermi)
//***                5) special case of hidrogen (Z<1.5), bbb,g1,g2
//***                6) expansions in 'xi' are simplified        (Jan.17th,1998)
//***
//***        Cross section for electron pair production by fast muon
//***        By R.P.Kokoulin, December 1997
//***        Formulae from Kokoulin & Petrukhin 1971, Hobart, Eqs.(1,8,9,10)
{
//    G4double Z,A,Tkin,EP;
    G4double crp_g4;
    G4double bbbtf,bbbh,g1tf,g2tf,g1h,g2h,e,z13,e1,alf,a3,bbb;
    G4double g1,g2,zeta1,zeta2,zeta,z2,screen0,a0,a1,bet,xi0,del;
    G4double tmn,sum,a4,a5,a6,a7,a9,xi,xii,xi1,screen,yeu,yed,ye1;
    G4double ale,cre,be,fe,ymu,ymd,ym1,alm_crm,a10,bm,fm;
//
    G4double ame=0.51099907e-3; // GeV
    G4double lamu=0.105658389; // GeV
    G4double re=2.81794092e-13; // cm
    G4double avno=6.022137e23;
    G4double lpi=3.14159265;
    G4double alpha=1./137.036;
    G4double rmass=lamu/ame; // "207"
    G4double coeff=4./(3.*lpi)*(alpha*re)*(alpha*re)*avno; // cm^2
    G4double sqrte=1.64872127; // sqrt(2.71828...)
    G4double c3=3.*sqrte*lamu/4.; // for limits
    G4double c7=4.*ame; // -"-
    G4double c8=6.*lamu*lamu; // -"-

    G4double xgi[8]={.0199,.1017,.2372,.4083,.5917,.7628,.8983,.9801}; // Gauss, 8
    G4double wgi[8]={.0506,.1112,.1569,.1813,.1813,.1569,.1112,.0506}; // Gauss, 8
    bbbtf=183.; // for the moment...
    bbbh=202.4; // for the moment...
    g1tf=1.95e-5;
    g2tf=5.3e-5;
    g1h=4.4e-5;
    g2h=4.8e-5;

        e=tkin+lamu;
        z13=pow(z,0.3333333);
        e1=e-ep;
        crp_g4=0.;
        if(e1 <= c3*z13) return crp_g4; // ep > max
        alf=c7/ep; // 4m/ep
        a3=1.-alf;
        if(a3 <= 0.) return crp_g4; // ep < min
//***                zeta calculation
        if(z <= 1.5) // special case of hidrogen
        {
          bbb=bbbh;
          g1=g1h;
          g2=g2h;
        }
        else
        {
          bbb=bbbtf;
          g1=g1tf;
          g2=g2tf;
        }
        zeta1=0.073*log(e/(lamu+g1*z13*z13*e))-0.26;
        if(zeta1 > 0.)
        {
          zeta2=0.058*log(e/(lamu+g2*z13*e))-0.14;
          zeta=zeta1/zeta2;
        }
        else
        {
          zeta=0.;
        }
        z2=z*(z+zeta); //
//***                just to check (for comparison with crp_m)
//        z2=z*(z+1.)
//        bbb=189.
//***
        screen0=2.*ame*sqrte*bbb/(z13*ep); // be careful with "ame"
        a0=e*e1;
        a1=ep*ep/a0; // 2*beta
        bet=0.5*a1; // beta
        xi0=0.25*rmass*rmass*a1; // xi0
        del=c8/a0; // 6mu^2/EE'
        tmn=log((alf+2.*del*a3)/(1.+(1.-del)*sqrt(a3))); // log(1-rmax)
        sum=0.;
        for(int i=0; i<=7; i++) // integration
        {
          a4=exp(tmn*xgi[i]); // 1-r
          a5=a4*(2.-a4); // 1-r2
          a6=1.-a5; // r2
          a7=1.+a6; // 1+r2
          a9=3.+a6; // 3+r2
          xi=xi0*a5;
          xii=1./xi;
          xi1=1.+xi;
          screen=screen0*xi1/a5;
          yeu=5.-a6+4.*bet*a7;
          yed=2.*(1.+3.*bet)*log(3.+xii)-a6-a1*(2.-a6);
          ye1=1.+yeu/yed;
          ale=log(bbb/z13*sqrt(xi1*ye1)/(1.+screen*ye1));
          cre=0.5*log(1.+(1.5/rmass*z13)*(1.5/rmass*z13)*xi1*ye1);
          if(xi <= 1e3) //
          {
            be=((2.+a6)*(1.+bet)+xi*a9)*log(1.+xii)+(a5-bet)/xi1-a9;
          }
          else
          {
            be=(3.-a6+a1*a7)/(2.*xi); // -(6.-5.*a6+3.*bet*a6)/(6.*xi*xi);
          }
          fe=max(0.,(ale-cre)*be);          
          ymu=4.+a6+3.*bet*a7;
          ymd=a7*(1.5+a1)*log(3.+xi)+1.-1.5*a6;
          ym1=1.+ymu/ymd;
          alm_crm=log(bbb*rmass/(1.5*z13*z13*(1.+screen*ym1)));
          if(xi >= 1e-3) //
          {
            a10=(1.+a1)*a5; // (1+2b)(1-r2)
            bm=(a7*(1.+1.5*bet)-a10*xii)*log(xi1)+xi*(a5-bet)/xi1+a10;
          }
          else
          {
            bm=(5.-a6+bet*a9)*(xi/2.); // -(11.-5.*a6+.5*bet*(5.+a6))*(xi*xi/6.)
          }
          fm=max(0.,(alm_crm)*bm);          
//***
          sum=sum+a4*(fe+fm/(rmass*rmass))*wgi[i];
        }
        crp_g4=-tmn*sum*(z2/a)*coeff*e1/(e*ep);
        return crp_g4;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


