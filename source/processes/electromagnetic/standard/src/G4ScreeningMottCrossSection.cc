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
//      G4ScreeningMottCrossSection.cc
//-------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4ScreeningMottCrossSection
//
// Author:      Cristina Consolandi
//
// Creation date: 20.10.2011  
//
// Modifications:
// 27-05-2012 Added Analytic Fitting to the Mott Cross Section by means of G4MottCoefficients class.
//
//
// Class Description:
//	Computation of electron Coulomb Scattering Cross Section.
//	Suitable for high energy electrons and light target materials. 
//
//      Reference:
//      M.J. Boschini et al.
//     "Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//      Proc. of the 13th International Conference on Particle Physics and Advanced Technology 
//      (13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//	Available at: http://arxiv.org/abs/1111.4042v4
//
//      1) Mott Differential Cross Section Approximation:  
//	   For Target material up to Z=92 (U):
//         As described in http://arxiv.org/abs/1111.4042v4 
//         par. 2.1 , eq. (16)-(17)
//         Else (Z>92):
//	   W. A. McKinley and H. Fashbach, Phys. Rev. 74, (1948) 1759.
//      2) Screening coefficient: 
//      vomn G. Moliere, Z. Naturforsh A2 (1947), 133-145; A3 (1948), 78.
//      3) Nuclear Form Factor: 
//      A.V. Butkevich et al. Nucl. Instr. and Meth. in Phys. Res. A 488 (2002), 282-294.
//
// -------------------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ScreeningMottCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MottCoefficients.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4LossTableManager.hh"
#include "G4NucleiProperties.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


using namespace std;

G4ScreeningMottCrossSection::G4ScreeningMottCrossSection():
   cosThetaMin(1.0),
   cosThetaMax(-1.0),
   alpha(fine_structure_const),
   htc2(hbarc_squared),
   e2(electron_mass_c2*classic_electr_radius) 
{
  TotalCross=0;

  fNistManager = G4NistManager::Instance();
  particle=0;

  spin = mass = mu_rel=0;
  tkinLab = momLab2 = invbetaLab2=0;
  tkin = mom2 = invbeta2=beta=gamma=0;

  Trec=targetZ = targetMass = As =0;
  etag = ecut = 0.0;

  targetA = 0;

  cosTetMinNuc=0;
  cosTetMaxNuc=0;

  for(G4int i=0 ; i<5; i++){
    for(G4int j=0; j< 6; j++){
      coeffb[i][j]=0;
    }
  }

  mottcoeff = new G4MottCoefficients();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ScreeningMottCrossSection::~G4ScreeningMottCrossSection()
{
  delete mottcoeff;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ScreeningMottCrossSection::Initialise(const G4ParticleDefinition* p,
                                          G4double CosThetaLim)
{
  SetupParticle(p);
  tkin = targetZ = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  particle = p;
  cosThetaMin = CosThetaLim; 

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4ScreeningMottCrossSection::SetScreeningCoefficient()
{
  G4double alpha2=alpha*alpha;
  //Bohr radius 
  G4double a0=  Bohr_radius  ;//0.529e-8*cm;
  //Thomas-Fermi screening length
  G4double aU=0.88534*a0/pow(targetZ,1./3.);
  G4double twoR2=aU*aU;

  G4double factor= 1.13 + 3.76*targetZ*targetZ*invbeta2*alpha2;
  As=0.25*(htc2)/(twoR2*mom2)*factor;
  //cout<<"0k .........................As  "<<As<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::GetScreeningAngle()
{
  SetScreeningCoefficient();

  //cout<<" As "<<As<<endl;
  if(As < 0.0) { As = 0.0; }
  else if(As > 1.0) { As = 1.0; }

  G4double screenangle=2.*asin(sqrt(As));
  //	cout<<"  screenangle  "<<  screenangle <<endl;
  if(screenangle>=pi) screenangle=pi;
	
  return screenangle;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ScreeningMottCrossSection::SetupKinematic(G4double ekin, G4double Z )
{
  //...Target
  G4int iz = G4int(Z);
  G4double A = fNistManager->GetAtomicMassAmu(iz);
  G4int ia = G4int(A);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(ia, iz);

  targetZ = Z;
  targetA = fNistManager->GetAtomicMassAmu(iz);
  targetMass= mass2;

  mottcoeff->SetMottCoeff(targetZ, coeffb);

  //cout<<"......... targetA "<< targetA <<endl;
  //cout<<"......... targetMass "<< targetMass/MeV <<endl;

  // incident particle lab
  tkinLab = ekin;
  momLab2 = tkinLab*(tkinLab + 2.0*mass);
  invbetaLab2 = 1.0 +  mass*mass/momLab2;

  G4double etot = tkinLab + mass;
  G4double ptot = sqrt(momLab2);
  G4double m12  = mass*mass;
                
  // relativistic reduced mass from publucation
  // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
        
  //incident particle & target nucleus
  G4double Ecm=sqrt(m12 + mass2*mass2 + 2.0*etot*mass2);
  mu_rel=mass*mass2/Ecm;
  G4double momCM= ptot*mass2/Ecm;
  // relative system
  mom2 = momCM*momCM;
  invbeta2 = 1.0 +  mu_rel*mu_rel/mom2;
  tkin = momCM*sqrt(invbeta2) - mu_rel;//Ekin of mu_rel
  G4double  beta2=1./invbeta2;
  beta=std::sqrt(beta2) ;
  G4double gamma2= 1./(1.-beta2);
  gamma=std::sqrt(gamma2);

  //.........................................................

  G4double screenangle=GetScreeningAngle()/10.;
  //cout<<" screenangle [rad] "<<screenangle/rad <<endl;

  cosTetMinNuc =min( cosThetaMin ,cos(screenangle));
  cosTetMaxNuc =cosThetaMax;
	
  //cout<<"ok..............mu_rel "<<mu_rel<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::FormFactor2ExpHof(G4double angle)
{
  G4double M=targetMass; 
  G4double E=tkinLab;
  G4double Etot=E+mass;
  G4double Tmax=2.*M*E*(E+2.*mass)/(mass*mass+M*M+2.*M*Etot);
  G4double T=Tmax*pow(sin(angle/2.),2.);
  G4double q2=T*(T+2.*M);
  q2/=htc2;//1/cm2
  G4double RN=1.27e-13*pow(targetA,0.27)*cm;
  G4double xN= (RN*RN*q2);
  G4double den=(1.+xN/12.);
  G4double FN=1./(den*den);
  G4double form2=(FN*FN);

  return form2;

  //cout<<"..................... form2 "<< form2<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::McFcorrection(G4double angle )
{
  G4double  beta2=1./invbeta2;
  G4double sintmezzi=std::sin(angle/2.);
  G4double sin2tmezzi = sintmezzi*sintmezzi;
  G4double R=1.-beta2*sin2tmezzi + targetZ*alpha*beta*pi*sintmezzi*(1.-sintmezzi);
  return R;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4ScreeningMottCrossSection::RatioMottRutherford(G4double angle)
{
  G4double R=0;
  G4double fcost=std::sqrt((1. -cos(angle)));
  G4double a[5];
  G4double shift=0.7181228;
  G4double beta0= beta -shift;

  for(G4int j=0 ;j<=4;j++){
    a[j]=0;
  }

  for(G4int j=0 ;j<=4;j++){
    for(G4int k=0;k<=5;k++ ){  
      a[j]+=coeffb[j][k]*pow(beta0,k);
    }
  }

  for(G4int j=0 ;j<=4 ;j++){
    R+=a[j]* pow(fcost,j);
  }
  return R;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::NuclearCrossSection()
{
  if(cosTetMaxNuc >= cosTetMinNuc) return 0.0;

  TotalCross=0;

  G4double anglemin =std::acos(cosTetMinNuc);
  G4double anglemax =std::acos(cosTetMaxNuc); 

  static const G4double limit = 1.e-9;
  if(anglemin < limit) {
    anglemin = GetScreeningAngle()/10.;
    if(anglemin < limit) { anglemin = limit; }
  }

  //cout<<" anglemin  "<< anglemin <<endl;

  G4double loganglemin=log10(anglemin);
  G4double loganglemax=log10(anglemax);
  G4double logdangle=0.01;

  G4int bins=(G4int)((loganglemax-loganglemin)/logdangle);

  vector<G4double> angle;
  vector<G4double> tet;
  vector<G4double> dangle;
  vector<G4double> cross;

  for(G4int i=0; i<=bins; i++ ){
    tet.push_back(0);
    dangle.push_back(0);
    angle.push_back(pow(10.,loganglemin+logdangle*i));
    cross.push_back(0);
  }

  G4int  dim = tet.size();
  //cout<<"dim--- "<<dim<<endl;

  for(G4int i=0; i<dim;i++){

    if(i!=dim-1){
      dangle[i]=(angle[i+1]-angle[i]);
      tet[i]=(angle[i+1]+angle[i])/2.;
    }else if(i==dim-1){
      break;
    }

    G4double R=0;
    G4double F2=FormFactor2ExpHof(tet[i]);
			
    if (coeffb[0][0]!=0){
      //cout<<" Mott....targetZ "<< targetZ<<endl;	
      R=RatioMottRutherford(tet[i]);
    } else if (coeffb[0][0]==0){
      // cout<<" McF.... targetZ "<< targetZ<<endl;
      R=McFcorrection(tet[i]);
    }

    //cout<<"----------------- R "<<R<<" F2 "<<F2<<endl;
    //                cout<<"angle "<<tet[i] << " F2 "<<F2<<endl;

    G4double e4=e2*e2;
    G4double den=2.*As+2.*pow(sin(tet[i]/2.),2.);
    G4double func=1./(den*den);

    G4double fatt= targetZ/(mu_rel*gamma*beta*beta);
    G4double sigma=e4*fatt*fatt*func;
    cross[i]=F2*R*sigma;
    G4double pi2sintet=2.*pi*sin(tet[i]);

    TotalCross+=pi2sintet*cross[i]*dangle[i];
  }//end integral

  //cout<< "ok ......... TotalCross "<<TotalCross<<endl;
  return TotalCross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::AngleDistribution(G4double anglein)
{
  G4double total=TotalCross ;
  G4double fatt= e2*targetZ/(mu_rel*gamma*beta*beta);
  G4double fatt2=fatt*fatt;
  total/=fatt2;

  G4double R=0;
  if (coeffb[0][0]!=0){
    //   cout<<" Mott....targetZ "<< targetZ<<endl;      
    R=RatioMottRutherford(anglein);
  } else if (coeffb[0][0]==0){
    // cout<<" McF.... targetZ "<< targetZ<<endl;
    R=McFcorrection(anglein);
  }

  G4double y=2.*pi*sin(anglein)*R*FormFactor2ExpHof(anglein)/
    ((2*As+2.*pow(sin(anglein/2.),2.))*(2.*As+2.*pow(sin(anglein/2.),2.) ));

  return y/total;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::GetScatteringAngle()
{
  //cout<<" G4ScreeningMottCrossSection::SampleCosineTheta ............."<<endl;	
  if(cosTetMaxNuc >= cosTetMinNuc) return 0.0;

  G4double anglemin=std::acos(cosTetMinNuc);	
  G4double anglemax= std::acos(cosTetMaxNuc);

  static const G4double limit = 1.e-9;
  if(anglemin < limit) {
    anglemin = GetScreeningAngle()/10.;
    if(anglemin < limit) { anglemin = limit; }
  }

  //	cout<<"................ tkinLab  "<< G4BestUnit(tkinLab,"Energy") << " anglemin=  "<<anglemin<<endl;
  //cout<<"anglemax=  "<<anglemax<<endl;
  G4double r =G4UniformRand();

  G4double loganglemin=log10(anglemin);
  G4double loganglemax=log10(anglemax);
  G4double logdangle=0.01;

  G4int bins=(G4int)((loganglemax-loganglemin)/logdangle);

  std::vector<G4double> angle;
  std::vector<G4double> tet;
  std::vector<G4double> dangle;
  
  for(G4int i=0; i<=bins; i++ ){
    tet.push_back(0);
    dangle.push_back(0);
    angle.push_back(pow(10.,loganglemin+logdangle*i));
  }

  G4int  dim = tet.size();
  G4double scattangle=0;
  G4double y=0;
  G4double dy=0;
  G4double area=0;

  for(G4int i=0; i<dim;i++){
    
    if(i!=dim-1){
      dangle[i]=(angle[i+1]-angle[i]);
      tet[i]=(angle[i+1]+angle[i])/2.;
    }else if(i==dim-1){
      break;
    }

    y+=AngleDistribution(tet[i])*dangle[i];
    dy= y-area ;
    area=y;

    if(r >=y-dy && r<=y+dy ){	
      scattangle= angle[i] +G4UniformRand()*dangle[i];
      //cout<<"y "<<y <<" r  "<< r  << " y/r "<<y/r << endl;
      break;
    }			
  }
  return scattangle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4ScreeningMottCrossSection::GetNewDirection(){

  G4ThreeVector dir(0.0,0.0,1.0);
	
  G4double z1=GetScatteringAngle();
  
  G4double sint = sin(z1);
  G4double cost = sqrt(1.0 - sint*sint);
  G4double phi  = twopi* G4UniformRand();
  G4double dirx = sint*cos(phi);
  G4double diry = sint*sin(phi);
  G4double dirz = cost;

  //.......set Trc
  G4double etot=tkinLab+mass;
  G4double mass2=targetMass;
  Trec=(1.0 - cost)* mass2*(etot*etot - mass*mass )/
    (mass*mass + mass2*mass2+ 2.*mass2*etot);
       
  dir.set(dirx,diry,dirz);

  return dir;
}


