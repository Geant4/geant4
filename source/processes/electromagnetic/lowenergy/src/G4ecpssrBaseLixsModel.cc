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

#include <iostream>
#include "G4ecpssrBaseLixsModel.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4LinLogInterpolation.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ecpssrBaseLixsModel::G4ecpssrBaseLixsModel()
{
  verboseLevel=0;

  // Storing FLi data needed for 0.2 to 3.0  velocities region
  const char* path = G4FindDataDir("G4LEDATA");

  if (!path) {
    G4Exception("G4ecpssrLCrossSection::G4ecpssrBaseLixsModel()","em0006", FatalException ,"G4LEDATA environment variable not set");
    return;
  }
  std::ostringstream fileName1;
  std::ostringstream fileName2;

  fileName1 << path << "/pixe/uf/FL1.dat";
  fileName2 << path << "/pixe/uf/FL2.dat";

  // Reading of FL1.dat
  std::ifstream FL1(fileName1.str().c_str());
  if (!FL1) G4Exception("G4ecpssrLCrossSection::G4ecpssrBaseLixsModel()","em0003",FatalException, "error opening FL1 data file");

  dummyVec1.push_back(0.);

  while(!FL1.eof())
    {
      double x1;
      double y1;

      FL1>>x1>>y1;

      //  Mandatory vector initialization
      if (x1 != dummyVec1.back())
        {
          dummyVec1.push_back(x1);
          aVecMap1[x1].push_back(-1.);
        }

      FL1>>FL1Data[x1][y1];

      if (y1 != aVecMap1[x1].back()) aVecMap1[x1].push_back(y1);
    }

  // Reading of FL2.dat

  std::ifstream FL2(fileName2.str().c_str());
  if (!FL2) G4Exception("G4ecpssrLCrossSection::G4ecpssrBaseLixsModel()","em0003", FatalException," error opening FL2 data file");

  dummyVec2.push_back(0.);
  
  while(!FL2.eof())
    {
      double x2;
      double y2;

      FL2>>x2>>y2;

      //  Mandatory vector initialization
      if (x2 != dummyVec2.back())
        {
          dummyVec2.push_back(x2);
          aVecMap2[x2].push_back(-1.);
        }

      FL2>>FL2Data[x2][y2];

      if (y2 != aVecMap2[x2].back()) aVecMap2[x2].push_back(y2);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ecpssrBaseLixsModel::~G4ecpssrBaseLixsModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::ExpIntFunction(G4int n,G4double x)

{
// this function allows fast evaluation of the n order exponential integral function En(x)
  G4int i;
  G4int ii;
  G4int nm1;
  G4double a;
  G4double b;
  G4double c;
  G4double d;
  G4double del;
  G4double fact;
  G4double h;
  G4double psi;
  G4double ans = 0;
  static const G4double euler= 0.5772156649;
  static const G4int maxit= 100;
  static const G4double fpmin = 1.0e-30;
  static const G4double eps = 1.0e-7;
  nm1=n-1;
  if (n<0 || x<0.0 || (x==0.0 && (n==0 || n==1)))
  G4cout << "*** WARNING in G4ecpssrBaseLixsModel::ExpIntFunction: bad arguments in ExpIntFunction"
         << G4endl;
  else {
    if (n==0) ans=G4Exp(-x)/x;
    else {
      if (x==0.0) ans=1.0/nm1;
      else {
	if (x > 1.0) {
	  b=x+n;
	  c=1.0/fpmin;
	  d=1.0/b;
	  h=d;
	  for (i=1;i<=maxit;i++) {
	    a=-i*(nm1+i);
	    b +=2.0;
	    d=1.0/(a*d+b);
	    c=b+a/c;
	    del=c*d;
	    h *=del;
	    if (std::fabs(del-1.0) < eps) {
	      ans=h*G4Exp(-x);
	      return ans;
	    }
	  }
	} else {
	  ans = (nm1!=0 ? 1.0/nm1 : -std::log(x)-euler);
	  fact=1.0;
	  for (i=1;i<=maxit;i++) {
	    fact *=-x/i;
	    if (i !=nm1) del = -fact/(i-nm1);
	    else {
	      psi = -euler;
	      for (ii=1;ii<=nm1;ii++) psi +=1.0/ii;
	      del=fact*(-std::log(x)+psi);
	    }
	    ans += del;
	    if (std::fabs(del) < std::fabs(ans)*eps) return ans;
	  }
	}
      }
    }
  }
  return ans;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::CalculateL1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  if (zTarget <=4) return 0.;

  //this L1-CrossSection calculation method is done according to Werner Brandt and Grzegorz Lapicki, Phys.Rev.A20 N2 (1979),
  //and using data tables of O. Benka et al. At.Data Nucl.Data Tables Vol.22 No.3 (1978).

  G4NistManager* massManager = G4NistManager::Instance();

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double  zIncident = 0;
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )
    {
      zIncident = (aProtone->GetPDGCharge())/eplus;
    }
  else
    {
      if (massIncident == aAlpha->GetPDGMass()) 
	{
	  zIncident  = (aAlpha->GetPDGCharge())/eplus;
	}
      else
	{
	  G4cout << "*** WARNING in G4ecpssrBaseLixsModel::CalculateL1CrossSection : Proton or Alpha incident particles only. " << G4endl;
	  G4cout << massIncident << ", " << aAlpha->GetPDGMass() << " (alpha)" << aProtone->GetPDGMass() << " (proton)" << G4endl;
	  return 0;
	}
    }

  G4double l1BindingEnergy = transitionManager->Shell(zTarget,1)->BindingEnergy(); //Observed binding energy of L1-subshell
  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;
  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2; //Mass of the system (projectile, target)
  static const G4double zlshell= 4.15;
  // *** see Benka, ADANDT 22, p 223
  G4double screenedzTarget = zTarget-zlshell; //Effective nuclear charge as seen by electrons in L1-sub shell
  static const G4double rydbergMeV= 13.6056923e-6;

  static const G4double nl= 2.;
  // *** see Benka, ADANDT 22, p 220, f3
  G4double tetal1 = (l1BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV); //Screening parameter
  // *** see Benka, ADANDT 22, p 220, f3

  if (verboseLevel>0) G4cout << "  tetal1=" <<  tetal1<< G4endl;

  G4double reducedEnergy = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget);
  // *** also called etaS
  // *** see Benka, ADANDT 22, p 220, f3

  static const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ; //Bohr radius of hydrogen

  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*std::pow(screenedzTarget,-4.);
  // *** see Benka, ADANDT 22, p 220, f2, for protons
  // *** see Basbas, Phys Rev A7, p 1000

  G4double velocityl1 = CalculateVelocity(1, zTarget, massIncident, energyIncident); // Scaled velocity

  if (verboseLevel>0) G4cout << "  velocityl1=" << velocityl1<< G4endl;

  static const G4double l1AnalyticalApproximation= 1.5;
  G4double x1 =(nl*l1AnalyticalApproximation)/velocityl1;
  // *** 1.5 is cK = cL1 (it is 1.25 for L2 & L3)
  // *** see Brandt, Phys Rev A20, p 469, f16 in expression of h

  if (verboseLevel>0) G4cout << "  x1=" << x1<< G4endl;

  G4double electrIonizationEnergyl1=0.;
  // *** see Basbas, Phys Rev A17, p1665, f27
  // *** see Brandt, Phys Rev A20, p469
  // *** see Liu, Comp Phys Comm 97, p325, f A5

  if ( x1<=0.035)  electrIonizationEnergyl1= 0.75*pi*(std::log(1./(x1*x1))-1.);
  else
    {
      if ( x1<=3.)
        electrIonizationEnergyl1 =G4Exp(-2.*x1)/(0.031+(0.213*std::pow(x1,0.5))+(0.005*x1)-(0.069*std::pow(x1,3./2.))+(0.324*x1*x1));
      else
	{if ( x1<=11.) electrIonizationEnergyl1 =2.*G4Exp(-2.*x1)/std::pow(x1,1.6);}
    }

  G4double hFunctionl1 =(electrIonizationEnergyl1*2.*nl)/(tetal1*std::pow(velocityl1,3)); //takes into account the polarization effect
  // *** see Brandt, Phys Rev A20, p 469, f16

  if (verboseLevel>0) G4cout << "  hFunctionl1=" << hFunctionl1<< G4endl;

  G4double gFunctionl1 = (1.+(9.*velocityl1)+(31.*velocityl1*velocityl1)+(49.*std::pow(velocityl1,3.))+(162.*std::pow(velocityl1,4.))+(63.*std::pow(velocityl1,5.))+(18.*std::pow(velocityl1,6.))+(1.97*std::pow(velocityl1,7.)))/std::pow(1.+velocityl1,9.);//takes into account the reduced binding effect
  // *** see Brandt, Phys Rev A20, p 469, f19

  if (verboseLevel>0) G4cout << "  gFunctionl1=" << gFunctionl1<< G4endl;

  G4double sigmaPSS_l1 = 1.+(((2.*zIncident)/(screenedzTarget*tetal1))*(gFunctionl1-hFunctionl1)); //Binding-polarization factor
  // *** also called dzeta
  // *** also called epsilon
  // *** see Basbas, Phys Rev A17, p1667, f45

  if (verboseLevel>0) G4cout << "sigmaPSS_l1 =" << sigmaPSS_l1<< G4endl;

  const G4double cNaturalUnit= 137.;

  G4double yl1Formula=0.4*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(nl*velocityl1/sigmaPSS_l1);
  // *** also called yS
  // *** see Brandt, Phys Rev A20, p467, f6
  // *** see Brandt, Phys Rev A23, p1728

  G4double l1relativityCorrection = std::pow((1.+(1.1*yl1Formula*yl1Formula)),0.5)+yl1Formula; // Relativistic correction parameter
  // *** also called mRS
  // *** see Brandt, Phys Rev A20, p467, f6

  //G4double reducedVelocity_l1 = velocityl1*std::pow(l1relativityCorrection,0.5); //Reduced velocity parameter

  G4double L1etaOverTheta2;

  G4double  universalFunction_l1 = 0.;

  G4double sigmaPSSR_l1;

  // low velocity formula
  // *****************
  if ( velocityl1 <20.  )
  {

    L1etaOverTheta2 =(reducedEnergy* l1relativityCorrection)/((tetal1*sigmaPSS_l1)*(tetal1*sigmaPSS_l1));
    // *** 1) RELATIVISTIC CORRECTION ADDED
    // *** 2) sigma_PSS_l1 ADDED
    // *** reducedEnergy is etaS, l1relativityCorrection is mRS
    // *** see Phys Rev A20, p468, top

    if ( ((tetal1*sigmaPSS_l1) >=0.2) && ((tetal1*sigmaPSS_l1) <=2.6670) && (L1etaOverTheta2>=0.1e-3) && (L1etaOverTheta2<=0.866e2) )
      universalFunction_l1 = FunctionFL1((tetal1*sigmaPSS_l1),L1etaOverTheta2);

    if (verboseLevel>0) G4cout << "at low velocity range, universalFunction_l1  =" << universalFunction_l1 << G4endl;

    sigmaPSSR_l1 = (sigma0/(tetal1*sigmaPSS_l1))*universalFunction_l1;// Plane-wave Born -Aproximation L1-subshell ionisation Cross Section
  // *** see Benka, ADANDT 22, p220, f1

    if (verboseLevel>0) G4cout << "  at low velocity range, sigma PWBA L1 CS  = " << sigmaPSSR_l1<< G4endl;
  }
  else
  {
    L1etaOverTheta2 = reducedEnergy/(tetal1*tetal1);
    // Medium & high velocity
    // *** 1) NO RELATIVISTIC CORRECTION
    // *** 2) NO sigma_PSS_l1
    // *** see Benka, ADANDT 22, p223

    if ( (tetal1 >=0.2) && (tetal1 <=2.6670) && (L1etaOverTheta2>=0.1e-3) && (L1etaOverTheta2<=0.866e2) )
      universalFunction_l1 = FunctionFL1(tetal1,L1etaOverTheta2);

    if (verboseLevel>0) G4cout << "at medium and high velocity range, universalFunction_l1  =" << universalFunction_l1 << G4endl;

    sigmaPSSR_l1 = (sigma0/tetal1)*universalFunction_l1;// Plane-wave Born -Aproximation L1-subshell ionisation Cross Section
  // *** see Benka, ADANDT 22, p220, f1

    if (verboseLevel>0) G4cout << "  sigma PWBA L1 CS at medium and high velocity range = " << sigmaPSSR_l1<< G4endl;
  }

  G4double pssDeltal1 = (4./(systemMass *sigmaPSS_l1*tetal1))*(sigmaPSS_l1/velocityl1)*(sigmaPSS_l1/velocityl1);
  // *** also called dzeta*delta
  // *** see Brandt, Phys Rev A23, p1727, f B2

  if (verboseLevel>0) G4cout << "  pssDeltal1=" << pssDeltal1<< G4endl;

  if (pssDeltal1>1) return 0.;

  G4double energyLossl1 = std::pow(1-pssDeltal1,0.5);
  // *** also called z
  // *** see Brandt, Phys Rev A23, p1727, after f B2

  if (verboseLevel>0) G4cout << "  energyLossl1=" << energyLossl1<< G4endl;

  G4double coulombDeflectionl1 =
   (8.*pi*zIncident/systemMass)*std::pow(tetal1*sigmaPSS_l1,-2.)*std::pow(velocityl1/sigmaPSS_l1,-3.)*(zTarget/screenedzTarget);
  // *** see Brandt, Phys Rev A20, v2s and f2 and B2
  // *** with factor n2 compared to Brandt, Phys Rev A23, p1727, f B3

  G4double cParameterl1 =2.* coulombDeflectionl1/(energyLossl1*(energyLossl1+1.));
  // *** see Brandt, Phys Rev A23, p1727, f B4

  G4double coulombDeflectionFunction_l1 = 9.*ExpIntFunction(10,cParameterl1); //Coulomb-deflection effect correction

  if (verboseLevel>0) G4cout << "  coulombDeflectionFunction_l1 =" << coulombDeflectionFunction_l1 << G4endl;

  G4double crossSection_L1 = coulombDeflectionFunction_l1 * sigmaPSSR_l1;

  //ECPSSR L1 -subshell cross section is estimated at perturbed-stationnairy-state(PSS)
  //and reduced by the energy-loss(E),the Coulomb deflection(C),and the relativity(R) effects

  if (verboseLevel>0) G4cout << "  crossSection_L1 =" << crossSection_L1 << G4endl;

  if (crossSection_L1 >= 0)
  {
    return crossSection_L1 * barn;
  }

  else {return 0;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::CalculateL2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)

{
  if (zTarget <=13 ) return 0.;

  // this L2-CrossSection calculation method is done according to Werner Brandt and Grzegorz Lapicki, Phys.Rev.A20 N2 (1979),
  // and using data tables of O. Benka et al. At.Data Nucl.Data Tables Vol.22 No.3 (1978).

  G4NistManager* massManager = G4NistManager::Instance();

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double  zIncident = 0;

  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )
   zIncident = (aProtone->GetPDGCharge())/eplus;

  else
    {
      if (massIncident == aAlpha->GetPDGMass())
	zIncident  = (aAlpha->GetPDGCharge())/eplus;

      else
	{
	  G4cout << "*** WARNING in G4ecpssrBaseLixsModel::CalculateL2CrossSection : Proton or Alpha incident particles only. " << G4endl;
	  G4cout << massIncident << ", " << aAlpha->GetPDGMass() << " (alpha)" << aProtone->GetPDGMass() << " (proton)" << G4endl;
	  return 0;
	}
    }

  G4double l2BindingEnergy = transitionManager->Shell(zTarget,2)->BindingEnergy(); //Observed binding energy of L2-subshell

  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;

  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2; //Mass of the system (projectile, target)

  const G4double zlshell= 4.15;

  G4double screenedzTarget = zTarget-zlshell; //Effective nuclear charge as seen by electrons in L2-subshell

  const G4double rydbergMeV= 13.6056923e-6;

  const G4double nl= 2.;

  G4double tetal2 = (l2BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV); //Screening parameter

  if (verboseLevel>0) G4cout << "  tetal2=" <<  tetal2<< G4endl;

  G4double reducedEnergy = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget);

  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ; //Bohr radius of hydrogen

  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*std::pow(screenedzTarget,-4.);

  G4double velocityl2 = CalculateVelocity(2,  zTarget, massIncident, energyIncident); // Scaled velocity

  if (verboseLevel>0) G4cout << "  velocityl2=" << velocityl2<< G4endl;

  const G4double l23AnalyticalApproximation= 1.25;

  G4double x2 = (nl*l23AnalyticalApproximation)/velocityl2;

  if (verboseLevel>0) G4cout << "  x2=" << x2<< G4endl;

  G4double electrIonizationEnergyl2=0.;

  if ( x2<=0.035)  electrIonizationEnergyl2= 0.75*pi*(std::log(1./(x2*x2))-1.);
  else
    {
      if ( x2<=3.)
        electrIonizationEnergyl2 =G4Exp(-2.*x2)/(0.031+(0.213*std::pow(x2,0.5))+(0.005*x2)-(0.069*std::pow(x2,3./2.))+(0.324*x2*x2));
      else
        {if ( x2<=11.) electrIonizationEnergyl2 =2.*G4Exp(-2.*x2)/std::pow(x2,1.6);  }
    }

  G4double hFunctionl2 =(electrIonizationEnergyl2*2.*nl)/(tetal2*std::pow(velocityl2,3)); //takes into account  the polarization effect

  if (verboseLevel>0) G4cout << "  hFunctionl2=" << hFunctionl2<< G4endl;

  G4double gFunctionl2 = (1.+(10.*velocityl2)+(45.*velocityl2*velocityl2)+(102.*std::pow(velocityl2,3.))+(331.*std::pow(velocityl2,4.))+(6.7*std::pow(velocityl2,5.))+(58.*std::pow(velocityl2,6.))+(7.8*std::pow(velocityl2,7.))+ (0.888*std::pow(velocityl2,8.)) )/std::pow(1.+velocityl2,10.);
  //takes into account the reduced binding effect

  if (verboseLevel>0) G4cout << "  gFunctionl2=" << gFunctionl2<< G4endl;

  G4double sigmaPSS_l2 = 1.+(((2.*zIncident)/(screenedzTarget*tetal2))*(gFunctionl2-hFunctionl2)); //Binding-polarization factor

  if (verboseLevel>0) G4cout << "  sigmaPSS_l2=" << sigmaPSS_l2<< G4endl;

  const G4double cNaturalUnit= 137.;

  G4double yl2Formula=0.15*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocityl2/sigmaPSS_l2);

  G4double l2relativityCorrection = std::pow((1.+(1.1*yl2Formula*yl2Formula)),0.5)+yl2Formula; // Relativistic correction parameter

  G4double L2etaOverTheta2;

  G4double  universalFunction_l2 = 0.;

  G4double sigmaPSSR_l2 ;

  if ( velocityl2 < 20. )
  {

    L2etaOverTheta2 = (reducedEnergy*l2relativityCorrection)/((sigmaPSS_l2*tetal2)*(sigmaPSS_l2*tetal2));

    if ( (tetal2*sigmaPSS_l2>=0.2) && (tetal2*sigmaPSS_l2<=2.6670) && (L2etaOverTheta2>=0.1e-3) && (L2etaOverTheta2<=0.866e2) )
      universalFunction_l2 = FunctionFL2((tetal2*sigmaPSS_l2),L2etaOverTheta2);

    sigmaPSSR_l2 = (sigma0/(tetal2*sigmaPSS_l2))*universalFunction_l2;

    if (verboseLevel>0) G4cout << "  sigma PWBA L2 CS at low velocity range = " << sigmaPSSR_l2<< G4endl;
  }
  else
  {

    L2etaOverTheta2 = reducedEnergy /(tetal2*tetal2);

    if ( (tetal2>=0.2) && (tetal2<=2.6670) && (L2etaOverTheta2>=0.1e-3) && (L2etaOverTheta2<=0.866e2) )
      universalFunction_l2 = FunctionFL2((tetal2),L2etaOverTheta2);

    sigmaPSSR_l2 = (sigma0/tetal2)*universalFunction_l2;

    if (verboseLevel>0) G4cout << "  sigma PWBA L2 CS at medium and high velocity range = " << sigmaPSSR_l2<< G4endl;

  }

  G4double pssDeltal2 = (4./(systemMass*sigmaPSS_l2*tetal2))*(sigmaPSS_l2/velocityl2)*(sigmaPSS_l2/velocityl2);

  if (pssDeltal2>1) return 0.;

  G4double energyLossl2 = std::pow(1-pssDeltal2,0.5);

    if (verboseLevel>0) G4cout << "  energyLossl2=" << energyLossl2<< G4endl;

  G4double coulombDeflectionl2
    =(8.*pi*zIncident/systemMass)*std::pow(tetal2*sigmaPSS_l2,-2.)*std::pow(velocityl2/sigmaPSS_l2,-3.)*(zTarget/screenedzTarget);

  G4double cParameterl2 = 2.*coulombDeflectionl2/(energyLossl2*(energyLossl2+1.));

  G4double coulombDeflectionFunction_l2 = 11.*ExpIntFunction(12,cParameterl2); //Coulomb-deflection effect correction
  // *** see Brandt, Phys Rev A10, p477, f25

  if (verboseLevel>0) G4cout << "  coulombDeflectionFunction_l2 =" << coulombDeflectionFunction_l2 << G4endl;

  G4double crossSection_L2 = coulombDeflectionFunction_l2 * sigmaPSSR_l2;
  //ECPSSR L2 -subshell cross section is estimated at perturbed-stationnairy-state(PSS)
  //and reduced by the energy-loss(E),the Coulomb deflection(C),and the relativity(R) effects

  if (verboseLevel>0) G4cout << "  crossSection_L2 =" << crossSection_L2 << G4endl;

  if (crossSection_L2 >= 0)
  {
     return crossSection_L2 * barn;
  }
  else {return 0;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4ecpssrBaseLixsModel::CalculateL3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)

{
  if (zTarget <=13) return 0.;

  //this L3-CrossSection calculation method is done according to Werner Brandt and Grzegorz Lapicki, Phys.Rev.A20 N2 (1979),
  //and using data tables of O. Benka et al. At.Data Nucl.Data Tables Vol.22 No.3 (1978).

  G4NistManager* massManager = G4NistManager::Instance();

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double  zIncident = 0;

  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )

   zIncident = (aProtone->GetPDGCharge())/eplus;

  else
    {
      if (massIncident == aAlpha->GetPDGMass())

	  zIncident  = (aAlpha->GetPDGCharge())/eplus;

      else
	{
	  G4cout << "*** WARNING in G4ecpssrBaseLixsModel::CalculateL3CrossSection : Proton or Alpha incident particles only. " << G4endl;
	  G4cout << massIncident << ", " << aAlpha->GetPDGMass() << " (alpha)" << aProtone->GetPDGMass() << " (proton)" << G4endl;
	  return 0;
	}
    }

  G4double l3BindingEnergy = transitionManager->Shell(zTarget,3)->BindingEnergy();

  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;

  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2;//Mass of the system (projectile, target)

  const G4double zlshell= 4.15;

  G4double screenedzTarget = zTarget-zlshell;//Effective nuclear charge as seen by electrons in L3-subshell

  const G4double rydbergMeV= 13.6056923e-6;

  const G4double nl= 2.;

  G4double tetal3 = (l3BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV);//Screening parameter

    if (verboseLevel>0) G4cout << "  tetal3=" <<  tetal3<< G4endl;

  G4double reducedEnergy = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget);

  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;//Bohr radius of hydrogen

  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*std::pow(screenedzTarget,-4.);

  G4double velocityl3 = CalculateVelocity(3, zTarget, massIncident, energyIncident);// Scaled velocity

    if (verboseLevel>0) G4cout << "  velocityl3=" << velocityl3<< G4endl;

  const G4double l23AnalyticalApproximation= 1.25;

  G4double x3 = (nl*l23AnalyticalApproximation)/velocityl3;

    if (verboseLevel>0) G4cout << "  x3=" << x3<< G4endl;

  G4double electrIonizationEnergyl3=0.;

  if ( x3<=0.035)  electrIonizationEnergyl3= 0.75*pi*(std::log(1./(x3*x3))-1.);
    else
    {
      if ( x3<=3.) electrIonizationEnergyl3 =G4Exp(-2.*x3)/(0.031+(0.213*std::pow(x3,0.5))+(0.005*x3)-(0.069*std::pow(x3,3./2.))+(0.324*x3*x3));
      else
	{
	  if ( x3<=11.) electrIonizationEnergyl3 =2.*G4Exp(-2.*x3)/std::pow(x3,1.6);}
    }

  G4double hFunctionl3 =(electrIonizationEnergyl3*2.*nl)/(tetal3*std::pow(velocityl3,3));//takes into account the polarization effect

    if (verboseLevel>0) G4cout << "  hFunctionl3=" << hFunctionl3<< G4endl;

  G4double gFunctionl3 = (1.+(10.*velocityl3)+(45.*velocityl3*velocityl3)+(102.*std::pow(velocityl3,3.))+(331.*std::pow(velocityl3,4.))+(6.7*std::pow(velocityl3,5.))+(58.*std::pow(velocityl3,6.))+(7.8*std::pow(velocityl3,7.))+ (0.888*std::pow(velocityl3,8.)) )/std::pow(1.+velocityl3,10.);
  //takes into account the reduced binding effect

    if (verboseLevel>0) G4cout << "  gFunctionl3=" << gFunctionl3<< G4endl;

  G4double sigmaPSS_l3 = 1.+(((2.*zIncident)/(screenedzTarget*tetal3))*(gFunctionl3-hFunctionl3));//Binding-polarization factor

    if (verboseLevel>0) G4cout << "sigmaPSS_l3 =" << sigmaPSS_l3<< G4endl;

  const G4double cNaturalUnit= 137.;

  G4double yl3Formula=0.15*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocityl3/sigmaPSS_l3);

  G4double l3relativityCorrection = std::pow((1.+(1.1*yl3Formula*yl3Formula)),0.5)+yl3Formula; // Relativistic correction parameter

  G4double L3etaOverTheta2;

  G4double  universalFunction_l3 = 0.;

  G4double sigmaPSSR_l3;

  if ( velocityl3 < 20. )
  {

    L3etaOverTheta2 = (reducedEnergy* l3relativityCorrection)/((sigmaPSS_l3*tetal3)*(sigmaPSS_l3*tetal3));

    if ( (tetal3*sigmaPSS_l3>=0.2) && (tetal3*sigmaPSS_l3<=2.6670) && (L3etaOverTheta2>=0.1e-3) && (L3etaOverTheta2<=0.866e2) )

      universalFunction_l3 = 2.*FunctionFL2((tetal3*sigmaPSS_l3), L3etaOverTheta2  );

    sigmaPSSR_l3 = (sigma0/(tetal3*sigmaPSS_l3))*universalFunction_l3;

    if (verboseLevel>0) G4cout << "  sigma PWBA L3 CS at low velocity range = " << sigmaPSSR_l3<< G4endl;

  }

  else

  {

    L3etaOverTheta2 = reducedEnergy/(tetal3*tetal3);

    if ( (tetal3>=0.2) && (tetal3<=2.6670) && (L3etaOverTheta2>=0.1e-3) && (L3etaOverTheta2<=0.866e2) )

      universalFunction_l3 = 2.*FunctionFL2(tetal3, L3etaOverTheta2  );

    sigmaPSSR_l3 = (sigma0/tetal3)*universalFunction_l3;

      if (verboseLevel>0) G4cout << "  sigma PWBA L3 CS at medium and high velocity range = " << sigmaPSSR_l3<< G4endl;
  }

  G4double pssDeltal3 = (4./(systemMass*sigmaPSS_l3*tetal3))*(sigmaPSS_l3/velocityl3)*(sigmaPSS_l3/velocityl3);

    if (verboseLevel>0) G4cout << "  pssDeltal3=" << pssDeltal3<< G4endl;

  if (pssDeltal3>1) return 0.;

  G4double energyLossl3 = std::pow(1-pssDeltal3,0.5);

  if (verboseLevel>0) G4cout << "  energyLossl3=" << energyLossl3<< G4endl;

  G4double coulombDeflectionl3 =
    (8.*pi*zIncident/systemMass)*std::pow(tetal3*sigmaPSS_l3,-2.)*std::pow(velocityl3/sigmaPSS_l3,-3.)*(zTarget/screenedzTarget);

  G4double cParameterl3 = 2.*coulombDeflectionl3/(energyLossl3*(energyLossl3+1.));

  G4double coulombDeflectionFunction_l3 = 11.*ExpIntFunction(12,cParameterl3);//Coulomb-deflection effect correction
  // *** see Brandt, Phys Rev A10, p477, f25

    if (verboseLevel>0) G4cout << "  coulombDeflectionFunction_l3 =" << coulombDeflectionFunction_l3 << G4endl;

  G4double crossSection_L3 =  coulombDeflectionFunction_l3 * sigmaPSSR_l3;
  //ECPSSR L3 -subshell cross section is estimated at perturbed-stationnairy-state(PSS)
  //and reduced by the energy-loss(E),the Coulomb deflection(C),and the relativity(R) effects

    if (verboseLevel>0) G4cout << "  crossSection_L3 =" << crossSection_L3 << G4endl;

  if (crossSection_L3 >= 0)
  {
    return crossSection_L3 * barn;
  }
  else {return 0;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::CalculateVelocity(G4int subShell, G4int zTarget, G4double massIncident,  G4double energyIncident)

{

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double liBindingEnergy = transitionManager->Shell(zTarget,subShell)->BindingEnergy();

  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (!((massIncident == aProtone->GetPDGMass()) || (massIncident == aAlpha->GetPDGMass())))
    {
      G4cout << "*** WARNING in G4ecpssrBaseLixsModel::CalculateVelocity : Proton or Alpha incident particles only. " << G4endl;
      G4cout << massIncident << ", " << aAlpha->GetPDGMass() << " (alpha)" << aProtone->GetPDGMass() << " (proton)" << G4endl;
      return 0;
    }

  constexpr G4double zlshell= 4.15;

  G4double screenedzTarget = zTarget- zlshell;

  constexpr G4double rydbergMeV= 13.6056923e-6;

  constexpr G4double nl= 2.;

  G4double tetali = (liBindingEnergy*nl*nl)/(screenedzTarget*screenedzTarget*rydbergMeV);

  G4double reducedEnergy = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget);

  G4double velocity = 2.*nl*std::pow(reducedEnergy,0.5)/tetali;
  // *** see Brandt, Phys Rev A10, p10, f4

  return velocity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::FunctionFL1(G4double k, G4double theta)
{
  G4double sigma = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;
  G4double xs12 = 0;
  G4double xs21 = 0;
  G4double xs22 = 0;

  // PROTECTION TO ALLOW INTERPOLATION AT MINIMUM AND MAXIMUM Eta/Theta2 values

  if (
       theta==8.66e-4 ||
       theta==8.66e-3 ||
       theta==8.66e-2 ||
       theta==8.66e-1 ||
       theta==8.66e+00 ||
       theta==8.66e+01
  ) theta=theta-1e-12;

  if (
       theta==1.e-4 ||
       theta==1.e-3 ||
       theta==1.e-2 ||
       theta==1.e-1 ||
       theta==1.e+00 ||
       theta==1.e+01
  ) theta=theta+1e-12;

  // END PROTECTION

  auto t2 = std::upper_bound(dummyVec1.begin(),dummyVec1.end(), k);
  auto t1 = t2-1;

  auto e12 = std::upper_bound(aVecMap1[(*t1)].begin(),aVecMap1[(*t1)].end(), theta);
  auto e11 = e12-1;

  auto e22 = std::upper_bound(aVecMap1[(*t2)].begin(),aVecMap1[(*t2)].end(), theta);
  auto e21 = e22-1;

  valueT1  =*t1;
  valueT2  =*t2;
  valueE21 =*e21;
  valueE22 =*e22;
  valueE12 =*e12;
  valueE11 =*e11;

  xs11 = FL1Data[valueT1][valueE11];
  xs12 = FL1Data[valueT1][valueE12];
  xs21 = FL1Data[valueT2][valueE21];
  xs22 = FL1Data[valueT2][valueE22];

  if (verboseLevel>0)
    G4cout
    << valueT1 << " "
    << valueT2 << " "
    << valueE11 << " "
    << valueE12 << " "
    << valueE21 << " "
    << valueE22 << " "
    << xs11 << " "
    << xs12 << " "
    << xs21 << " "
    << xs22 << " "
    << G4endl;

  G4double xsProduct = xs11 * xs12 * xs21 * xs22;

  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);

  if (xsProduct != 0.)
  {
    sigma = QuadInterpolator(  valueE11, valueE12,
    			       valueE21, valueE22,
			       xs11, xs12,
			       xs21, xs22,
			       valueT1, valueT2,
			       k, theta );
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::FunctionFL2(G4double k, G4double theta)
{

  G4double sigma = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;
  G4double xs12 = 0;
  G4double xs21 = 0;
  G4double xs22 = 0;

  // PROTECTION TO ALLOW INTERPOLATION AT MINIMUM AND MAXIMUM Eta/Theta2 values

  if (
       theta==8.66e-4 ||
       theta==8.66e-3 ||
       theta==8.66e-2 ||
       theta==8.66e-1 ||
       theta==8.66e+00 ||
       theta==8.66e+01
  ) theta=theta-1e-12;

  if (
       theta==1.e-4 ||
       theta==1.e-3 ||
       theta==1.e-2 ||
       theta==1.e-1 ||
       theta==1.e+00 ||
       theta==1.e+01
  ) theta=theta+1e-12;

  // END PROTECTION

  auto t2 = std::upper_bound(dummyVec2.begin(),dummyVec2.end(), k);
  auto t1 = t2-1;
  auto e12 = std::upper_bound(aVecMap2[(*t1)].begin(),aVecMap2[(*t1)].end(), theta);
  auto e11 = e12-1;
  auto e22 = std::upper_bound(aVecMap2[(*t2)].begin(),aVecMap2[(*t2)].end(), theta);
  auto e21 = e22-1;

  valueT1  =*t1;
  valueT2  =*t2;
  valueE21 =*e21;
  valueE22 =*e22;
  valueE12 =*e12;
  valueE11 =*e11;

  xs11 = FL2Data[valueT1][valueE11];
  xs12 = FL2Data[valueT1][valueE12];
  xs21 = FL2Data[valueT2][valueE21];
  xs22 = FL2Data[valueT2][valueE22];

  if (verboseLevel>0)
    G4cout
    << valueT1 << " "
    << valueT2 << " "
    << valueE11 << " "
    << valueE12 << " "
    << valueE21 << " "
    << valueE22 << " "
    << xs11 << " "
    << xs12 << " "
    << xs21 << " "
    << xs22 << " "
    << G4endl;

  G4double xsProduct = xs11 * xs12 * xs21 * xs22;

  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);

  if (xsProduct != 0.)
  {
    sigma = QuadInterpolator(  valueE11, valueE12,
    			       valueE21, valueE22,
			       xs11, xs12,
			       xs21, xs22,
			       valueT1, valueT2,
			       k, theta );
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::LinLinInterpolate(G4double e1,
						        G4double e2,
						        G4double e,
						        G4double xs1,
						        G4double xs2)
{
  G4double value = xs1 + (xs2 - xs1)*(e - e1)/ (e2 - e1);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::LinLogInterpolate(G4double e1,
						        G4double e2,
						        G4double e,
						        G4double xs1,
						        G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = G4Exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::LogLogInterpolate(G4double e1,
						        G4double e2,
						        G4double e,
						        G4double xs1,
						        G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseLixsModel::QuadInterpolator(G4double e11, G4double e12,
						       G4double e21, G4double e22,
						       G4double xs11, G4double xs12,
						       G4double xs21, G4double xs22,
						       G4double t1, G4double t2,
						       G4double t, G4double e)
{
  // Log-Log
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;

}
