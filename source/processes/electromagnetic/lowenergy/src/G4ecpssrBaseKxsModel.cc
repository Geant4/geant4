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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <cmath>
#include <iostream>
#include "G4ecpssrBaseKxsModel.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ecpssrBaseKxsModel::G4ecpssrBaseKxsModel()
{
    verboseLevel=0;

    // Storing C coefficients for high velocity formula
    G4String fileC1("pixe/uf/c1");
    tableC1 = new G4CrossSectionDataSet(new G4SemiLogInterpolation, 1.,1.);

    G4String fileC2("pixe/uf/c2");
    tableC2 = new G4CrossSectionDataSet(new G4SemiLogInterpolation, 1.,1.);

    G4String fileC3("pixe/uf/c3");
    tableC3 = new G4CrossSectionDataSet(new G4SemiLogInterpolation, 1.,1.);

    // Storing FK data needed for medium velocities region
    const char* path = G4FindDataDir("G4LEDATA");

    if (!path) {
      G4Exception("G4ecpssrBaseKxsModel::G4ecpssrBaseKxsModel()", "em0006", FatalException,"G4LEDATA environment variable not set" );
      return;
    }

    std::ostringstream fileName;
    fileName << path << "/pixe/uf/FK.dat";
    std::ifstream FK(fileName.str().c_str());

    if (!FK)
      G4Exception("G4ecpssrBaseKxsModel::G4ecpssrBaseKxsModel()", "em0003", FatalException,"error opening FK data file" );

    dummyVec.push_back(0.);

    while(!FK.eof())
    {
	double x;
	double y;

	FK>>x>>y;

	//  Mandatory vector initialization
        if (x != dummyVec.back())
        {
          dummyVec.push_back(x);
          aVecMap[x].push_back(-1.);
        }

        FK>>FKData[x][y];

        if (y != aVecMap[x].back()) aVecMap[x].push_back(y);

    }

    tableC1->LoadData(fileC1);
    tableC2->LoadData(fileC2);
    tableC3->LoadData(fileC3);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void print (G4double elem)
{
  G4cout << elem << " ";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ecpssrBaseKxsModel::~G4ecpssrBaseKxsModel()
{
  delete tableC1;
  delete tableC2;
  delete tableC3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseKxsModel::ExpIntFunction(G4int n,G4double x)

{
// this "ExpIntFunction" function allows fast evaluation of the n order exponential integral function En(x)
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
  const G4double euler= 0.5772156649;
  const G4int maxit= 100;
  const G4double fpmin = 1.0e-30;
  const G4double eps = 1.0e-7;
  nm1=n-1;
  if (n<0 || x<0.0 || (x==0.0 && (n==0 || n==1))) {
    G4cout << "*** WARNING in G4ecpssrBaseKxsModel::ExpIntFunction: bad arguments in ExpIntFunction" << G4endl;
    G4cout << n << ", " << x << G4endl;
  }
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

G4double G4ecpssrBaseKxsModel::CalculateCrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)

{
  // this K-CrossSection calculation method is done according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)//
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
	  G4cout << "*** WARNING in G4ecpssrBaseKxsModel::CalculateCrossSection : we can treat only Proton or Alpha incident particles " << G4endl;
	  return 0;
	}
  }

    if (verboseLevel>0) G4cout << "  massIncident=" << massIncident<< G4endl;

  G4double kBindingEnergy = transitionManager->Shell(zTarget,0)->BindingEnergy();

    if (verboseLevel>0) G4cout << "  kBindingEnergy=" << kBindingEnergy/eV<< G4endl;

  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;

    if (verboseLevel>0) G4cout << "  massTarget=" <<  massTarget<< G4endl;

  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2; //the mass of the system (projectile, target)

    if (verboseLevel>0) G4cout << "  systemMass=" <<  systemMass<< G4endl;

  constexpr G4double zkshell= 0.3;
  // *** see Brandt, Phys Rev A23, p 1727

  G4double screenedzTarget = zTarget-zkshell; // screenedzTarget is the screened nuclear charge of the target
  // *** see Brandt, Phys Rev A23, p 1727

  constexpr G4double rydbergMeV= 13.6056923e-6;

  G4double tetaK = kBindingEnergy/((screenedzTarget*screenedzTarget)*rydbergMeV);  //tetaK denotes the reduced binding energy of the electron
  // *** see Rice, ADANDT 20, p 504, f 2

    if (verboseLevel>0) G4cout << "  tetaK=" <<  tetaK<< G4endl;

  G4double velocity =(2./(tetaK*screenedzTarget))*std::pow(((energyIncident*electron_mass_c2)/(massIncident*rydbergMeV)),0.5);
  // *** also called xiK
  // *** see Brandt, Phys Rev A23, p 1727
  // *** see Basbas, Phys Rev A17, p 1656, f4

    if (verboseLevel>0) G4cout << "  velocity=" << velocity<< G4endl;

  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;

    if (verboseLevel>0) G4cout << "  bohrPow2Barn=" <<  bohrPow2Barn<< G4endl;

  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*std::pow(screenedzTarget,-4.);     //sigma0 is the initial cross section of K shell at stable state
  // *** see Benka, ADANDT 22, p 220, f2, for protons
  // *** see Basbas, Phys Rev A7, p 1000

  if (verboseLevel>0) G4cout << "  sigma0=" <<  sigma0<< G4endl;

  const G4double kAnalyticalApproximation= 1.5;
  G4double x = kAnalyticalApproximation/velocity;
  // *** see Brandt, Phys Rev A23, p 1727
  // *** see Brandt, Phys Rev A20, p 469, f16 in expression of h

    if (verboseLevel>0) G4cout << "  x=" << x<< G4endl;

  G4double electrIonizationEnergy;
  // *** see Basbas, Phys Rev A17, p1665, f27
  // *** see Brandt, Phys Rev A20, p469
  // *** see Liu, Comp Phys Comm 97, p325, f A5

  if ((0.< x) && (x <= 0.035))
    {
      electrIonizationEnergy= 0.75*pi*(std::log(1./(x*x))-1.);
    }
  else
    {
      if ( (0.035 < x) && (x <=3.))
	{
	  electrIonizationEnergy =G4Exp(-2.*x)/(0.031+(0.213*std::pow(x,0.5))+(0.005*x)-(0.069*std::pow(x,3./2.))+(0.324*x*x));
	}

      else
	{
	  if ( (3.< x) && (x<=11.))
	   {
	     electrIonizationEnergy =2.*G4Exp(-2.*x)/std::pow(x,1.6);
	   }

	  else electrIonizationEnergy =0.;
	}
    }

    if (verboseLevel>0) G4cout << "  electrIonizationEnergy=" << electrIonizationEnergy<< G4endl;

  G4double hFunction =(electrIonizationEnergy*2.)/(tetaK*std::pow(velocity,3)); //hFunction represents the correction for polarization effet
  // *** see Brandt, Phys Rev A20, p 469, f16

    if (verboseLevel>0) G4cout << "  hFunction=" << hFunction<< G4endl;

  G4double gFunction = (1.+(9.*velocity)+(31.*velocity*velocity)+(98.*std::pow(velocity,3.))+(12.*std::pow(velocity,4.))+(25.*std::pow(velocity,5.))
			+(4.2*std::pow(velocity,6.))+(0.515*std::pow(velocity,7.)))/std::pow(1.+velocity,9.); //gFunction represents the correction for binding effet
  // *** see Brandt, Phys Rev A20, p 469, f19

    if (verboseLevel>0) G4cout << "  gFunction=" << gFunction<< G4endl;

  //-----------------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSS = 1.+(((2.*zIncident)/(screenedzTarget*tetaK))*(gFunction-hFunction)); //describes the perturbed stationnairy state of the affected atomic electon
  // *** also called dzeta
  // *** also called epsilon
  // *** see Basbas, Phys Rev A17, p1667, f45

    if (verboseLevel>0) G4cout << "  sigmaPSS=" << sigmaPSS<< G4endl;

    if (verboseLevel>0) G4cout << "  sigmaPSS*tetaK=" << sigmaPSS*tetaK<< G4endl;

  //----------------------------------------------------------------------------------------------------------------------------

  const G4double cNaturalUnit= 1/fine_structure_const;  // it's the speed of light according to Atomic-Unit-System

    if (verboseLevel>0) G4cout << "  cNaturalUnit=" << cNaturalUnit<< G4endl;

  G4double ykFormula=0.4*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocity/sigmaPSS);
  // *** also called yS
  // *** see Brandt, Phys Rev A20, p467, f6
  // *** see Brandt, Phys Rev A23, p1728

    if (verboseLevel>0) G4cout << "  ykFormula=" << ykFormula<< G4endl;

  G4double relativityCorrection = std::pow((1.+(1.1*ykFormula*ykFormula)),0.5)+ykFormula;// the relativistic correction parameter
  // *** also called mRS
  // *** see Brandt, Phys Rev A20, p467, f6

    if (verboseLevel>0) G4cout << "  relativityCorrection=" << relativityCorrection<< G4endl;

  G4double reducedVelocity = velocity*std::pow(relativityCorrection,0.5);  // presents the reduced collision velocity parameter
  // *** also called xiR
  // *** see Brandt, Phys Rev A20, p468, f7
  // *** see Brandt, Phys Rev A23, p1728

    if (verboseLevel>0) G4cout << "  reducedVelocity=" << reducedVelocity<< G4endl;

  G4double etaOverTheta2 = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget)
                           /(sigmaPSS*tetaK)/(sigmaPSS*tetaK);
  // *** see Benka, ADANDT 22, p220, f4 for eta
  // then we use sigmaPSS*tetaK == epsilon*tetaK

    if (verboseLevel>0) G4cout << "  etaOverTheta2=" << etaOverTheta2<< G4endl;

  G4double universalFunction = 0;

  // low velocity formula
  // *****************
   if ( velocity < 1. )
  // OR
  //if ( reducedVelocity/sigmaPSS < 1.)
  // *** see Brandt, Phys Rev A23, p1727
  // *** reducedVelocity/sigmaPSS is also called xiR/dzeta
  // *****************
    {
      if (verboseLevel>0) G4cout << "  Notice : FK is computed from low velocity formula" << G4endl;

      universalFunction = (std::pow(2.,9.)/45.)*std::pow(reducedVelocity/sigmaPSS,8.)*std::pow((1.+(1.72*(reducedVelocity/sigmaPSS)*(reducedVelocity/sigmaPSS))),-4.);// is the reduced universal cross section
      // *** see Brandt, Phys Rev A23, p1728

      if (verboseLevel>0) G4cout << "  universalFunction by Brandt 1981 =" << universalFunction<< G4endl;

    }
  else
  {

    if ( etaOverTheta2 > 86.6 && (sigmaPSS*tetaK) > 0.4 && (sigmaPSS*tetaK) < 2.9996 )
    {
      // High and medium energies. Method from Rice ADANDT 20, p506, 1977 on tables from Benka 1978

      if (verboseLevel>0) G4cout << "  Notice : FK is computed from high velocity formula" << G4endl;

      if (verboseLevel>0) G4cout << "  sigmaPSS*tetaK=" << sigmaPSS*tetaK << G4endl;

      G4double C1= tableC1->FindValue(sigmaPSS*tetaK);
      G4double C2= tableC2->FindValue(sigmaPSS*tetaK);
      G4double C3= tableC3->FindValue(sigmaPSS*tetaK);

        if (verboseLevel>0) G4cout << "  C1=" << C1 << G4endl;
        if (verboseLevel>0) G4cout << "  C2=" << C2 << G4endl;
        if (verboseLevel>0) G4cout << "  C3=" << C3 << G4endl;

      G4double etaK = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget);
      // *** see Benka, ADANDT 22, p220, f4 for eta

        if (verboseLevel>0) G4cout << "  etaK=" << etaK << G4endl;

      G4double etaT = (sigmaPSS*tetaK)*(sigmaPSS*tetaK)*(86.6); // at any theta, the largest tabulated etaOverTheta2 is 86.6
      // *** see Rice, ADANDT 20, p506

        if (verboseLevel>0) G4cout << "  etaT=" << etaT << G4endl;

      G4double fKT = FunctionFK((sigmaPSS*tetaK),86.6)*(etaT/(sigmaPSS*tetaK));
      // *** see Rice, ADANDT 20, p506

      if (FunctionFK((sigmaPSS*tetaK),86.6)<=0.)
	{
        G4cout <<
        "*** WARNING in G4ecpssrBaseKxsModel::CalculateCrossSection : unable to interpolate FK function in high velocity region ! ***" << G4endl;
	return 0;
	}

        if (verboseLevel>0) G4cout << "  FunctionFK=" << FunctionFK((sigmaPSS*tetaK),86.6) << G4endl;

        if (verboseLevel>0) G4cout << "  fKT=" << fKT << G4endl;

      G4double GK = C2/(4*etaK) + C3/(32*etaK*etaK);

	if (verboseLevel>0) G4cout << "  GK=" << GK << G4endl;

      G4double GT = C2/(4*etaT) + C3/(32*etaT*etaT);

        if (verboseLevel>0) G4cout << "  GT=" << GT << G4endl;

      G4double DT = fKT - C1*std::log(etaT) + GT;

        if (verboseLevel>0) G4cout << "  DT=" << DT << G4endl;

      G4double fKK = C1*std::log(etaK) + DT - GK;

        if (verboseLevel>0) G4cout << "  fKK=" << fKK << G4endl;

      G4double universalFunction3= fKK/(etaK/tetaK);
      // *** see Rice, ADANDT 20, p505, f7

        if (verboseLevel>0) G4cout << "  universalFunction3=" << universalFunction3 << G4endl;

      universalFunction=universalFunction3;

    }
    else if ( etaOverTheta2 >= 1.e-3 && etaOverTheta2 <= 86.6 && (sigmaPSS*tetaK) >= 0.4 && (sigmaPSS*tetaK) <= 2.9996 )
    {
      // From Benka 1978

      if (verboseLevel>0) G4cout << "  Notice : FK is computed from INTERPOLATED data" << G4endl;

      G4double universalFunction2 = FunctionFK((sigmaPSS*tetaK),etaOverTheta2);

      if (universalFunction2<=0)
      {
        G4cout <<
        "*** WARNING : G4ecpssrBaseKxsModel::CalculateCrossSection is unable to interpolate FK function in medium velocity region ! ***" << G4endl;
	return 0;
      }

      if (verboseLevel>0) G4cout << "  universalFunction2=" << universalFunction2 << " for theta=" << sigmaPSS*tetaK << " and etaOverTheta2=" << etaOverTheta2 << G4endl;

      universalFunction=universalFunction2;
    }

  }

  //----------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSSR = (sigma0/(sigmaPSS*tetaK))*universalFunction; //sigmaPSSR is the straight-line K-shell ionization cross section
  // *** see Benka, ADANDT 22, p220, f1

    if (verboseLevel>0) G4cout << "  sigmaPSSR=" << sigmaPSSR<< G4endl;

  //-----------------------------------------------------------------------------------------------------------------------

  G4double pssDeltaK = (4./(systemMass*sigmaPSS*tetaK))*(sigmaPSS/velocity)*(sigmaPSS/velocity);
  // *** also called dzetaK*deltaK
  // *** see Brandt, Phys Rev A23, p1727, f B2

    if (verboseLevel>0) G4cout << "  pssDeltaK=" << pssDeltaK<< G4endl;

  if (pssDeltaK>1) return 0.;

  G4double energyLoss = std::pow(1-pssDeltaK,0.5); //energyLoss incorporates the straight-line energy-loss
  // *** also called zK
  // *** see Brandt, Phys Rev A23, p1727, after f B2

    if (verboseLevel>0) G4cout << "  energyLoss=" << energyLoss<< G4endl;

  G4double energyLossFunction = (std::pow(2.,-9)/8.)*((((9.*energyLoss)-1.)*std::pow(1.+energyLoss,9.))+(((9.*energyLoss)+1.)*std::pow(1.-energyLoss,9.)));//energy loss function
  // *** also called fs
  // *** see Brandt, Phys Rev A23, p1718, f7

    if (verboseLevel>0) G4cout << "  energyLossFunction=" <<  energyLossFunction<< G4endl;

  //----------------------------------------------------------------------------------------------------------------------------------------------

  G4double coulombDeflection = (4.*pi*zIncident/systemMass)*std::pow(tetaK*sigmaPSS,-2.)*std::pow(velocity/sigmaPSS,-3.)*(zTarget/screenedzTarget); //incorporates Coulomb deflection parameter
  // *** see Brandt, Phys Rev A23, p1727, f B3

    if (verboseLevel>0) G4cout << "  cParameter-short=" << coulombDeflection<< G4endl;

  G4double cParameter = 2.*coulombDeflection/(energyLoss*(energyLoss+1.));
  // *** see Brandt, Phys Rev A23, p1727, f B4

    if (verboseLevel>0) G4cout << "  cParameter-full=" << cParameter<< G4endl;

  G4double coulombDeflectionFunction = 9.*ExpIntFunction(10,cParameter);                         //this function describes Coulomb-deflection effect
  // *** see Brandt, Phys Rev A23, p1727

    if (verboseLevel>0) G4cout << "  ExpIntFunction(10,cParameter) =" << ExpIntFunction(10,cParameter) << G4endl;

    if (verboseLevel>0) G4cout << "  coulombDeflectionFunction =" << coulombDeflectionFunction << G4endl;

  //--------------------------------------------------------------------------------------------------------------------------------------------------

  G4double crossSection =  0;

  crossSection = energyLossFunction* coulombDeflectionFunction*sigmaPSSR;  //this ECPSSR cross section is estimated at perturbed-stationnairy-state(PSS)
                                                                           //and it's reduced by the energy-loss(E),the Coulomb deflection(C),
                                                                           //and the relativity(R) effects

  //--------------------------------------------------------------------------------------------------------------------------------------------------

  if (crossSection >= 0) {
    return crossSection * barn;
  }
  else {return 0;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrBaseKxsModel::FunctionFK(G4double k, G4double theta)
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

  // PROTECTION TO ALLOW INTERPOLATION AT MINIMUM AND MAXIMUM EtaK/Theta2 values
  // (in particular for FK computation at 8.66EXX for high velocity formula)

  if (
  theta==8.66e-3 ||
  theta==8.66e-2 ||
  theta==8.66e-1 ||
  theta==8.66e+0 ||
  theta==8.66e+1
  ) theta=theta-1e-12;

  if (
  theta==1.e-3 ||
  theta==1.e-2 ||
  theta==1.e-1 ||
  theta==1.e+00 ||
  theta==1.e+01
  ) theta=theta+1e-12;

  // END PROTECTION

  auto t2 = std::upper_bound(dummyVec.begin(),dummyVec.end(), k);
  auto t1 = t2-1;

  auto e12 = std::upper_bound(aVecMap[(*t1)].begin(),aVecMap[(*t1)].end(), theta);
  auto e11 = e12-1;

  auto e22 = std::upper_bound(aVecMap[(*t2)].begin(),aVecMap[(*t2)].end(), theta);
  auto e21 = e22-1;

  valueT1  =*t1;
  valueT2  =*t2;
  valueE21 =*e21;
  valueE22 =*e22;
  valueE12 =*e12;
  valueE11 =*e11;

  xs11 = FKData[valueT1][valueE11];
  xs12 = FKData[valueT1][valueE12];
  xs21 = FKData[valueT2][valueE21];
  xs22 = FKData[valueT2][valueE22];

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

G4double G4ecpssrBaseKxsModel::LinLogInterpolate(G4double e1,
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

G4double G4ecpssrBaseKxsModel::LogLogInterpolate(G4double e1,
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

G4double G4ecpssrBaseKxsModel::QuadInterpolator(G4double e11, G4double e12,
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
