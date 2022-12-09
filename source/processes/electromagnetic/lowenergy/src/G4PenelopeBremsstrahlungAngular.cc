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
// --------------------------------------------------------------
//
// File name:     G4PenelopeBremsstrahlungAngular
//
// Author:        Luciano Pandola
//
// Creation date: November 2010
//
// History:
// -----------
// 23 Nov 2010  L. Pandola       1st implementation
// 24 May 2011  L. Pandola       Renamed (make v2008 as default Penelope)
// 13 Mar 2012  L. Pandola       Made a derived class of G4VEmAngularDistribution
//                               and update the interface accordingly
// 18 Jul 2012  L. Pandola       Migrated to the new basic interface of G4VEmAngularDistribution
//                               Now returns a G4ThreeVector and takes care of the rotation
// 03 Oct 2013  L. Pandola       Migrated to MT: only the master model handles tables
// 17 Oct 2013  L. Pandola       Partially revert MT migration. The angular generator is kept as
//                                thread-local, and each worker has full access to it.
//
//----------------------------------------------------------------

#include "G4PenelopeBremsstrahlungAngular.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4Exp.hh"

G4PenelopeBremsstrahlungAngular::G4PenelopeBremsstrahlungAngular() :
  G4VEmAngularDistribution("Penelope"), fEffectiveZSq(nullptr),
  fLorentzTables1(nullptr),fLorentzTables2(nullptr)

{
  fDataRead = false;
  fVerbosityLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungAngular::~G4PenelopeBremsstrahlungAngular()
{
  ClearTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungAngular::Initialize()
{
  ClearTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungAngular::ClearTables()
{
  if (fLorentzTables1)
    {
      for (auto j=fLorentzTables1->begin(); j != fLorentzTables1->end(); ++j)
        {
	  G4PhysicsTable* tab = j->second;
          tab->clearAndDestroy();
          delete tab;
        }
      fLorentzTables1->clear();
      delete fLorentzTables1;
      fLorentzTables1 = nullptr;
    }

  if (fLorentzTables2)
    {
      for (auto j=fLorentzTables2->begin(); j != fLorentzTables2->end(); ++j)
        {
	  G4PhysicsTable* tab = j->second;
	  tab->clearAndDestroy();
          delete tab;
        }
      fLorentzTables2->clear();
      delete fLorentzTables2;
      fLorentzTables2 = nullptr;
    }
  if (fEffectiveZSq)
    {
      delete fEffectiveZSq;
      fEffectiveZSq = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungAngular::ReadDataFile()
{
   //Read information from DataBase file
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    {
      G4String excep =
	"G4PenelopeBremsstrahlungAngular - G4LEDATA environment variable not set!";
      G4Exception("G4PenelopeBremsstrahlungAngular::ReadDataFile()",
		  "em0006",FatalException,excep);
      return;
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/bremsstrahlung/pdbrang.p08";
  std::ifstream file(pathFile);

  if (!file.is_open())
    {
      G4String excep = "G4PenelopeBremsstrahlungAngular - data file " + pathFile + " not found!";
      G4Exception("G4PenelopeBremsstrahlungAngular::ReadDataFile()",
		  "em0003",FatalException,excep);
      return;
    }
  G4int i=0,j=0,k=0; // i=index for Z, j=index for E, k=index for K

  for (k=0;k<fNumberofKPoints;k++)
    for (i=0;i<fNumberofZPoints;i++)
      for (j=0;j<fNumberofEPoints;j++)
	{
	  G4double a1,a2;
	  G4int ik1,iz1,ie1;
	  G4double zr,er,kr;
	  file >> iz1 >> ie1 >> ik1 >> zr >> er >> kr >> a1 >> a2;
	  //check the data are correct
	  if ((iz1-1 == i) && (ik1-1 == k) && (ie1-1 == j))
	    {
	      fQQ1[i][j][k]=a1;
	      fQQ2[i][j][k]=a2;
	    }
	  else
	    {
	      G4ExceptionDescription ed;
	      ed << "Corrupted data file " << pathFile << "?" << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungAngular::ReadDataFile()",
		  "em0005",FatalException,ed);
	    }
	}
  file.close();
  fDataRead = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungAngular::PrepareTables(const G4Material* material,G4bool /*isMaster*/ )
{
  //Unused at the moment: the G4PenelopeBremsstrahlungAngular is thread-local, so each worker
  //builds its own version of the tables.
  
  //Check if data file has already been read
  if (!fDataRead)
    {
      ReadDataFile();
      if (!fDataRead)
	G4Exception("G4PenelopeBremsstrahlungAngular::PrepareInterpolationTables()",
		    "em2001",FatalException,"Unable to build interpolation table");
    }

  if (!fLorentzTables1)
      fLorentzTables1 = new std::map<G4double,G4PhysicsTable*>;
  if (!fLorentzTables2)
    fLorentzTables2 = new std::map<G4double,G4PhysicsTable*>;

  G4double Zmat = CalculateEffectiveZ(material);

  const G4int reducedEnergyGrid=21;
  //Support arrays.
  G4double betas[fNumberofEPoints]; //betas for interpolation
  //tables for interpolation
  G4double Q1[fNumberofEPoints][fNumberofKPoints];
  G4double Q2[fNumberofEPoints][fNumberofKPoints];
  //expanded tables for interpolation
  G4double Q1E[fNumberofEPoints][reducedEnergyGrid];
  G4double Q2E[fNumberofEPoints][reducedEnergyGrid];
  G4double pZ[fNumberofZPoints] = {2.0,8.0,13.0,47.0,79.0,92.0};

  G4int i=0,j=0,k=0; // i=index for Z, j=index for E, k=index for K
  //Interpolation in Z
  for (i=0;i<fNumberofEPoints;i++)
    {
      for (j=0;j<fNumberofKPoints;j++)
	{
	  G4PhysicsFreeVector* QQ1vector = 
	    new G4PhysicsFreeVector(fNumberofZPoints,/*spline =*/true);
	  G4PhysicsFreeVector* QQ2vector = 
	    new G4PhysicsFreeVector(fNumberofZPoints,/*spline =*/true);

	  //fill vectors
	  for (k=0;k<fNumberofZPoints;k++)
	    {
	      QQ1vector->PutValues(k,pZ[k],G4Log(fQQ1[k][i][j]));
	      QQ2vector->PutValues(k,pZ[k],fQQ2[k][i][j]);
	    }
	  //Filled: now calculate derivatives
	  QQ1vector->FillSecondDerivatives();
	  QQ2vector->FillSecondDerivatives();

	  Q1[i][j]= G4Exp(QQ1vector->Value(Zmat));
	  Q2[i][j]=QQ2vector->Value(Zmat);
	  delete QQ1vector;
	  delete QQ2vector;
	}
    }
  G4double pE[fNumberofEPoints] = {1.0e-03*MeV,5.0e-03*MeV,1.0e-02*MeV,5.0e-02*MeV,
				  1.0e-01*MeV,5.0e-01*MeV};
  G4double pK[fNumberofKPoints] = {0.0,0.6,0.8,0.95};
  G4double ppK[reducedEnergyGrid];

  for(i=0;i<reducedEnergyGrid;i++)
    ppK[i]=((G4double) i) * 0.05;


  for(i=0;i<fNumberofEPoints;i++)
    betas[i]=std::sqrt(pE[i]*(pE[i]+2*electron_mass_c2))/(pE[i]+electron_mass_c2);


  for (i=0;i<fNumberofEPoints;i++)
    {
      for (j=0;j<fNumberofKPoints;j++)
	Q1[i][j]=Q1[i][j]/Zmat;
    }

  //Expanded table of distribution parameters
  for (i=0;i<fNumberofEPoints;i++)
    {
      G4PhysicsFreeVector* Q1vector = new G4PhysicsFreeVector(fNumberofKPoints);
      G4PhysicsFreeVector* Q2vector = new G4PhysicsFreeVector(fNumberofKPoints);

      for (j=0;j<fNumberofKPoints;j++)
	{
	  Q1vector->PutValues(j,pK[j],G4Log(Q1[i][j])); //logarithmic
	  Q2vector->PutValues(j,pK[j],Q2[i][j]);
	}

      for (j=0;j<reducedEnergyGrid;j++)
	{
	  Q1E[i][j]=Q1vector->Value(ppK[j]);
	  Q2E[i][j]=Q2vector->Value(ppK[j]);
	}
      delete Q1vector;
      delete Q2vector;
    }
  //
  //TABLES to be stored
  //
  G4PhysicsTable* theTable1 = new G4PhysicsTable();
  G4PhysicsTable* theTable2 = new G4PhysicsTable();
  // the table will contain reducedEnergyGrid G4PhysicsFreeVectors with different
  // values of k,
  // Each of the G4PhysicsFreeVectors has a profile of
  // y vs. E
  //
  //reserve space of the vectors.
  for (j=0;j<reducedEnergyGrid;j++)
    {
      G4PhysicsFreeVector* thevec = new G4PhysicsFreeVector(fNumberofEPoints,/*spline=*/true);
      theTable1->push_back(thevec);
      G4PhysicsFreeVector* thevec2 = new G4PhysicsFreeVector(fNumberofEPoints,/*spline=*/true);
      theTable2->push_back(thevec2);
    }

  for (j=0;j<reducedEnergyGrid;j++)
    {
      G4PhysicsFreeVector* thevec = (G4PhysicsFreeVector*) (*theTable1)[j];
      G4PhysicsFreeVector* thevec2 = (G4PhysicsFreeVector*) (*theTable2)[j];
      for (i=0;i<fNumberofEPoints;i++)
	{
	  thevec->PutValues(i,betas[i],Q1E[i][j]);
	  thevec2->PutValues(i,betas[i],Q2E[i][j]);
	}
      //Vectors are filled: calculate derivatives
      thevec->FillSecondDerivatives();
      thevec2->FillSecondDerivatives();
    }

  if (fLorentzTables1 && fLorentzTables2)
    {
      fLorentzTables1->insert(std::make_pair(Zmat,theTable1));
      fLorentzTables2->insert(std::make_pair(Zmat,theTable2));
    }
  else
    {
      G4ExceptionDescription ed;
      ed << "Unable to create tables of Lorentz coefficients for " << G4endl;
      ed << "<Z>= "  << Zmat << " in G4PenelopeBremsstrahlungAngular" << G4endl;
      delete theTable1;
      delete theTable2;
      G4Exception("G4PenelopeBremsstrahlungAngular::PrepareInterpolationTables()",
		  "em2005",FatalException,ed);
    }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector& G4PenelopeBremsstrahlungAngular::SampleDirection(const G4DynamicParticle* dp,
								G4double eGamma,
								G4int,
								const G4Material* material)
{
  if (!material)
    {
      G4Exception("G4PenelopeBremsstrahlungAngular::SampleDirection()",
		  "em2040",FatalException,"The pointer to G4Material* is nullptr");
      return fLocalDirection;
    }

  //Retrieve the effective Z
  G4double Zmat = 0;

  if (!fEffectiveZSq)
    {
      G4Exception("G4PenelopeBremsstrahlungAngular::SampleDirection()",
		  "em2040",FatalException,"EffectiveZ table not available");
      return fLocalDirection;
    }

  //found in the table: return it
  if (fEffectiveZSq->count(material))
    Zmat = fEffectiveZSq->find(material)->second;
  else
    {
      G4Exception("G4PenelopeBremsstrahlungAngular::SampleDirection()",
		  "em2040",FatalException,"Material not found in the effectiveZ table");
      return fLocalDirection;
    }

  if (fVerbosityLevel > 0)
    {
      G4cout << "Effective <Z> for material : " << material->GetName() <<
	" = " << Zmat << G4endl;
    }

  G4double ePrimary = dp->GetKineticEnergy();

  G4double beta = std::sqrt(ePrimary*(ePrimary+2*electron_mass_c2))/
    (ePrimary+electron_mass_c2);
  G4double cdt = 0;
  G4double sinTheta = 0;
  G4double phi = 0;

  //Use a pure dipole distribution for energy above 500 keV
  if (ePrimary > 500*keV)
    {
      cdt = 2.0*G4UniformRand() - 1.0;
      if (G4UniformRand() > 0.75)
	{
	  if (cdt<0)
	    cdt = -1.0*std::pow(-cdt,1./3.);
	  else
	    cdt = std::pow(cdt,1./3.);
	}
      cdt = (cdt+beta)/(1.0+beta*cdt);
      //Get primary kinematics
      sinTheta = std::sqrt(1. - cdt*cdt);
      phi  = twopi * G4UniformRand();
      fLocalDirection.set(sinTheta* std::cos(phi),
		      sinTheta* std::sin(phi),
		      cdt);
      //rotate
      fLocalDirection.rotateUz(dp->GetMomentumDirection());
      //return
      return fLocalDirection;
    }

  if (!(fLorentzTables1->count(Zmat)) || !(fLorentzTables2->count(Zmat)))
    {
      G4ExceptionDescription ed;
      ed << "Unable to retrieve Lorentz tables for Z= " << Zmat << G4endl;
      G4Exception("G4PenelopeBremsstrahlungAngular::SampleDirection()",
		  "em2006",FatalException,ed);
    }

  //retrieve actual tables
  const G4PhysicsTable* theTable1 = fLorentzTables1->find(Zmat)->second;
  const G4PhysicsTable* theTable2 = fLorentzTables2->find(Zmat)->second;

  G4double RK=20.0*eGamma/ePrimary;
  G4int ik=std::min((G4int) RK,19);

  G4double P10=0,P11=0,P1=0;
  G4double P20=0,P21=0,P2=0;

  //First coefficient
  const G4PhysicsFreeVector* v1 = (G4PhysicsFreeVector*) (*theTable1)[ik];
  const G4PhysicsFreeVector* v2 = (G4PhysicsFreeVector*) (*theTable1)[ik+1];
  P10 = v1->Value(beta);
  P11 = v2->Value(beta);
  P1=P10+(RK-(G4double) ik)*(P11-P10);

  //Second coefficient
  const G4PhysicsFreeVector* v3 = (G4PhysicsFreeVector*) (*theTable2)[ik];
  const G4PhysicsFreeVector* v4 = (G4PhysicsFreeVector*) (*theTable2)[ik+1];
  P20=v3->Value(beta);
  P21=v4->Value(beta);
  P2=P20+(RK-(G4double) ik)*(P21-P20);

  //Sampling from the Lorenz-trasformed dipole distributions
  P1=std::min(G4Exp(P1)/beta,1.0);
  G4double betap = std::min(std::max(beta*(1.0+P2/beta),0.0),0.9999);

  G4double testf=0;

  if (G4UniformRand() < P1)
    {
      do{
	cdt = 2.0*G4UniformRand()-1.0;
	testf=2.0*G4UniformRand()-(1.0+cdt*cdt);
      }while(testf>0);
    }
  else
    {
      do{
	cdt = 2.0*G4UniformRand()-1.0;
	testf=G4UniformRand()-(1.0-cdt*cdt);
      }while(testf>0);
    }
  cdt = (cdt+betap)/(1.0+betap*cdt);

  //Get primary kinematics
  sinTheta = std::sqrt(1. - cdt*cdt);
  phi  = twopi * G4UniformRand();
  fLocalDirection.set(sinTheta* std::cos(phi),
		      sinTheta* std::sin(phi),
		      cdt);
  //rotate
  fLocalDirection.rotateUz(dp->GetMomentumDirection());
  //return
  return fLocalDirection;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungAngular::CalculateEffectiveZ(const G4Material* material)
{
  if (!fEffectiveZSq)
    fEffectiveZSq = new std::map<const G4Material*,G4double>;

  //Already exists: return it
  if (fEffectiveZSq->count(material))
    return fEffectiveZSq->find(material)->second;

  //Helper for the calculation
  std::vector<G4double> *StechiometricFactors = new std::vector<G4double>;
  G4int nElements = (G4int)material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();
  const G4double* fractionVector = material->GetFractionVector();
  for (G4int i=0;i<nElements;++i)
    {
      G4double fraction = fractionVector[i];
      G4double atomicWeigth = (*elementVector)[i]->GetA()/(g/mole);
      StechiometricFactors->push_back(fraction/atomicWeigth);
    }
  //Find max
  G4double MaxStechiometricFactor = 0.;
  for (G4int i=0;i<nElements;++i)
    {
      if ((*StechiometricFactors)[i] > MaxStechiometricFactor)
        MaxStechiometricFactor = (*StechiometricFactors)[i];
    }
  //Normalize
  for (G4int i=0;i<nElements;++i)
    (*StechiometricFactors)[i] /=  MaxStechiometricFactor;

  G4double sumz2 = 0;
  G4double sums = 0;
  for (G4int i=0;i<nElements;++i)
    {
      G4double Z = (*elementVector)[i]->GetZ();
      sumz2 += (*StechiometricFactors)[i]*Z*Z;
      sums  += (*StechiometricFactors)[i];
    }
  delete StechiometricFactors;

  G4double ZBR = std::sqrt(sumz2/sums);
  fEffectiveZSq->insert(std::make_pair(material,ZBR));

  return ZBR;
}
