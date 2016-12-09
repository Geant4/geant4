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
// $Id: G4PenelopeBremsstrahlungFS.cc 101820 2016-12-01 08:34:05Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// --------
// 23 Nov 2010   L Pandola    First complete implementation
// 02 May 2011   L.Pandola    Remove dependency on CLHEP::HepMatrix
// 24 May 2011   L.Pandola    Renamed (make v2008 as default Penelope)
// 03 Oct 2013   L.Pandola    Migration to MT
// 30 Oct 2013   L.Pandola    Use G4Cache to avoid new/delete of the
//                             data vector on the fly in SampleGammaEnergy()
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include "G4PenelopeBremsstrahlungFS.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4AutoDelete.hh"
#include "G4Exp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungFS::G4PenelopeBremsstrahlungFS(G4int verbosity) :
  theReducedXSTable(0),theEffectiveZSq(0),theSamplingTable(0),
  thePBcut(0),fVerbosity(verbosity)
{
  fCache.Put(0);
  G4double tempvector[nBinsX] =
    {1.0e-12,0.025e0,0.05e0,0.075e0,0.1e0,0.15e0,0.2e0,0.25e0,
    0.3e0,0.35e0,0.4e0,0.45e0,0.5e0,0.55e0,0.6e0,0.65e0,0.7e0,
    0.75e0,0.8e0,0.85e0,0.9e0,0.925e0,0.95e0,0.97e0,0.99e0,
    0.995e0,0.999e0,0.9995e0,0.9999e0,0.99995e0,0.99999e0,1.0e0};

  for (size_t ix=0;ix<nBinsX;ix++)
    theXGrid[ix] = tempvector[ix];

  for (size_t i=0;i<nBinsE;i++)
    theEGrid[i] = 0.;

  theElementData = new std::map<G4int,G4DataVector*>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungFS::~G4PenelopeBremsstrahlungFS()
{
  ClearTables();

  //The G4Physics*Vector pointers contained in the fCache are automatically deleted by
  //the G4AutoDelete so there is no need to take care of them manually

  //Clear manually theElementData
  if (theElementData)
    {
      for (auto& item : (*theElementData)) 
	delete item.second;
      delete theElementData;
      theElementData = nullptr;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...


void G4PenelopeBremsstrahlungFS::ClearTables(G4bool isMaster)
{
  //Just to check
  if (!isMaster)
    G4Exception("G4PenelopeBremsstrahlungFS::ClearTables()",
		"em0100",FatalException,"Worker thread in this method");


  if (theReducedXSTable)
    {
      for (auto& item : (*theReducedXSTable))
	{
	  //G4PhysicsTable* tab = item.second;
	  //tab->clearAndDestroy();
	  delete item.second;
	}
      delete theReducedXSTable;
      theReducedXSTable = nullptr;
    }

  if (theSamplingTable)
    {
      for (auto& item : (*theSamplingTable))
	{
	  //G4PhysicsTable* tab = item.second;
	  // tab->clearAndDestroy();
          delete item.second;
	}
      delete theSamplingTable;
      theSamplingTable = nullptr;
    }
  if (thePBcut)
    {
      /*
	std::map< std::pair<const G4Material*,G4double> ,G4PhysicsFreeVector*>::iterator kk;
	for (kk=thePBcut->begin(); kk != thePBcut->end(); kk++)
	delete kk->second;
      */
      delete thePBcut;
      thePBcut = nullptr;
    }


  if (theEffectiveZSq)
    {
      delete theEffectiveZSq;
      theEffectiveZSq = nullptr;
    }

 return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungFS::GetEffectiveZSquared(const G4Material* material) const
{
  if (!theEffectiveZSq)
    {
      G4ExceptionDescription ed;
      ed << "The container for the <Z^2> values is not initialized" << G4endl;
      G4Exception("G4PenelopeBremsstrahlungFS::GetEffectiveZSquared()",
		  "em2007",FatalException,ed);
      return 0;
    }
  //found in the table: return it
  if (theEffectiveZSq->count(material))
    return theEffectiveZSq->find(material)->second;
  else
    {
      G4ExceptionDescription ed;
      ed << "The value of  <Z^2> is not properly set for material " <<
	material->GetName() << G4endl;
      //requires running of BuildScaledXSTable()
      G4Exception("G4PenelopeBremsstrahlungFS::GetEffectiveZSquared()",
		  "em2008",FatalException,ed);
    }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungFS::BuildScaledXSTable(const G4Material* material,
						    G4double cut,G4bool isMaster)
{
  //Corresponds to subroutines EBRaW and EBRaR of PENELOPE
  /*
    This method generates the table of the scaled energy-loss cross section from
    bremsstrahlung emission for the given material. Original data are read from
    file. The table is normalized according to the Berger-Seltzer cross section.
  */

  //Just to check
  if (!isMaster)
    G4Exception("G4PenelopeBremsstrahlungFS::BuildScaledXSTable()",
		"em0100",FatalException,"Worker thread in this method");

  if (fVerbosity > 2)
    {
      G4cout << "Entering in G4PenelopeBremsstrahlungFS::BuildScaledXSTable for " <<
	material->GetName() << G4endl;
      G4cout << "Threshold = " << cut/keV << " keV, isMaster= " << isMaster <<
	G4endl;
    }

  //This method should be accessed by the master only
  if (!theSamplingTable)
    theSamplingTable =
      new std::map< std::pair<const G4Material*,G4double> , G4PhysicsTable*>;
  if (!thePBcut)
    thePBcut =
      new std::map< std::pair<const G4Material*,G4double> , G4PhysicsFreeVector* >;

  //check if the container exists (if not, create it)
  if (!theReducedXSTable)
    theReducedXSTable = new std::map< std::pair<const G4Material*,G4double> ,
    G4PhysicsTable*>;
  if (!theEffectiveZSq)
    theEffectiveZSq = new std::map<const G4Material*,G4double>;



 //*********************************************************************
  //Determine the equivalent atomic number <Z^2>
  //*********************************************************************
  std::vector<G4double> *StechiometricFactors = new std::vector<G4double>;
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();
  const G4double* fractionVector = material->GetFractionVector();
  for (G4int i=0;i<nElements;i++)
    {
      G4double fraction = fractionVector[i];
      G4double atomicWeigth = (*elementVector)[i]->GetA()/(g/mole);
      StechiometricFactors->push_back(fraction/atomicWeigth);
    }
  //Find max
  G4double MaxStechiometricFactor = 0.;
  for (G4int i=0;i<nElements;i++)
    {
      if ((*StechiometricFactors)[i] > MaxStechiometricFactor)
        MaxStechiometricFactor = (*StechiometricFactors)[i];
    }
  //Normalize
  for (G4int i=0;i<nElements;i++)
    (*StechiometricFactors)[i] /=  MaxStechiometricFactor;

  G4double sumz2 = 0;
  G4double sums = 0;
  for (G4int i=0;i<nElements;i++)
    {
      G4double Z = (*elementVector)[i]->GetZ();
      sumz2 += (*StechiometricFactors)[i]*Z*Z;
      sums  += (*StechiometricFactors)[i];
    }
  G4double ZBR2 = sumz2/sums;

  theEffectiveZSq->insert(std::make_pair(material,ZBR2));


  //*********************************************************************
  // loop on elements and read data files
  //*********************************************************************
  G4DataVector* tempData = new G4DataVector(nBinsE);
  G4DataVector* tempMatrix = new G4DataVector(nBinsE*nBinsX,0.);

  for (G4int iel=0;iel<nElements;iel++)
    {
      G4double Z = (*elementVector)[iel]->GetZ();
      G4int iZ = (G4int) Z;
      G4double wgt = (*StechiometricFactors)[iel]*Z*Z/ZBR2;
      //

      //the element is not already loaded
      if (!theElementData->count(iZ))
	{
	  ReadDataFile(iZ);
	  if (!theElementData->count(iZ))
	    {
	      G4ExceptionDescription ed;
	      ed << "Error in G4PenelopeBremsstrahlungFS::BuildScaledXSTable" << G4endl;
	      ed << "Unable to retrieve data for element " << iZ << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungFS::BuildScaledXSTable()",
			  "em2009",FatalException,ed);
	    }
	}

      G4DataVector* atomData = theElementData->find(iZ)->second;

      for (size_t ie=0;ie<nBinsE;ie++)
	{
	  (*tempData)[ie] += wgt*(*atomData)[ie*(nBinsX+1)+nBinsX]; //last column contains total XS
	  for (size_t ix=0;ix<nBinsX;ix++)
	    (*tempMatrix)[ie*nBinsX+ix] += wgt*(*atomData)[ie*(nBinsX+1)+ix];
	}
    }

  //*********************************************************************
  // the total energy loss spectrum is re-normalized to reproduce the total
  // scaled cross section of Berger and Seltzer
  //*********************************************************************
  for (size_t ie=0;ie<nBinsE;ie++)
    {
      //for each energy, calculate integral of dSigma/dx over dx
      G4double* tempData2 = new G4double[nBinsX];
      for (size_t ix=0;ix<nBinsX;ix++)
	tempData2[ix] = (*tempMatrix)[ie*nBinsX+ix];
      G4double rsum = GetMomentumIntegral(tempData2,1.0,0);
      delete[] tempData2;
      G4double fact = millibarn*(theEGrid[ie]+electron_mass_c2)*(1./fine_structure_const)/
	(classic_electr_radius*classic_electr_radius*(theEGrid[ie]+2.0*electron_mass_c2));
      G4double fnorm = (*tempData)[ie]/(rsum*fact);
      G4double TST = 100.*std::fabs(fnorm-1.0);
      if (TST > 1.0)
	{
	  G4ExceptionDescription ed;
	  ed << "G4PenelopeBremsstrahlungFS. Corrupted data files?" << G4endl;
	  G4cout << "TST= " << TST << "; fnorm = " << fnorm << G4endl;
	  G4cout << "rsum = " << rsum << G4endl;
	  G4cout << "fact = " << fact << G4endl;
	  G4cout << ie << " " << theEGrid[ie]/keV << " " << (*tempData)[ie]/barn << G4endl;
	  G4Exception("G4PenelopeBremsstrahlungFS::BuildScaledXSTable()",
		      "em2010",FatalException,ed);
	}
      for (size_t ix=0;ix<nBinsX;ix++)
	(*tempMatrix)[ie*nBinsX+ix] *= fnorm;
    }

  //*********************************************************************
  // create and fill the tables
  //*********************************************************************
  G4PhysicsTable* thePhysicsTable = new G4PhysicsTable();
  // the table will contain 32 G4PhysicsFreeVectors with different
  // values of x. Each of the G4PhysicsFreeVectors has a profile of
  // log(XS) vs. log(E)

  //reserve space of the vectors. Everything is log-log
  //I add one extra "fake" point at low energy, since the Penelope
  //table starts at 1 keV
  for (size_t i=0;i<nBinsX;i++)
    thePhysicsTable->push_back(new G4PhysicsFreeVector(nBinsE+1));

  for (size_t ix=0;ix<nBinsX;ix++)
    {
      G4PhysicsFreeVector* theVec =
	(G4PhysicsFreeVector*) ((*thePhysicsTable)[ix]);
      for (size_t ie=0;ie<nBinsE;ie++)
	{
	  G4double logene = std::log(theEGrid[ie]);
	  G4double aValue = (*tempMatrix)[ie*nBinsX+ix];
	  if (aValue < 1e-20*millibarn) //protection against log(0)
	    aValue = 1e-20*millibarn;
	  theVec->PutValue(ie+1,logene,std::log(aValue));
	}
      //Add fake point at 1 eV using an extrapolation with the derivative
      //at the first valid point (Penelope approach)
      G4double derivative = ((*theVec)[2]-(*theVec)[1])/(theVec->Energy(2) - theVec->Energy(1));
      G4double log1eV = std::log(1*eV);
      G4double val1eV = (*theVec)[1]+derivative*(log1eV-theVec->Energy(1));
      //fake point at very low energy
      theVec->PutValue(0,log1eV,val1eV);
    }
  std::pair<const G4Material*,G4double> theKey = std::make_pair(material,cut);
  theReducedXSTable->insert(std::make_pair(theKey,thePhysicsTable));

  delete StechiometricFactors;
  delete tempData;
  delete tempMatrix;

  //Do here also the initialization of the energy sampling
  if (!(theSamplingTable->count(theKey)))
    InitializeEnergySampling(material,cut);

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungFS::ReadDataFile(G4int Z)
{

  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeBremsstrahlungFS - G4LEDATA environment variable not set!";
      G4Exception("G4PenelopeBremsstrahlungFS::ReadDataFile()",
		  "em0006",FatalException,excep);
      return;
    }
  /*
    Read the cross section file
  */
  std::ostringstream ost;
  if (Z>9)
    ost << path << "/penelope/bremsstrahlung/pdebr" << Z << ".p08";
  else
    ost << path << "/penelope/bremsstrahlung/pdebr0" << Z << ".p08";
  std::ifstream file(ost.str().c_str());
  if (!file.is_open())
    {
      G4String excep = "G4PenelopeBremsstrahlungFS - data file " +
	G4String(ost.str()) + " not found!";
      G4Exception("G4PenelopeBremsstrahlungFS::ReadDataFile()",
		  "em0003",FatalException,excep);
      return;
    }

  G4int readZ =0;
  file >> readZ;

  //check the right file is opened.
  if (readZ != Z)
    {
      G4ExceptionDescription ed;
      ed << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception("G4PenelopeBremsstrahlungFS::ReadDataFile()",
		  "em0005",FatalException,ed);
      return;
    }

  G4DataVector* theMatrix = new G4DataVector(nBinsE*(nBinsX+1),0.); //initialized with zeros

  for (size_t ie=0;ie<nBinsE;ie++)
    {
      G4double myDouble = 0;
      file >> myDouble; //energy (eV)
      if (!theEGrid[ie]) //fill only the first time
	theEGrid[ie] = myDouble*eV;
      //
      for (size_t ix=0;ix<nBinsX;ix++)
	{
	  file >> myDouble;
	  (*theMatrix)[ie*(nBinsX+1)+ix] = myDouble*millibarn;
	}
      file >> myDouble; //total cross section
      (*theMatrix)[ie*(nBinsX+1)+nBinsX] = myDouble*millibarn;
    }

  if (theElementData)
    theElementData->insert(std::make_pair(Z,theMatrix));
  else
    delete theMatrix;
  file.close();
  return;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungFS::GetMomentumIntegral(G4double* y,
							   G4double xup,G4int momOrder) const
//x is always the gridX
{
  //Corresponds to the function RLMOM of Penelope
  //This method performs the calculation of the integral of (x^momOrder)*y over the interval
  //from x[0] to xup, obtained by linear interpolation on a table of y.
  //The independent variable is assumed to take positive values only.
  //
  size_t size = nBinsX;
  const G4double eps = 1e-35;

  //Check that the call is valid
  if (momOrder<-1 || size<2 || theXGrid[0]<0)
    {
      G4Exception("G4PenelopeBremsstrahlungFS::GetMomentumIntegral()",
		  "em2011",FatalException,"Invalid call");
    }

  for (size_t i=1;i<size;i++)
    {
      if (theXGrid[i]<0 || theXGrid[i]<theXGrid[i-1])
	{
	  G4ExceptionDescription ed;
	  ed << "Invalid call for bin " << i << G4endl;
	  G4Exception("G4PenelopeBremsstrahlungFS::GetMomentumIntegral()",
		  "em2012",FatalException,ed);
	}
    }

  //Compute the integral
  G4double result = 0;
  if (xup < theXGrid[0])
    return result;
  bool loopAgain = true;
  G4double xt = std::min(xup,theXGrid[size-1]);
  G4double xtc = 0;
  for (size_t i=0;i<size-1;i++)
    {
      G4double x1 = std::max(theXGrid[i],eps);
      G4double y1 = y[i];
      G4double x2 = std::max(theXGrid[i+1],eps);
      G4double y2 = y[i+1];
      if (xt < x2)
	{
	  xtc = xt;
	  loopAgain = false;
	}
      else
	xtc = x2;
      G4double dx = x2-x1;
      G4double dy = y2-y1;
      G4double ds = 0;
      if (std::fabs(dx)>1e-14*std::fabs(dy))
	{
	  G4double b=dy/dx;
	  G4double a=y1-b*x1;
	  if (momOrder == -1)
	    ds = a*std::log(xtc/x1)+b*(xtc-x1);
	  else if (momOrder == 0) //speed it up, not using pow()
	    ds = a*(xtc-x1) + 0.5*b*(xtc*xtc-x1*x1);
	  else
	    ds = a*(std::pow(xtc,momOrder+1)-std::pow(x1,momOrder+1))/((G4double) (momOrder + 1))
	      + b*(std::pow(xtc,momOrder+2)-std::pow(x1,momOrder+2))/((G4double) (momOrder + 2));
	}
      else
	ds = 0.5*(y1+y2)*(xtc-x1)*std::pow(xtc,momOrder);
      result += ds;
      if (!loopAgain)
	return result;
    }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4PhysicsTable* G4PenelopeBremsstrahlungFS::GetScaledXSTable(const G4Material* mat,
								   const G4double cut) const
{
  //check if it already contains the entry
  std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);

  if (!(theReducedXSTable->count(theKey)))
    {
      G4Exception("G4PenelopeBremsstrahlungFS::GetScaledXSTable()",
		  "em2013",FatalException,"Unable to retrieve the cross section table");
    }

  return theReducedXSTable->find(theKey)->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungFS::InitializeEnergySampling(const G4Material* material,
							    G4double cut)
{
  if (fVerbosity > 2)
    G4cout << "Entering in G4PenelopeBremsstrahlungFS::InitializeEnergySampling() for " <<
      material->GetName() << G4endl;

  //This method should be accessed by the master only
  std::pair<const G4Material*,G4double> theKey = std::make_pair(material,cut);

  G4PhysicsTable* thePhysicsTable = new G4PhysicsTable();
  // the table will contain 57 G4PhysicsFreeVectors with different
  // values of E.
  G4PhysicsFreeVector* thePBvec = new G4PhysicsFreeVector(nBinsE);

  //I reserve space of the vectors.
  for (size_t i=0;i<nBinsE;i++)
    thePhysicsTable->push_back(new G4PhysicsFreeVector(nBinsX));

  //Retrieve the table. Must already exist at this point, because this
  //method is invoked by GetScaledXSTable()
  if (!(theReducedXSTable->count(theKey)))
    G4Exception("G4PenelopeBremsstrahlungFS::InitializeEnergySampling()",
		"em2013",FatalException,"Unable to retrieve the cross section table");
  G4PhysicsTable* theTableReduced = theReducedXSTable->find(theKey)->second;

  for (size_t ie=0;ie<nBinsE;ie++)
    {
      G4PhysicsFreeVector* theVec =
	(G4PhysicsFreeVector*) ((*thePhysicsTable)[ie]);
      //Fill the table
      G4double value = 0; //first value
      theVec->PutValue(0,theXGrid[0],value);
      for (size_t ix=1;ix<nBinsX;ix++)
	{
	  //Here calculate the cumulative distribution
	  // int_{0}^{x} dSigma(x',E)/dx' (1/x') dx'
	  G4PhysicsFreeVector* v1 = (G4PhysicsFreeVector*) (*theTableReduced)[ix-1];
	  G4PhysicsFreeVector* v2 = (G4PhysicsFreeVector*) (*theTableReduced)[ix];

	  G4double x1=std::max(theXGrid[ix-1],1.0e-35);
	  //Remember: the table theReducedXSTable has a fake first point in energy
	  //so, it contains one more bin than nBinsE.
	  G4double y1=G4Exp((*v1)[ie+1]);
	  G4double x2=std::max(theXGrid[ix],1.0e-35);
	  G4double y2=G4Exp((*v2)[ie+1]);
	  G4double B = (y2-y1)/(x2-x1);
	  G4double A = y1-B*x1;
	  G4double dS = A*std::log(x2/x1)+B*(x2-x1);
	  value += dS;
	  theVec->PutValue(ix,theXGrid[ix],value);
	}

      //fill the PB vector
      G4double xc = cut/theEGrid[ie];
      //Fill a temp data vector
      G4double* tempData = new G4double[nBinsX];
      for (size_t ix=0;ix<nBinsX;ix++)
	{
	  G4PhysicsFreeVector* vv = (G4PhysicsFreeVector*) (*theTableReduced)[ix];
	  tempData[ix] = G4Exp((*vv)[ie+1]);
	}
      G4double pbval = (xc<=1) ?
	GetMomentumIntegral(tempData,xc,-1) :
	GetMomentumIntegral(tempData,1.0,-1);
      thePBvec->PutValue(ie,theEGrid[ie],pbval);
      delete[] tempData;
    }

  theSamplingTable->insert(std::make_pair(theKey,thePhysicsTable));
  thePBcut->insert(std::make_pair(theKey,thePBvec));
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungFS::SampleGammaEnergy(G4double energy,const G4Material* mat,
							     const G4double cut) const
{
  std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
  if (!(theSamplingTable->count(theKey)) || !(thePBcut->count(theKey)))
    {
      G4ExceptionDescription ed;
      ed << "Unable to retrieve the SamplingTable: " <<
	theSamplingTable->count(theKey) << " " <<
	thePBcut->count(theKey) << G4endl;
      G4Exception("G4PenelopeBremsstrahlungFS::SampleGammaEnergy()",
		  "em2014",FatalException,ed);
    }


  const G4PhysicsTable* theTableInte = theSamplingTable->find(theKey)->second;
  const G4PhysicsTable* theTableRed = theReducedXSTable->find(theKey)->second;

  //Find the energy bin using bi-partition
  size_t eBin = 0;
  G4bool firstOrLastBin = false;

  if (energy < theEGrid[0]) //below first bin
    {
      eBin = 0;
      firstOrLastBin = true;
    }
  else if (energy > theEGrid[nBinsE-1]) //after last bin
    {
      eBin = nBinsE-1;
      firstOrLastBin = true;
    }
  else
    {
      size_t i=0;
      size_t j=nBinsE-1;
      while ((j-i)>1)
	{
	  size_t k = (i+j)/2;
	  if (energy > theEGrid[k])
	    i = k;
	  else
	    j = k;
	}
      eBin = i;
    }

  //Get the appropriate physics vector
  const G4PhysicsFreeVector* theVec1 = (G4PhysicsFreeVector*) (*theTableInte)[eBin];

  //Use a "temporary" vector which contains the linear interpolation of the x spectra
  //in energy. The temporary vector is thread-local, so that there is no conflict.
  //This is achieved via G4Cache. The theTempVect is allocated only once per thread
  //(member variable), but it is overwritten at every call of this method
  //(because the interpolation factors change!)
  G4PhysicsFreeVector* theTempVec = fCache.Get();
  if (!theTempVec) //First time this thread gets the cache
    {
      theTempVec = new G4PhysicsFreeVector(nBinsX);   
      fCache.Put(theTempVec);
      // The G4AutoDelete takes care here to clean up the vectors
      G4AutoDelete::Register(theTempVec);
      if (fVerbosity > 4)
	G4cout << "Creating new instance of G4PhysicsFreeVector() on the worker" << G4endl;
    }

  //theTempVect is allocated only once (member variable), but it is overwritten at
  //every call of this method (because the interpolation factors change!)
  if (!firstOrLastBin)
    {
      const G4PhysicsFreeVector* theVec2 = (G4PhysicsFreeVector*) (*theTableInte)[eBin+1];
      for (size_t iloop=0;iloop<nBinsX;iloop++)
	{
	  G4double val = (*theVec1)[iloop]+(((*theVec2)[iloop]-(*theVec1)[iloop]))*
	    (energy-theEGrid[eBin])/(theEGrid[eBin+1]-theEGrid[eBin]);
	  theTempVec->PutValue(iloop,theXGrid[iloop],val);
	}
    }
  else //first or last bin, no interpolation
    {
      for (size_t iloop=0;iloop<nBinsX;iloop++)
	theTempVec->PutValue(iloop,theXGrid[iloop],(*theVec1)[iloop]);
    }

  //Start the game
  G4double pbcut = (*(thePBcut->find(theKey)->second))[eBin];

  if (!firstOrLastBin) //linear interpolation on pbcut as well
    {
      pbcut = (*(thePBcut->find(theKey)->second))[eBin] +
	((*(thePBcut->find(theKey)->second))[eBin+1]-(*(thePBcut->find(theKey)->second))[eBin])*
	(energy-theEGrid[eBin])/(theEGrid[eBin+1]-theEGrid[eBin]);
    }

  G4double pCumulative = (*theTempVec)[nBinsX-1]; //last value

  G4double eGamma = 0;
  G4int nIterations = 0;
  do
    {
      G4double pt = pbcut + G4UniformRand()*(pCumulative - pbcut);
      nIterations++;

      //find where it is
      size_t ibin = 0;
      if (pt < (*theTempVec)[0])
	ibin = 0;
      else if (pt > (*theTempVec)[nBinsX-1])
	{
	  //We observed problems due to numerical rounding here (STT).
	  //delta here is a tiny positive number
	  G4double delta = pt-(*theTempVec)[nBinsX-1];
	  if (delta < pt*1e-10) // very small! Numerical rounding only
	    {
	      ibin = nBinsX-2;
	      G4ExceptionDescription ed;
	      ed << "Found that (pt > (*theTempVec)[nBinsX-1]) with pt = " << pt <<
		" , (*theTempVec)[nBinsX-1] = " << (*theTempVec)[nBinsX-1] << " and delta = " <<
		(pt-(*theTempVec)[nBinsX-1]) << G4endl;
	      ed << "Possible symptom of problem with numerical precision" << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungFS::SampleGammaEnergy()",
			  "em2015",JustWarning,ed);
	    }
	  else //real problem
	    {
	      G4ExceptionDescription ed;
	      ed << "Crash at (pt > (*theTempVec)[nBinsX-1]) with pt = " << pt <<
		" , (*theTempVec)[nBinsX-1]=" << (*theTempVec)[nBinsX-1] << " and nBinsX = " <<
		nBinsX << G4endl;
	      ed << "Material: " << mat->GetName() << ", energy = " << energy/keV << " keV" <<
		G4endl;
	      G4Exception("G4PenelopeBremsstrahlungFS::SampleGammaEnergy()",
			  "em2015",FatalException,ed);
	    }
	}
      else
	{
	  size_t i=0;
	  size_t j=nBinsX-1;
	  while ((j-i)>1)
	    {
	      size_t k = (i+j)/2;
	      if (pt > (*theTempVec)[k])
		i = k;
	      else
		j = k;
	    }
	  ibin = i;
	}

      G4double w1 = theXGrid[ibin];
      G4double w2 = theXGrid[ibin+1];

      const G4PhysicsFreeVector* v1 = (G4PhysicsFreeVector*) (*theTableRed)[ibin];
      const G4PhysicsFreeVector* v2 = (G4PhysicsFreeVector*) (*theTableRed)[ibin+1];
      //Remember: the table theReducedXSTable has a fake first point in energy
      //so, it contains one more bin than nBinsE.
      G4double pdf1 = G4Exp((*v1)[eBin+1]);
      G4double pdf2 = G4Exp((*v2)[eBin+1]);
      G4double deltaW = w2-w1;
      G4double dpdfb = pdf2-pdf1;
      G4double B = dpdfb/deltaW;
      G4double A = pdf1-B*w1;
      //I already made an interpolation in energy, so I can use the actual value for the
      //calculation of the wbcut, instead of the grid values (except for the last bin)
      G4double wbcut  = (cut < energy) ? cut/energy : 1.0;
      if (firstOrLastBin) //this is an particular case: no interpolation available
	wbcut  = (cut < theEGrid[eBin]) ? cut/theEGrid[eBin] : 1.0;

      if (w1 < wbcut)
	w1 = wbcut;
      if (w2 < w1)
	{
	  //This configuration can happen if initially wbcut > w2 > w1. Due to the previous
	  //statement, (w1 = wbcut), it becomes wbcut = w1 > w2. In this case, it is not a
	  //real problem. It becomes a problem if w2 < w1 before the w1 = wbcut statement. Issue
	  //a warning only in this specific case.
	  if (w2 > wbcut)
	    {
	      G4ExceptionDescription ed;
	      ed << "Warning in G4PenelopeBremsstrahlungFS::SampleX()" << G4endl;
	      ed << "Conflicting end-point values: w1=" << w1 << "; w2 = " << w2 << G4endl;
	      ed << "wbcut = " << wbcut << " energy= " << energy/keV << " keV" << G4endl;
	      ed << "cut = " << cut/keV << " keV" << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungFS::SampleGammaEnergy()","em2015",
			  JustWarning,ed);
	    }
	  return w1*energy;
	}

      G4double pmax = std::max(A+B*w1,A+B*w2);
      G4bool loopAgain = false;
      do
	{
	  loopAgain = false;
	  eGamma = w1* std::pow((w2/w1),G4UniformRand());
	  if  (G4UniformRand()*pmax > (A+B*eGamma))
	    loopAgain = true;
	}while(loopAgain);
      eGamma *= energy;
      if (nIterations > 100) //protection against infinite loops
	return eGamma;
    }while(eGamma < cut); //repeat if sampled sub-cut!

  return eGamma;
}
