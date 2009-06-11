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
// --------------------------------------------------------------------
//
// $Id: G4PenelopeRayleigh.cc,v 1.19 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L. Pandola (luciano.pandola@cern.ch)
//
// History:
// --------
// 14 Feb 2003   MG Pia       Corrected compilation errors and warnings
//                            from SUN
// 10 Mar 2003 V.Ivanchenko   Remove CutPerMaterial warning
// 12 Mar 2003 L.Pandola      Code "cleaned" - Cuts per region 
// 17 Mar 2004 L.Pandola      Removed unnecessary calls to std::pow(a,b)
// 18 Mar 2004 M.Mendenhall   Introduced SamplingTable (performance improvement)
// 03 Sep 2007 L.Pandola      Bug fix for the filling of physics table for 
//                            compounds defined by the mass fraction (bug #965)
// 02 Apr 2009 L.Pandola      Bux fixed in the calculation of mfp for compound 
//                            materials defined as fraction of mass 
//                            (reported by Zhang Qiwei)
// --------------------------------------------------------------------

#include "G4PenelopeRayleigh.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4PenelopeIntegrator.hh"
#include "G4MaterialCutsCouple.hh"

G4PenelopeRayleigh::G4PenelopeRayleigh(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),             
    highEnergyLimit(100*GeV),
    samplingConstant(0.0),
    nBins(200),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV)

{
  if (lowEnergyLimit < intrinsicLowEnergyLimit ||
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4PenelopeRayleigh::G4PenelopeRayleigh - energy limit outside intrinsic process validity range");
    }

  samplingFunction_x = 0;
  samplingFunction_y = 0;
  meanFreePathTable = 0;
  material = 0;

   if (verboseLevel > 0)
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: "
	      << lowEnergyLimit / keV << " keV - "
	      << highEnergyLimit / GeV << " GeV"
	      << G4endl;
     }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4PenelopeRayleigh is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

G4PenelopeRayleigh::~G4PenelopeRayleigh()
{
	delete meanFreePathTable;
	
	SamplingTablePair::iterator i;
	for(i=SamplingTables.begin(); i != SamplingTables.end(); i++) {
		delete (*i).second.first;
		delete (*i).second.second;
	}
}


void G4PenelopeRayleigh::BuildPhysicsTable(const G4ParticleDefinition& )
{

  G4DataVector energyVector;
  G4double dBin = std::log10(highEnergyLimit/lowEnergyLimit)/nBins;
  for (G4int i=0;i<nBins;i++)
    {
      energyVector.push_back(std::pow(10.,std::log10(lowEnergyLimit)+i*dBin));
    }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();

  size_t nOfBins = energyVector.size();
  size_t bin=0;

  G4VDataSetAlgorithm* algo = new G4LogLogInterpolation();
  G4VEMDataSet* materialSet = new G4CompositeEMDataSet(algo,1.,1.);
  std::vector<G4VEMDataSet*> matCrossSections;

  G4int m;
  for (m=0; m<nMaterials; m++)
    {
      G4DataVector* energies = new G4DataVector;
      G4DataVector* data = new G4DataVector;
      material = (*materialTable)[m];

      G4int nElements = material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double k1=849.3315;


      G4int IZZ=0;
      G4int iright=0;
      for (G4int i=0; i<nElements; i++) {
        G4int Z = (G4int) (*elementVector)[i]->GetZ();
	if (Z>IZZ){
	  IZZ = Z;
	  iright=i;
	}
      }
      for (bin=0; bin<nOfBins; bin++)
	{
	  energies->push_back(energyVector[bin]);
	  G4double ec=std::min(energyVector[bin],0.5*IZZ);
	  G4double energyRatio = ec/electron_mass_c2;
	  facte = k1*energyRatio*energyRatio;
	  G4double cs=0;
	  G4PenelopeIntegrator<G4PenelopeRayleigh,G4double(G4PenelopeRayleigh::*)(G4double)> theIntegrator;
	  cs =
	    theIntegrator.Calculate(this,&G4PenelopeRayleigh::DifferentialCrossSection,-1.0,0.90,1e-06);
	  cs += theIntegrator.Calculate(this,&G4PenelopeRayleigh::DifferentialCrossSection,0.90,0.9999999,1e-06);
	  cs = cs*(ec/energyVector[bin])*(ec/energyVector[bin])*pi*classic_electr_radius*classic_electr_radius;
	  const G4double* vector_of_atoms = material->GetVecNbOfAtomsPerVolume();
	  const G4int* stechiometric = material->GetAtomsVector();
	  //cs is the cross section _per atom_ in the case of compounds, while it is 
	  //_per molecule_ in the case of molecules
	  G4double density = material->GetTotNbOfAtomsPerVolume();  //non-bound molecules (default)
	  if (stechiometric)
	    {
	      if (stechiometric[iright])
		density = vector_of_atoms[iright]/stechiometric[iright]; //number of molecules per volume
	    }
	  G4double cross = density*cs;
	  data->push_back(cross);
	}
      G4VEMDataSet* elSet = new G4EMDataSet(0,energies,data,algo); 
      G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo);
      setForMat->AddComponent(elSet); 
      matCrossSections.push_back(setForMat);
    }
  G4double matCS = 0.0;

  for (m=0; m<nMaterials; m++)
    {
      G4DataVector* energies = new G4DataVector;
      G4DataVector* data = new G4DataVector;
      for (bin=0;bin<nOfBins;bin++){
	energies->push_back(energyVector[bin]);
	matCS = (matCrossSections[m]->GetComponent(0))->FindValue(energyVector[bin]);
	if (matCS > 0.){
	  data->push_back(1./matCS);
	}
	else
	  {
	    data->push_back(DBL_MAX);
	  }
      }
      G4VEMDataSet* dataSet = new G4EMDataSet(m,energies,data,algo,1.,1.); 
      materialSet->AddComponent(dataSet);
    }
  meanFreePathTable = materialSet;
}

G4VParticleChange* G4PenelopeRayleigh::PostStepDoIt(const G4Track& aTrack,
						     const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();

  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }


  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  material = couple->GetMaterial();
  
  // Sampling inizialitation (build internal table)
 
  InizialiseSampling();

  // Sample the angle of the scattered photon
  const G4double xpar=41.2148;
  G4double x2max = 2.0*std::log(xpar*photonEnergy0/electron_mass_c2);
  G4int jm;
  G4int asize = samplingFunction_x->size();
  if (x2max<(*samplingFunction_x)[1]) 
    {
      jm=0;
    }
  else if(x2max>(*samplingFunction_x)[asize-2])
    {
    jm=asize-2;
    }
  else 
    {
      jm=(G4int) ((x2max-(*samplingFunction_x)[0])/samplingConstant);
    }
  G4double rumax = (*samplingFunction_y)[jm]+((*samplingFunction_y)[jm+1]-(*samplingFunction_y)[jm])*
    (x2max-(*samplingFunction_x)[jm])/((*samplingFunction_x)[jm+1]-(*samplingFunction_x)[jm]); 
  G4int j,ju,jt;
  G4double ru,denomin,x2rat;
  G4double CDT,G,rand;
  do{
    ru = rumax + std::log(G4UniformRand());
    j=0;
    ju=jm+1;
    do{
     jt=(j+ju)/2; //bipartition
     if (ru > (*samplingFunction_y)[jt])
      {
        j=jt;
      }
     else
      {
        ju=jt;
      }
    }while ((ju-j)>1);
    denomin = (*samplingFunction_y)[j+1]-(*samplingFunction_y)[j];
    if (denomin > 1e-12) 
     {
       x2rat = (*samplingFunction_x)[j]+(((*samplingFunction_x)[j+1]-(*samplingFunction_x)[j])*
					 (ru-(*samplingFunction_y)[j])/denomin)-x2max;
     }
    else
     {
       x2rat = (*samplingFunction_x)[j]-x2max;
     }
    CDT = 1.0-2.0*std::exp(x2rat);
    G = 0.5*(1.0+CDT*CDT);
    rand = G4UniformRand();
   }while (rand>G);
  
  G4double cosTheta = CDT;
  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
 



  // Scattered photon angles. ( Z - axis along the parent photon)
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*std::cos(phi);
  G4double dirY = sinTheta*std::sin(phi);
  G4double dirZ = cosTheta;
  
  // Update G4VParticleChange for the scattered photon 
  G4ThreeVector photonDirection1(dirX, dirY, dirZ);

  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.ProposeEnergy(photonEnergy0);
  aParticleChange.ProposeMomentumDirection(photonDirection1);
  
  aParticleChange.SetNumberOfSecondaries(0);

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4bool G4PenelopeRayleigh::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4PenelopeRayleigh::GetMeanFreePath(const G4Track& track,
					     G4double, // previousStepSize
					      G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy(); 
  material = track.GetMaterial();
  size_t materialIndex = material->GetIndex();

  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else
    {
      meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
    }
  return meanFreePath;
}

void G4PenelopeRayleigh::InizialiseSampling()
{
  SamplingTablePair::iterator theTable=SamplingTables.find(material);
  const G4int points=241;
  G4double Xlow = 0.;
  G4double Xhigh=1e-04;
  G4double fact = std::pow((1e06/Xhigh),(1/240.0));
  samplingConstant=std::log(fact);
  
  if (theTable==SamplingTables.end()) { //material not inizialized yet
    samplingFunction_x = new G4DataVector();
    samplingFunction_y = new G4DataVector();
    G4double sum = 0.0;
    G4PenelopeIntegrator<G4PenelopeRayleigh,G4double(G4PenelopeRayleigh::*)(G4double)> theIntegrator;
    sum = theIntegrator.Calculate(this,&G4PenelopeRayleigh::MolecularFormFactor,
				  Xlow,Xhigh,1e-10); 
    samplingFunction_x->push_back(Xhigh);
    samplingFunction_y->push_back(sum);
    G4int i;
    for (i=1;i<points;i++){
      Xlow=Xhigh;
      Xhigh=Xhigh*fact;
      sum = theIntegrator.Calculate(this,
				    &G4PenelopeRayleigh::MolecularFormFactor,
				    Xlow,Xhigh,1e-10);
      samplingFunction_x->push_back(Xhigh);
      samplingFunction_y->push_back(sum+(*samplingFunction_y)[i-1]);
    }
    for (i=0;i<points;i++){
      (*samplingFunction_x)[i]=std::log((*samplingFunction_x)[i]);
      (*samplingFunction_y)[i]=std::log((*samplingFunction_y)[i]);
    }
    SamplingTables[material] = std::pair<G4DataVector*,G4DataVector*> (samplingFunction_x,samplingFunction_y);
  }
  else { //material already inizialized
    samplingFunction_x=(*theTable).second.first;
    samplingFunction_y=(*theTable).second.second;
  }
}


G4double G4PenelopeRayleigh::MolecularFormFactor(G4double y)
{
  const G4int ntot=95;
  G4double RA1[ntot] = {0.0e0, 3.9265e+0, 4.3100e+1, 5.2757e+1, 2.5021e+1,
		      1.2211e+1, 9.3229e+0, 3.2455e+0, 2.4197e+0, 1.5985e+0,
		      3.0926e+1, 1.5315e+1, 7.7061e+0, 3.9493e+0, 2.2042e+0,
		      1.9453e+1, 1.9354e+1, 8.0374e+0, 8.3779e+1, 5.7370e+1,
		      5.2310e+1, 4.7514e+1, 4.3785e+1, 4.2048e+1, 3.6972e+1,
		      3.3849e+1, 3.1609e+1, 2.8763e+1, 2.7217e+1, 2.4263e+1,
		      2.2403e+1, 1.8606e+1, 1.5143e+1, 1.4226e+1, 1.1792e+1,
		      9.7574e+0, 1.2796e+1, 1.2854e+1, 1.2368e+1, 1.0208e+1,
		      8.2823e+0, 7.4677e+0, 7.6028e+0, 6.1090e+0, 5.5346e+0,
		      4.2340e+0, 4.0444e+0, 4.2905e+0, 4.7950e+0, 5.1112e+0,
		      5.2407e+0, 5.2153e+0, 5.1639e+0, 4.8814e+0, 5.8054e+0,
		      6.6724e+0, 6.5104e+0, 6.3364e+0, 6.2889e+0, 6.3028e+0,
		      6.3853e+0, 6.3475e+0, 6.5779e+0, 6.8486e+0, 7.0993e+0,
		      7.6122e+0, 7.9681e+0, 8.3481e+0, 6.3875e+0, 8.0042e+0,
		      8.0820e+0, 7.6940e+0, 7.1927e+0, 6.6751e+0, 6.1623e+0,
		      5.8335e+0, 5.5599e+0, 4.6551e+0, 4.4327e+0, 4.7601e+0,
		      5.2872e+0, 5.6084e+0, 5.7680e+0, 5.8041e+0, 5.7566e+0,
		      5.6541e+0, 6.3932e+0, 6.9313e+0, 7.0027e+0, 6.8796e+0,
		      6.4739e+0, 6.2405e+0, 6.0081e+0, 5.5708e+0, 5.3680e+0};


  G4double  RA2[ntot] = {0.0e0, 1.3426e-1, 9.4875e+1,-1.0896e+2,-4.5494e+1,
		       -1.9572e+1,-1.2382e+1,-3.6827e+0,-2.4542e+0,-1.4453e+0,
		       1.3401e+2, 7.9717e+1, 6.2164e+1, 4.0300e+1, 3.1682e+1,
		       -1.3639e+1,-1.5950e+1,-5.1523e+0, 1.8351e+2, 1.2205e+2,
		       1.0007e+2, 8.5632e+1, 7.9145e+1, 6.3675e+1, 6.2954e+1,
		       5.6601e+1, 5.4171e+1, 4.8752e+1, 3.8062e+1, 3.9933e+1,
		       4.8343e+1, 4.2137e+1, 3.4617e+1, 2.9430e+1, 2.4010e+1,
		       1.9744e+1, 4.0009e+1, 5.1614e+1, 5.0456e+1, 3.9088e+1,
		       2.6824e+1, 2.2953e+1, 2.4773e+1, 1.6893e+1, 1.4548e+1,
		       9.7226e+0, 1.0192e+1, 1.1153e+1, 1.3188e+1, 1.4733e+1,
		       1.5644e+1, 1.5939e+1, 1.5923e+1, 1.5254e+1, 2.0748e+1,
		       2.6901e+1, 2.7032e+1, 2.4938e+1, 2.1528e+1, 2.0362e+1,
		       1.9474e+1, 1.8238e+1, 1.7898e+1, 1.9174e+1, 1.9023e+1,
		       1.8194e+1, 1.8504e+1, 1.8955e+1, 1.4276e+1, 1.7558e+1,
		       1.8651e+1, 1.7984e+1, 1.6793e+1, 1.5469e+1, 1.4143e+1,
		       1.3149e+1, 1.2255e+1, 9.2352e+0, 8.6067e+0, 9.7460e+0,
		       1.1749e+1, 1.3281e+1, 1.4326e+1, 1.4920e+1, 1.5157e+1,
		       1.5131e+1, 1.9489e+1, 2.3649e+1, 2.4686e+1, 2.4760e+1,
		       2.1519e+1, 2.0099e+1, 1.8746e+1, 1.5943e+1, 1.4880e+1};

  G4double RA3[ntot] =  {0.0e0, 2.2648e+0, 1.0579e+3, 8.6177e+2, 2.4422e+2,
		       7.8788e+1, 3.8293e+1, 1.2564e+1, 6.9091e+0, 3.7926e+0,
		       0.0000e+0, 0.0000e+0, 1.6759e-9, 1.3026e+1, 3.0569e+0,
		       1.5521e+2, 1.2815e+2, 4.7378e+1, 9.2802e+2, 4.7508e+2,
		       3.6612e+2, 2.7582e+2, 2.1008e+2, 1.5903e+2, 1.2322e+2,
		       9.2898e+1, 7.1345e+1, 5.1651e+1, 3.8474e+1, 2.7410e+1,
		       1.9126e+1, 1.0889e+1, 5.3479e+0, 8.2223e+0, 5.0837e+0,
		       2.8905e+0, 2.7457e+0, 6.7082e-1, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 1.7264e-1, 2.7322e-1,
		       3.9444e-1, 4.5648e-1, 6.2286e-1, 7.2468e-1, 8.4296e-1,
		       1.1698e+0, 1.2994e+0, 1.4295e+0, 0.0000e+0, 8.1570e-1,
		       6.9349e-1, 4.9536e-1, 3.1211e-1, 1.5931e-1, 2.9512e-2,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0,
		       0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0, 0.0000e+0};

 G4double RA4[ntot] =  {1.1055e1,6.3519e0,4.7367e+1, 3.9402e+1, 2.2896e+1,
		      1.3979e+1, 1.0766e+1, 6.5252e+0, 5.1631e+0, 4.0524e+0,
		      2.7145e+1, 1.8724e+1, 1.4782e+1, 1.1608e+1, 9.7750e+0,
		      1.6170e+1, 1.5249e+1, 9.1916e+0, 5.4499e+1, 4.1381e+1,
		      3.7395e+1, 3.3815e+1, 3.1135e+1, 2.8273e+1, 2.6140e+1,
		      2.3948e+1, 2.2406e+1, 2.0484e+1, 1.8453e+1, 1.7386e+1,
		      1.7301e+1, 1.5388e+1, 1.3411e+1, 1.2668e+1, 1.1133e+1,
		      9.8081e+0, 1.3031e+1, 1.4143e+1, 1.3815e+1, 1.2077e+1,
		      1.0033e+1, 9.2549e+0, 9.5338e+0, 7.9076e+0, 7.3263e+0,
		      5.9996e+0, 6.0087e+0, 6.2660e+0, 6.7914e+0, 7.1501e+0,
		      7.3367e+0, 7.3729e+0, 7.3508e+0, 7.1465e+0, 8.2731e+0,
		      9.3745e+0, 9.3508e+0, 8.9897e+0, 8.4566e+0, 8.2690e+0,
		      8.1398e+0, 7.9183e+0, 7.9123e+0, 8.1677e+0, 8.1871e+0,
		      8.1766e+0, 8.2881e+0, 8.4227e+0, 7.0273e+0, 8.0002e+0,
		      8.1440e+0, 7.9104e+0, 7.5685e+0, 7.1970e+0, 6.8184e+0,
		      6.5469e+0, 6.3056e+0, 5.4844e+0, 5.2832e+0, 5.5889e+0,
		      6.0919e+0, 6.4340e+0, 6.6426e+0, 6.7428e+0, 6.7636e+0,
		      6.7281e+0, 7.5729e+0, 8.2808e+0, 8.4400e+0, 8.4220e+0,
		      7.8662e+0, 7.5993e+0, 7.3353e+0, 6.7829e+0, 6.5520e+0};

  G4double RA5[ntot] = {0.0e0, 4.9828e+0, 5.5674e+1, 3.0902e+1, 1.1496e+1,
		      4.8936e+0, 2.5506e+0, 1.2236e+0, 7.4698e-1, 4.7042e-1,
		      4.7809e+0, 4.6315e+0, 4.3677e+0, 4.9269e+0, 2.6033e+0,
		      9.6229e+0, 7.2592e+0, 4.1634e+0, 1.3999e+1, 8.6975e+0,
		      6.9630e+0, 5.4681e+0, 4.2653e+0, 3.2848e+0, 2.7354e+0,
		      2.1617e+0, 1.7030e+0, 1.2826e+0, 9.7080e-1, 7.2227e-1,
		      5.0874e-1, 3.1402e-1, 1.6360e-1, 3.2918e-1, 2.3570e-1,
		      1.5868e-1, 1.5146e-1, 9.7662e-2, 7.3151e-2, 6.4206e-2,
		      4.8945e-2, 4.3189e-2, 4.4368e-2, 3.3976e-2, 3.0466e-2,
		      2.4477e-2, 3.7202e-2, 3.7093e-2, 3.8161e-2, 3.8576e-2,
		      3.8403e-2, 3.7806e-2, 3.4958e-2, 3.6029e-2, 4.3087e-2,
		      4.7069e-2, 4.6452e-2, 4.2486e-2, 4.1517e-2, 4.1691e-2,
		      4.2813e-2, 4.2294e-2, 4.5287e-2, 4.8462e-2, 4.9726e-2,
		      5.5097e-2, 5.6568e-2, 5.8069e-2, 1.2270e-2, 3.8006e-2,
		      3.5048e-2, 3.0050e-2, 2.5069e-2, 2.0485e-2, 1.6151e-2,
		      1.4631e-2, 1.4034e-2, 1.1978e-2, 1.1522e-2, 1.2375e-2,
		      1.3805e-2, 1.4954e-2, 1.5832e-2, 1.6467e-2, 1.6896e-2,
		      1.7166e-2, 1.9954e-2, 2.2497e-2, 2.1942e-2, 2.1965e-2,
		      2.0005e-2, 1.8927e-2, 1.8167e-2, 1.6314e-2, 1.5522e-2};

  G4double x=std::sqrt(y);
  G4double gradx1=0.0;
  G4double fa=0.0;

  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();
  const G4int* stechiometric = material->GetAtomsVector();
  const G4double* vector_of_atoms = material->GetVecNbOfAtomsPerVolume();
  const G4double tot_atoms = material->GetTotNbOfAtomsPerVolume();
  for (G4int i=0;i<nElements;i++){
    G4int Z = (G4int) (*elementVector)[i]->GetZ();
    if (Z>ntot) Z=95;
    G4double denomin = 1+y*(RA4[Z-1]+y*RA5[Z-1]);
    fa=Z*(1+y*(RA1[Z-1]+x*(RA2[Z-1]+x*RA3[Z-1])))/(denomin*denomin);
    G4bool a = ((Z>10) && (fa>2.0));
    if (a) 
      {
	G4double Pa,Pg,Pq,fb;
	G4double k1=0.3125;
	G4double k2=2.426311e-02;
	Pa=(Z-k1)*fine_structure_const;
	Pg=std::sqrt(1-(Pa*Pa));
	Pq=k2*x/Pa;
	fb=std::sin(2*Pg*std::atan(Pq))/(Pg*Pq*std::pow((1+Pq*Pq),Pg));
	fa=std::max(fa,fb);
      }
    if (stechiometric && stechiometric[i]!=0)
      {
	gradx1=gradx1+stechiometric[i]*(fa*fa); //sum on the molecule
      }
    else
      {
	gradx1=gradx1+(vector_of_atoms[i]/tot_atoms)*(fa*fa); //weighted mean
      }
  }
  return gradx1;
}


G4double G4PenelopeRayleigh::DifferentialCrossSection(G4double x)
{
  G4double x2=facte*(1-x);
  G4double gradx = (1+x*x)*MolecularFormFactor(x2);
  return gradx;
}


