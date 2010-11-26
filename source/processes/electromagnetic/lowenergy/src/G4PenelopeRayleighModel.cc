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
// $Id: G4PenelopeRayleighModel.cc,v 1.8 2010-11-26 11:51:11 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 14 Oct 2008   L Pandola    Migration from process to model 
// 17 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            PrepareConstants(), since they might be checked later on
// 18 Dec 2009   L Pandola    Added a dummy ComputeCrossSectionPerAtom() method issueing a 
//                            warning if users try to access atomic cross sections via 
//                            G4EmCalculator
//

#include "G4PenelopeRayleighModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4PenelopeIntegrator.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopeRayleighModel::G4PenelopeRayleighModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),samplingFunction_x(0),samplingFunction_xNoLog(0),
   theMaterial(0),
   isInitialised(false)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  PrepareConstants();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeRayleighModel::~G4PenelopeRayleighModel()
{
  std::map <const G4Material*,G4DataVector*>::iterator i;
  for(i=SamplingTable.begin(); i != SamplingTable.end(); i++) {
    delete (*i).second;
  }
  if (samplingFunction_x) delete samplingFunction_x;
  if (samplingFunction_xNoLog) delete samplingFunction_xNoLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::Initialise(const G4ParticleDefinition* ,
					 const G4DataVector& )
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeRayleighModel::Initialise()" << G4endl;


  if (verboseLevel > 0) {
    G4cout << "Penelope Rayleigh model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / keV << " keV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4PenelopeRayleighModel::CrossSectionPerVolume(const G4Material* material,
					       const G4ParticleDefinition* p,
					       G4double ekin,
					       G4double,
					       G4double)
{
  // Penelope model to calculate the Rayleigh scattering inverse mean 
  // free path. 
  //
  // The basic method is from
  //  M. Born, Atomic physics, Ed. Blackie and Sons (1969)
  // using numerical approximations developed in
  //  J. Baro' et al., Radiat. Phys. Chem. 44 (1994) 531
  // Data used for the form factor used in the calculation (numerical integral of 
  // dSigma/dOmega) are been derived by fitting the atomic forn factor tables 
  // tabulated in 
  //  J.H. Hubbel et al., J. Phys. Chem. Ref. Data 4 (1975) 471; erratum ibid. 
  //   6 (1977) 615.
  // The numerical integration of the differential cross section dSigma/dOmega, 
  // which is implemented in the DifferentialCrossSection() method, is performed 
  // by 20-point Gaussian method (managed by G4PenelopeIntegrator). 
  //

  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4PenelopeRayleighModel" << G4endl;
  SetupForMaterial(p, material, ekin);
  
  //Assign local variable "material" to private member "theMaterial", because 
  //this information is necessary to calculate the cross section
  theMaterial = material;

  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  G4int maxZ=0;
  for (G4int i=0; i<nElements; i++) 
    {
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      if (Z>maxZ)
	maxZ = Z;
    }
      
  G4double ec=std::min(ekin,0.5*maxZ);
  factorE = 849.3315*(ec/electron_mass_c2)*(ec/electron_mass_c2);
  G4double cs=0;
  //
  //Integrate the Differential Cross Section dSigma/dCosTheta between -1 and 1.
  //
  G4PenelopeIntegrator<G4PenelopeRayleighModel,G4double(G4PenelopeRayleighModel::*)(G4double)> 
    theIntegrator;
  cs =
    theIntegrator.Calculate(this,&G4PenelopeRayleighModel::DifferentialCrossSection,
			    -1.0,0.90,1e-06);
  cs += theIntegrator.Calculate(this,&G4PenelopeRayleighModel::DifferentialCrossSection,
				0.90,0.9999999,1e-06);
  cs = cs*(ec/ekin)*(ec/ekin)*pi*classic_electr_radius*classic_electr_radius;
  //
  // Here cs represents the cross section per molecule for materials defined with chemical 
  // formulas and the average cross section per atom in compounds (defined with the mass 
  // fraction)
  //
  const G4int* stechiometric = material->GetAtomsVector();
  //This is the total density of atoms in the material
  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double moleculeDensity = 0;

  //Default case: the material is a compound. In this case, cs is the average cross section 
  // _per_atom_ and one has to multiply for the atom density.
  G4double cross = atomDensity*cs;
  G4bool isAMolecule = false;
  
  //Alternative case: the material is a molecule. In this case cs is the cross section 
  // _per_molecule_ and one has to multiply for the molecule density
  if (stechiometric)
    {
      //Calculate the total number of atoms per molecule
      G4int atomsPerMolecule = 0;
      for (G4int k=0;k<nElements;k++)
       {
	atomsPerMolecule += stechiometric[k];
        if (verboseLevel > 2)
	   {
	     G4cout << "Element: " << (G4int) (*elementVector)[k]->GetZ() << " has " << 
	        stechiometric[k] << " atoms/molecule" << G4endl;	
           }
       }
      if (atomsPerMolecule)
	{
	  isAMolecule = true;
	  if (verboseLevel > 3)
	    {
	      G4cout << "Material " << material->GetName() << " is a molecule composed by " << 
		atomsPerMolecule << " atoms" << G4endl;
	    }
	  moleculeDensity = atomDensity/((G4double) atomsPerMolecule);
	  cross = cs*moleculeDensity;
	}
    }
  
  if (verboseLevel > 2)
    {
      if (isAMolecule)
	{
	  G4cout << "Rayleigh cross section at " << ekin/keV << " keV for material " << 
	    material->GetName() << " (molecule) = " << cs/barn << " barn/molecule." << G4endl;
	  G4cout << "Mean free path: " << (1./cross)/mm << " mm" << G4endl;
	}
      else
	{
	  G4cout << "Rayleigh cross section at " << ekin/keV << " keV for material " << 
	    material->GetName() << " (compound) = " << cs/barn << " barn/atom." << G4endl;
	  G4cout << "Mean free path: " << (1./cross)/mm << " mm" << G4endl;
	}
    }
  return cross;
}


//This is a dummy method. Never inkoved by the tracking, it just issues
//a warning if one tries to get Cross Sections per Atom via the
//G4EmCalculator.
G4double G4PenelopeRayleighModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
							     G4double,
							     G4double,
							     G4double,
							     G4double,
							     G4double)
{
  G4cout << "*** G4PenelopeRayleighModel -- WARNING ***" << G4endl;
  G4cout << "Penelope Rayleigh model does not calculate cross section _per atom_ " << G4endl;
  G4cout << "so the result is always zero. For physics values, please invoke " << G4endl;
  G4cout << "GetCrossSectionPerVolume() or GetMeanFreePath() via the G4EmCalculator" << G4endl;
  return 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* ,
						const G4MaterialCutsCouple* couple,
						const G4DynamicParticle* aDynamicGamma,
						G4double,
						G4double)
{
  // Penelope model to sample the Rayleigh scattering final state.
  //
  // The angular deflection of the scattered photon is sampled according to the 
  // differential cross section dSigma/dOmega used for the numerical integration, 
  // and implemented in the DifferentialCrossSection() method. See comments in 
  // method CrossSectionPerVolume() for more details on the original references 
  // of the model.
  //

  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeRayleighModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();
 
  if (photonEnergy0 <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();
   
  // Sampling inizialitation (build internal table)
  theMaterial = couple->GetMaterial();
  InitialiseSampling();

  G4DataVector* samplingFunction_y = SamplingTable.find(theMaterial)->second;

   // Sample the angle of the scattered photon
  G4double x2max = 2.0*std::log(41.2148*photonEnergy0/electron_mass_c2);
  G4int jm=0;
  G4int asize = samplingFunction_x->size();
  if (x2max < (*samplingFunction_x)[1]) 
    jm=0;
  else if(x2max>(*samplingFunction_x)[asize-2])
    jm=asize-2;
  else 
    {
      //spacing in the logTable
      G4double logScalingFactor = (*samplingFunction_x)[1]-(*samplingFunction_x)[0];
      jm=(G4int) ((x2max-(*samplingFunction_x)[0])/logScalingFactor);
    }

  G4double rumax = (*samplingFunction_y)[jm]+((*samplingFunction_y)[jm+1]-(*samplingFunction_y)[jm])*
    (x2max-(*samplingFunction_x)[jm])/((*samplingFunction_x)[jm+1]-(*samplingFunction_x)[jm]); 
  G4double cosTheta=0;
  G4double rejectionValue = 0;
  do{
    G4double ru = rumax + std::log(G4UniformRand());
    G4int j=0;
    G4int ju=jm+1;
    do{
      G4int jt=(j+ju)/2; //bipartition
      if (ru > (*samplingFunction_y)[jt])
	j=jt;
      else
	ju=jt;
    }while ((ju-j)>1);
    G4double x2rat = 0;
    G4double denomin = (*samplingFunction_y)[j+1]-(*samplingFunction_y)[j];
    if (denomin > 1e-12) 
     {
       x2rat = (*samplingFunction_x)[j]+(((*samplingFunction_x)[j+1]-(*samplingFunction_x)[j])*
					 (ru-(*samplingFunction_y)[j])/denomin)-x2max;
     }
    else
     {
       x2rat = (*samplingFunction_x)[j]-x2max;
     }
    cosTheta = 1.0-2.0*std::exp(x2rat);
    rejectionValue = 0.5*(1.0+cosTheta*cosTheta);
   }while (G4UniformRand() > rejectionValue);

  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
 
  // Scattered photon angles. ( Z - axis along the parent photon)
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*std::cos(phi);
  G4double dirY = sinTheta*std::sin(phi);
  G4double dirZ = cosTheta;
  
  // Update G4VParticleChange for the scattered photon 
  G4ThreeVector photonDirection1(dirX, dirY, dirZ);

  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;
  fParticleChange->SetProposedKineticEnergy(photonEnergy0) ;
   
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeRayleighModel::DifferentialCrossSection(G4double cosTheta)
{
  //Differential cross section for Rayleigh scattering
  G4double x2 = factorE*(1-cosTheta);
  G4double gradx = (1+cosTheta*cosTheta)*MolecularFormFactor(x2);
  return gradx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeRayleighModel::MolecularFormFactor(G4double x2)
{
  //Squared molecular form factor (additivity rule)
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

  G4double x=std::sqrt(x2);
  G4double gradx1=0.0;

  G4int nElements = theMaterial->GetNumberOfElements();
  const G4ElementVector* elementVector = theMaterial->GetElementVector();
  const G4int* stechiometric = theMaterial->GetAtomsVector();
  const G4double* vector_of_atoms = theMaterial->GetVecNbOfAtomsPerVolume();
  const G4double tot_atoms = theMaterial->GetTotNbOfAtomsPerVolume();
  for (G4int i=0;i<nElements;i++)
    {
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      if (Z>ntot) Z=95;
      G4double denomin = 1.+x2*(RA4[Z-1]+x2*RA5[Z-1]);
      G4double fa=Z*(1+x2*(RA1[Z-1]+x*(RA2[Z-1]+x*RA3[Z-1])))/(denomin*denomin);
      if (Z>10 && fa>2.0)
	{
	  G4double k1=0.3125;
	  G4double k2=2.426311e-02;
	  G4double Pa=(Z-k1)*fine_structure_const;
	  G4double Pg=std::sqrt(1-(Pa*Pa));
	  G4double Pq=k2*x/Pa;
	  G4double fb=std::sin(2.0*Pg*std::atan(Pq))/(Pg*Pq*std::pow((1+Pq*Pq),Pg));
	  fa=std::max(fa,fb);
	}
      if (stechiometric && stechiometric[i]!=0)
	gradx1 += stechiometric[i]*(fa*fa); //sum on the molecule
      else
	gradx1 += (vector_of_atoms[i]/tot_atoms)*(fa*fa); //weighted mean
    }
  return gradx1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::InitialiseSampling()
{
  if (!samplingFunction_x || !samplingFunction_xNoLog)
    {
      G4cout << "G4PenelopeRayleighModel::InitialiseSampling(), something wrong" << G4endl;
      G4cout << "It looks like G4PenelopeRayleighModel::PrepareConstants() has not been called" 
	     << G4endl;
      G4Exception();
      return;
    }
  if (!SamplingTable.count(theMaterial)) //material not defined yet
    { 
      G4double XlowInte = 0.;
      G4double XhighInte=(*samplingFunction_xNoLog)[0];
      //material not inizialized yet: initialize it now
      G4DataVector* samplingFunction_y = new G4DataVector();
      G4PenelopeIntegrator<G4PenelopeRayleighModel,G4double(G4PenelopeRayleighModel::*)(G4double)> 
	theIntegrator;
      G4double sum = theIntegrator.Calculate(this,&G4PenelopeRayleighModel::MolecularFormFactor,
				    XlowInte,XhighInte,1e-10); 
      samplingFunction_y->push_back(sum);
      for (G4int i=1;i<nPoints;i++)
	{
	  XlowInte=(*samplingFunction_xNoLog)[i-1];
	  XhighInte=(*samplingFunction_xNoLog)[i];
	  sum += theIntegrator.Calculate(this,
					&G4PenelopeRayleighModel::MolecularFormFactor,
					XlowInte,XhighInte,1e-10);
	  samplingFunction_y->push_back(sum);
	}
      for (G4int i=0;i<nPoints;i++)	
	(*samplingFunction_y)[i]=std::log((*samplingFunction_y)[i]);
	
      //
      /*
      G4String nnn = theMaterial->GetName()+".ntab";
      std::ofstream file(nnn);
      for (size_t k=0;k<samplingFunction_x->size();k++)
	file << (*samplingFunction_x)[k] << " " << (*samplingFunction_y)[k] << G4endl;
      file.close();
      */
      //
      SamplingTable[theMaterial] = samplingFunction_y;
  }
}

void G4PenelopeRayleighModel::PrepareConstants()
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeRayleighModel::PrepareConstants()" << G4endl;
  nPoints=241;
  Xlow=1e-04;
  Xhigh=1e06;
  if (samplingFunction_x)
    {
      delete samplingFunction_x;
      samplingFunction_x = 0;
    }
  if (samplingFunction_xNoLog)
    {
      delete samplingFunction_xNoLog;
      samplingFunction_xNoLog = 0;
    }

  samplingFunction_x = new G4DataVector();
  samplingFunction_xNoLog = new G4DataVector();
 
  G4double scalingFactor = std::pow((Xhigh/Xlow),((1/((G4double)(nPoints-1)))));
  G4double logScalingFactor=(1/((G4double)(nPoints-1)))*std::log(Xhigh/Xlow);
  //Logarithmic table between log(Xlow) and log(Xhigh)
  samplingFunction_x->push_back(std::log(Xlow));
  //Table between Xlow and Xhigh with log spacement: needed for InitialiseSampling()
  samplingFunction_xNoLog->push_back(Xlow);

  for (G4int i=1;i<nPoints;i++)
    {
      G4double nextx = (*samplingFunction_x)[i-1]+logScalingFactor;
      samplingFunction_x->push_back(nextx);
      nextx = (*samplingFunction_xNoLog)[i-1]*scalingFactor;
      samplingFunction_xNoLog->push_back(nextx);
    }
}
