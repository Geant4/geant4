//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LowEnergyPolarizedCompton.cc,v 1.12 2002-05-31 18:08:35 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//

// --------- G4LowEnergyPolarizedCompton class -----
//
//           by G.Depaola & F.Longo (21 may 2001)
//
// 21 May 2001 - MGP      Modified to inherit from G4VDiscreteProcess
//                        Applies same algorithm as LowEnergyCompton
//                        if the incoming photon is not polarised
//                        Temporary protection to avoid crash in the case 
//                        of polarisation || incident photon direction
//
// 17 October 2001 - F.Longo - Revised according to a design iteration
//
// 21 February 2002 - F.Longo Revisions with A.Zoglauer and G.Depaola
//                            - better description of parallelism
//                            - system of ref change method improved
//
//
// ************************************************************
//
// Corrections by Rui Curado da Silva (2000)
// New Implementation by G.Depaola & F.Longo      
//
// - sampling of phi
// - polarization of scattered photon
//
// --------------------------------------------------------------

#include "G4LowEnergyPolarizedCompton.hh"
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
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VRangeTest.hh"
#include "G4RangeTest.hh"

#include "G4CutsPerMaterialWarning.hh"

// constructor

G4LowEnergyPolarizedCompton::G4LowEnergyPolarizedCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit (250*eV),              // initialization
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyPolarizedCompton::G4LowEnergyPolarizedCompton - energy outside intrinsic process validity range");
    }
  
  crossSectionHandler = new G4CrossSectionHandler;
  
  
  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new 
    G4CompositeEMDataSet(scatterFile,scatterInterpolation,1.,1.);
  
  meanFreePathTable = 0;
  
  rangeTest = new G4RangeTest;
  
   if (verboseLevel > 0) 
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: " 
	      << lowEnergyLimit / keV << " keV - "
	      << highEnergyLimit / GeV << " GeV" 
	      << G4endl;
     }
}

// destructor

G4LowEnergyPolarizedCompton::~G4LowEnergyPolarizedCompton()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete scatterFunctionData;
  delete rangeTest;
}
 

void G4LowEnergyPolarizedCompton::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  
  G4CutsPerMaterialWarning warning;
  warning.PrintWarning(&photon);

  crossSectionHandler->Clear();
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);
  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyPolarizedCompton::PostStepDoIt(const G4Track& aTrack,
							     const G4Step&  aStep)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
  // GEANT4 internal units
  //
  // Note : Effects due to binding of atomic electrons are negliged.

  aParticleChange.Initialize(aTrack);

  // Dynamic particle quantities  
  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double gammaEnergy0 = incidentPhoton->GetKineticEnergy();
  G4ThreeVector gammaPolarization0 = incidentPhoton->GetPolarization(); 

  //  gammaPolarization0 = gammaPolarization0.unit(); // 

  // Protection: a polarisation parallel to the 
  // direction causes problems; 
  // in that case find a random polarization
  
  G4ThreeVector gammaDirection0 = incidentPhoton->GetMomentumDirection();
  G4double scalarproduct = gammaPolarization0.dot(gammaDirection0);
  G4double angle = gammaPolarization0.angle(gammaDirection0);

  // Make sure that the polarization vector is perpendicular to the 
  // gamma direction. If not 
  
  if(!(gammaPolarization0.isOrthogonal(gammaDirection0, 1e-6))||(gammaPolarization0.mag()==0)) 
    { // only for testing now
      gammaPolarization0 = GetRandomPolarization(gammaDirection0);
    }
  else
    {
      if ( gammaPolarization0.howOrthogonal(gammaDirection0) != 0)
	{
	  gammaPolarization0 = GetPerpendicularPolarization(gammaDirection0, gammaPolarization0);
	}
    }

  // End of Protection

  // Within energy limit?
  
  if(gammaEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(gammaEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }
  
  G4double E0_m = gammaEnergy0 / electron_mass_c2 ;

  // Select randomly one element in the current material
  
  G4Material* material = aTrack.GetMaterial();
  G4int Z = crossSectionHandler->SelectRandomAtom(material,gammaEnergy0);

  // Sample the energy and the polarization of the scattered photon 
  
  G4double epsilon, epsilonSq, onecost, sinThetaSqr, greject ;
  
  G4double epsilon0 = 1./(1. + 2*E0_m);
  G4double epsilon0Sq = epsilon0*epsilon0;
  G4double alpha1   = - log(epsilon0);
  G4double alpha2 = 0.5*(1.- epsilon0Sq);

  G4double wlGamma = h_Planck*c_light/gammaEnergy0;
  G4double gammaEnergy1;
  G4ThreeVector gammaDirection1;
  
  do {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand() )
      { 
	epsilon   = exp(-alpha1*G4UniformRand());  
	epsilonSq = epsilon*epsilon; 
      }
    else 
      {
	epsilonSq = epsilon0Sq + (1.- epsilon0Sq)*G4UniformRand();
	epsilon   = sqrt(epsilonSq);
      }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sinThetaSqr   = onecost*(2.-onecost);
    
    // Protection
    if (sinThetaSqr > 1.)
      {
	if (verboseLevel>0) G4cout 
	  << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	  << "sin(theta)**2 = " 
	  << sinThetaSqr
	  << "; set to 1" 
	  << G4endl;
	sinThetaSqr = 1.;
      }
    if (sinThetaSqr < 0.)
      {
	if (verboseLevel>0) G4cout 
	  << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	  << "sin(theta)**2 = " 
	  << sinThetaSqr
	  << "; set to 0" 
	  << G4endl;
	sinThetaSqr = 0.;
      }
    // End protection

    G4double x =  sqrt(onecost/2.) / (wlGamma/cm);;
    G4double scatteringFunction = scatterFunctionData->FindValue(x,Z-1);
    greject = (1. - epsilon*sinThetaSqr/(1.+ epsilonSq))*scatteringFunction;    
        
  } while(greject < G4UniformRand()*Z);
    
  
  // ****************************************************
  //		Phi determination
  // ****************************************************
  
  G4double phi = SetPhi(epsilon,sinThetaSqr);
  
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //
  
  G4double cosTheta = 1. - onecost;
  
  // Protection
  
  if (cosTheta > 1.)
    {
      if (verboseLevel>0) G4cout 
	<< " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	<< "cosTheta = " 
	<< cosTheta
	<< "; set to 1" 
	<< G4endl;
      cosTheta = 1.;
    }
  if (cosTheta < -1.)
    {
      if (verboseLevel>0) G4cout 
	<< " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	<< "cosTheta = " 
	<< cosTheta
	<< "; set to -1" 
	<< G4endl;
      cosTheta = -1.;
    }
  // End protection      
  
  
  G4double sinTheta = sqrt (sinThetaSqr);
  
  // Protection
  if (sinTheta > 1.)
    {
      if (verboseLevel>0) G4cout 
	<< " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	<< "sinTheta = " 
	<< sinTheta
	<< "; set to 1" 
	<< G4endl;
      sinTheta = 1.;
    }
  if (sinTheta < -1.)
    {
      if (verboseLevel>0) G4cout 
	<< " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	<< "sinTheta = " 
	<< sinTheta
	<< "; set to -1" 
	<< G4endl;
      sinTheta = -1.;
    }
  // End protection
  
      
  G4double dirx = sinTheta*cos(phi);
  G4double diry = sinTheta*sin(phi);
  G4double dirz = cosTheta ;
  
  //
  // update G4VParticleChange for the scattered photon 
  //
  
  gammaEnergy1 = epsilon*gammaEnergy0;
      
  // New polarization
  
  G4ThreeVector gammaPolarization1 = SetNewPolarization(epsilon,
							sinThetaSqr,
							phi,
							cosTheta);
  
  // Set new direction 
  G4ThreeVector tmpDirection1( dirx,diry,dirz );
  gammaDirection1 = tmpDirection1;
      
  // Change reference frame.
  
  SystemOfRefChange(gammaDirection0,gammaDirection1,
		    gammaPolarization0,gammaPolarization1);   
  
  if (gammaEnergy1 > 0.)
    {
      aParticleChange.SetEnergyChange( gammaEnergy1 ) ;
      aParticleChange.SetMomentumChange( gammaDirection1 );
      aParticleChange.SetPolarizationChange( gammaPolarization1 );
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
    }
  
  //
  // kinematic of the scattered electron
  //
  
  G4double ElecKineEnergy = gammaEnergy0 - gammaEnergy1 ;
  
  // Generate the electron only if with large enough range w.r.t. cuts and safety
  
  G4double safety = aStep.GetPostStepPoint()->GetSafety();
  
  if (rangeTest->Escape(G4Electron::Electron(),material,ElecKineEnergy,safety))
    {
      G4double ElecMomentum = sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));
      G4ThreeVector ElecDirection((gammaEnergy0 * gammaDirection0 - 
				   gammaEnergy1 * gammaDirection1) * (1./ElecMomentum));  
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),ElecDirection.unit(),ElecKineEnergy) ;
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(electron);
      aParticleChange.SetLocalEnergyDeposit(0.); 
    }
  else
    {
      aParticleChange.SetNumberOfSecondaries(0);
      aParticleChange.SetLocalEnergyDeposit(ElecKineEnergy);
    }
  
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
  
}


G4double G4LowEnergyPolarizedCompton::SetPhi(G4double energyRate,
					     G4double sinSqrTh)
{
  G4double rand1;
  G4double rand2;
  G4double phiProbability;
  G4double phi;
  G4double a, b;
  
  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      phiProbability=0.;
      phi = twopi*rand1;
      
      a = 2*sinSqrTh;
      b = energyRate + 1/energyRate;
      
      phiProbability = 1 - (a/b)*(cos(phi)*cos(phi));

      
 
    }
  while ( rand2 > phiProbability );
  return phi;
}


G4ThreeVector G4LowEnergyPolarizedCompton::SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}

G4ThreeVector G4LowEnergyPolarizedCompton::GetRandomPolarization(G4ThreeVector& direction0)
{
  G4ThreeVector d0 = direction0.unit();
  G4ThreeVector a1 = SetPerpendicularVector(d0); //different orthogonal
  G4ThreeVector a0 = a1.unit(); // unit vector

  G4double rand1 = G4UniformRand();
  
  G4double angle = twopi*rand1; // random polar angle
  G4ThreeVector b0 = d0.cross(a0); // cross product
  
  G4ThreeVector c;
  
  c.setX(cos(angle)*(a0.x())+sin(angle)*b0.x());
  c.setY(cos(angle)*(a0.y())+sin(angle)*b0.y());
  c.setZ(cos(angle)*(a0.z())+sin(angle)*b0.z());
  
  G4ThreeVector c0 = c.unit();

  return c0;
  
}


G4ThreeVector G4LowEnergyPolarizedCompton::GetPerpendicularPolarization
(const G4ThreeVector& gammaDirection, const G4ThreeVector& gammaPolarization) const
{

  // 
  // The polarization of a photon is always perpendicular to its momentum direction.
  // Therefore this function removes those vector component of gammaPolarization, which
  // points in direction of gammaDirection
  //
  // Mathematically we search the projection of the vector a on the plane E, where n is the
  // plains normal vector.
  // The basic equation can be found in each geometry book (e.g. Bronstein):
  // p = a + (a o n)/(n o n)*n
  
  return gammaPolarization + gammaPolarization.dot(gammaDirection)/gammaDirection.dot(gammaDirection) * gammaDirection;
}


G4ThreeVector G4LowEnergyPolarizedCompton::SetNewPolarization(G4double epsilon, 
							      G4double sinSqrTh, 
							      G4double phi,
							      G4double costheta) 
{
  G4double rand1;
  G4double rand2;
  G4double cosPhi = cos(phi);
  G4double sinPhi = sin(phi);
  G4double sinTheta = sqrt(sinSqrTh);
  G4double cosSqrPhi = cosPhi*cosPhi;
  //  G4double cossqrth = 1.-sinSqrTh;
  //  G4double sinsqrphi = sinPhi*sinPhi;
  G4double normalisation = sqrt(1. - cosSqrPhi*sinSqrTh);
  

  // Determination of Theta 
  
  G4double thetaProbability;
  G4double theta;
  G4double a, b;
  G4double cosTheta;

  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      thetaProbability=0.;
      theta = twopi*rand1;
      a = 4*normalisation*normalisation;
      b = (epsilon + 1/epsilon) - 2;
      thetaProbability = (b + a*cos(theta)*cos(theta))/(a+b);
      cosTheta = cos(theta);       
    }
  while ( rand2 > thetaProbability );
  
  G4double cosBeta = cosTheta;
  G4double sinBeta = sqrt(1-cosBeta*cosBeta);
  
  G4ThreeVector gammaPolarization1;
  
  G4double xParallel = normalisation*cosBeta;
  G4double yParallel = -(sinSqrTh*cosPhi*sinPhi)*cosBeta/normalisation;
  G4double zParallel = -(costheta*sinTheta*cosPhi)*cosBeta/normalisation;
  G4double xPerpendicular = 0.;
  G4double yPerpendicular = (costheta)*sinBeta/normalisation;
  G4double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;
  
  G4double xTotal = (xParallel + xPerpendicular);
  G4double yTotal = (yParallel + yPerpendicular);
  G4double zTotal = (zParallel + zPerpendicular);
  
  gammaPolarization1.setX(xTotal);
  gammaPolarization1.setY(yTotal);
  gammaPolarization1.setZ(zTotal);
  
  return gammaPolarization1;

}


void G4LowEnergyPolarizedCompton::SystemOfRefChange
(G4ThreeVector& direction0,G4ThreeVector& direction1,
 G4ThreeVector& polarization0,G4ThreeVector& polarization1)
  
{
  // direction0 is the original photon direction ---> z
  // polarization0 is the original photon polarization ---> x
  // need to specify y axis in the real reference frame ---> y 
  G4ThreeVector Axis_Z0 = direction0.unit();
  G4ThreeVector Axis_X0 = polarization0.unit();
  G4ThreeVector Axis_Y0 = (Axis_Z0.cross(Axis_X0)).unit(); // to be confirmed;

  G4double direction_x = direction1.getX();
  G4double direction_y = direction1.getY();
  G4double direction_z = direction1.getZ();
  
  direction1 = (direction_x*Axis_X0 + direction_y*Axis_Y0 +  direction_z*Axis_Z0).unit();
  
  G4double polarization_x = polarization1.getX();
  G4double polarization_y = polarization1.getY();
  G4double polarization_z = polarization1.getZ();
  
  polarization1 =(polarization_x*Axis_X0+polarization_y*Axis_Y0+polarization_z*Axis_Z0).unit();

}


G4bool G4LowEnergyPolarizedCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}


G4double G4LowEnergyPolarizedCompton::GetMeanFreePath(const G4Track& track, 
						      G4double previousStepSize,
						      G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  size_t materialIndex = material->GetIndex();
  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
  return meanFreePath;
}














