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
// --------------------------------------------------------------------
//
// $Id: G4LowEnergyPolarizedRayleigh.cc,v 1.1 2003-05-16 08:58:56 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//         Gerardo O. Depaola 
//
// History:
// -------- 
// Added Livermore data table construction methods A. Forti
// Added BuildMeanFreePath A. Forti
// Added PostStepDoIt A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements  A.Forti
// 24.04.01 V.Ivanchenko remove RogueWave 
// 11.08.2001 MGP - Major revision according to a design iteration
// 06.10.2001 MGP - Added strategy to test range for secondary generation
// 05.06.2002 F.Longo and G.Depaola  - bug fixed in angular distribution
// 20.10.2002 G. Depaola - Change sampling method of theta
// 08.02.2003 G. Depaola - Introduce Polarization from G4LowEnergyRayleigh
// 25.04.2003 G. Depaola & F.Longo - cuts per region 
//
// --------------------------------------------------------------------

#include "G4LowEnergyPolarizedRayleigh.hh"
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
#include "G4MaterialCutsCouple.hh"

G4LowEnergyPolarizedRayleigh::G4LowEnergyPolarizedRayleigh(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),            // initialization 
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyPolarizedRayleigh::G4LowEnergyPolarizedRayleigh - energy limit outside intrinsic process validity range");
    }
  
  crossSectionHandler = new G4CrossSectionHandler();

  G4VDataSetAlgorithm* ffInterpolation = new G4LogLogInterpolation;
  G4String formFactorFile = "rayl/re-ff-";
  formFactorData = new G4CompositeEMDataSet(formFactorFile,ffInterpolation,1.,1.);

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
 
G4LowEnergyPolarizedRayleigh::~G4LowEnergyPolarizedRayleigh()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete formFactorData;
  delete rangeTest;
}

void G4LowEnergyPolarizedRayleigh::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  
  crossSectionHandler->Clear();
  G4String crossSectionFile = "rayl/re-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyPolarizedRayleigh::PostStepDoIt(const G4Track& aTrack, 
						     const G4Step& aStep)
{

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();
  G4ThreeVector gammaPolarization0 = incidentPhoton->GetPolarization(); 

  // Protection: a polarisation parallel to the direction causes problems; 
  // in that case find a random polarization
  
  G4ThreeVector photonDirection0 = incidentPhoton->GetMomentumDirection();

  // Make sure that the polarization vector is perpendicular to the 
  // gamma direction. If not 
  
  if(!(gammaPolarization0.isOrthogonal(photonDirection0, 1e-6))||(gammaPolarization0.mag()==0)) 
    { // only for testing now
      gammaPolarization0 = GetRandomPolarization(photonDirection0);
    }
  else
    {
      if ( gammaPolarization0.howOrthogonal(photonDirection0) != 0)
	{
	  gammaPolarization0 = GetPerpendicularPolarization(photonDirection0, gammaPolarization0);
	}
    }
  
  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  //  G4double e0m = photonEnergy0 / electron_mass_c2 ;
//  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  // Select randomly one element in the current material

  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);

  // Sample the energy of the scattered photon 

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  G4double gReject,x,dataFormFactor;
  G4double randomFormFactor;
  G4double cosTheta;
  G4double sinTheta;
  G4double costheta;
  G4double fcostheta;
  
  do
    {
      do
      {
      costheta = 2. * G4UniformRand() - 1.;
      fcostheta = ( 1. + costheta*costheta)/2.;
      } while (fcostheta < G4UniformRand());

      G4double sinThetaHalf = sqrt((1. - costheta) / 2.);
      x = sinThetaHalf / (wlPhoton/cm);
      if (x > 1.e+005) 
         dataFormFactor = formFactorData->FindValue(x,Z-1);
      else
         dataFormFactor = formFactorData->FindValue(0.,Z-1);
      randomFormFactor = G4UniformRand() * Z * Z;
      cosTheta = costheta;
      sinTheta = sqrt(1. - costheta*costheta);
      gReject = dataFormFactor * dataFormFactor;
      
    } while( gReject < randomFormFactor);

  // Scattered photon angles. ( Z - axis along the parent photon)

  G4double sinThetaSqr = sinTheta*sinTheta;
  G4double phi = SetPhi(sinThetaSqr);

  
  G4double dirX = sinTheta*cos(phi);
  G4double dirY = sinTheta*sin(phi);
  G4double dirZ = cosTheta;
  
  G4ThreeVector gammaPolarization1 = SetNewPolarization(sinThetaSqr,
							phi,
							cosTheta);
  
  // Set new direction 
  G4ThreeVector tmpDirection1( dirX,dirY,dirZ );
  G4ThreeVector gammaDirection1 = tmpDirection1;
      
  // Change reference frame.
  
  SystemOfRefChange(photonDirection0,gammaDirection1,
		    gammaPolarization0,gammaPolarization1);   
  
  if (photonEnergy0 > 0.)
    {
      aParticleChange.SetEnergyChange( photonEnergy0 ) ;
      aParticleChange.SetMomentumChange( gammaDirection1 );
      aParticleChange.SetPolarizationChange( gammaPolarization1 );
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
    }

  
  aParticleChange.SetNumberOfSecondaries(0);

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4bool G4LowEnergyPolarizedRayleigh::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyPolarizedRayleigh::GetMeanFreePath(const G4Track& track, 
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
G4double G4LowEnergyPolarizedRayleigh::SetPhi(G4double sinSqrTh)
{
  G4double rand1;
  G4double rand2;
  G4double phiProbability;
  G4double phi;
  
  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      phiProbability=0.;
      phi = twopi*rand1;
      
      phiProbability = 1 - (sinSqrTh*sinSqrTh)*(cos(phi)*cos(phi)); 
    }
  while ( rand2 > phiProbability );
  return phi;
}


G4ThreeVector G4LowEnergyPolarizedRayleigh::SetPerpendicularVector(G4ThreeVector& a)
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

G4ThreeVector G4LowEnergyPolarizedRayleigh::GetRandomPolarization(G4ThreeVector& direction0)
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

G4ThreeVector G4LowEnergyPolarizedRayleigh::GetPerpendicularPolarization
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


G4ThreeVector G4LowEnergyPolarizedRayleigh::SetNewPolarization(G4double sinSqrTh, 
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
  G4double cosTheta;

  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      thetaProbability=0.;
      theta = twopi*rand1;
      thetaProbability = cos(theta)*cos(theta);
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


void G4LowEnergyPolarizedRayleigh::SystemOfRefChange
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







