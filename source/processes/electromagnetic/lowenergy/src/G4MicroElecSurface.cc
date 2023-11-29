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
// G4MicroElecSurface.cc, 
//                   	    	2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//				       	   Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//				           M. Raine and D. Lambert are with CEA [a]
//
// A part of this work has been funded by the French space agency(CNES[c])
// [a] CEA, DAM, DIF - 91297 ARPAJON, France
// [b] ONERA - DPHY, 2 avenue E.Belin, 31055 Toulouse, France
// [c] CNES, 18 av.E.Belin, 31401 Toulouse CEDEX, France
//
// Based on the following publications
//
//	- Q.Gibaru, C.Inguimbert, P.Caron, M.Raine, D.Lambert, J.Puech, 
//	      Geant4 physics processes for microdosimetry and secondary electron emission simulation : 
//	      Extension of MicroElec to very low energies and new materials
//	      NIM B, 2020, in review.
//
//
//      - Modèle de transport d'électrons à basse énergie (10 eV- 2 keV) pour
//        applications spatiales (OSMOSEE, GEANT4), PhD dissertation, 2017. 
//      
//
////////////////////////////////////////////////////////////////////////

#include "G4MicroElecSurface.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4EmProcessSubType.hh"
#include "G4GeometryTolerance.hh"
#include "G4EmProcessSubType.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecSurface::G4MicroElecSurface(const G4String& processName,G4ProcessType type)
  : G4VDiscreteProcess(processName, type),
    oldMomentum(0.,0.,0.), previousMomentum(0.,0.,0.),
    theGlobalNormal(0.,0.,0.), theFacetNormal(0.,0.,0.)
{
  if ( verboseLevel > 0) 
    {
      G4cout << GetProcessName() << " is created " << G4endl;
    }
  
  isInitialised=false;
  SetProcessSubType(fSurfaceReflection);
  
  theStatus = UndefinedSurf;
  material1 = nullptr;
  material2 = nullptr;
  
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(); 
  theParticleMomentum = 0.;
  
  flag_franchissement_surface = false;
  flag_normal = false;
  flag_reflexion = false;
  teleportToDo = teleportDone = false;

  ekint = thetat = thetaft = energyThreshold = crossingProbability = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecSurface::~G4MicroElecSurface()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool 
G4MicroElecSurface::IsApplicable(const G4ParticleDefinition& aParticleType) 
{ 
  return ( aParticleType.GetPDGEncoding() == 11 ); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecSurface::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (isInitialised) { return; }
  
  G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  G4cout << "G4MicroElecSurface::Initialise: Ncouples= " 
         << numOfCouples << G4endl;
  
  for (G4int i = 0; i < numOfCouples; ++i) {
    const G4Material* material =
      theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
    
    G4cout << "G4Surface, Material " << i + 1 << " / " << numOfCouples << " : " << material->GetName() << G4endl;
    if (material->GetName() == "Vacuum") { tableWF[material->GetName()] = 0; continue; }
    G4String mat = material->GetName();
    G4MicroElecMaterialStructure str = G4MicroElecMaterialStructure(mat);
    tableWF[mat] = str.GetWorkFunction();
  }
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4MicroElecSurface::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  theStatus = UndefinedSurf;
  
  //Definition of the parameters for the particle
  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeVelocity(aTrack.GetVelocity());
  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  
  material1 = pPreStepPoint  -> GetMaterial();
  material2 = pPostStepPoint -> GetMaterial();
  
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  
  theParticleMomentum = aParticle->GetTotalMomentum();
  previousMomentum = oldMomentum;
  oldMomentum = aParticle->GetMomentumDirection(); 
  	
  //First case: not a boundary
  if (pPostStepPoint->GetStepStatus() != fGeomBoundary || 
      pPostStepPoint->GetPhysicalVolume() == pPreStepPoint->GetPhysicalVolume())
    {     
      theStatus = NotAtBoundarySurf;
      flag_franchissement_surface = false;
      flag_reflexion = false;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      
    }
  theStatus = UndefinedSurf;
  
  //Third case: same material  
  if (material1 == material2)
    {      
      theStatus = SameMaterialSurf;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    
    }
  if (verboseLevel > 0)
    {
      G4cout << G4endl << " Electron at Boundary! " << G4endl;
      G4VPhysicalVolume* thePrePV = pPreStepPoint->GetPhysicalVolume();
      G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume();
      if (thePrePV)  G4cout << " thePrePV:  " << thePrePV->GetName() << G4endl;
      if (thePostPV) G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
      G4cout << " Old Momentum Direction: " << oldMomentum << G4endl;
    }
  
  //Definition of the parameters for the surface
  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
  
  G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();
  
  G4bool valid;
  
  theGlobalNormal = theNavigator->GetGlobalExitNormal(theGlobalPoint, &valid);
  //	G4cout << "Global exit normal = " << theGlobalNormal << " valid = " << valid << G4endl;
  
  if (valid)
    {
      theGlobalNormal = -theGlobalNormal;
    }
  else
    {
      G4ExceptionDescription ed;
      ed << " G4MicroElecSurface/PostStepDoIt(): "
	 << " The Navigator reports that it returned an invalid normal.\n"
	 << "PV: " << pPreStepPoint->GetPhysicalVolume()->GetName()
	 << " TrackID= " << aTrack.GetTrackID()
	 << " Ekin(MeV)= " << aTrack.GetKineticEnergy()
	 << " position: " << theGlobalPoint
	 << " direction: " << oldMomentum
	 << G4endl;
      G4Exception("G4MuElecSurf::PostStepDoIt", "OpBoun01",
		  FatalException, ed,
		  "Invalid Surface Normal - Geometry must return valid surface normal");
      return 0;
    }
  
  //Exception: the particle is not in the right direction
  if (oldMomentum * theGlobalNormal > 0.0)
    {
      theGlobalNormal = -theGlobalNormal;
    }
  
  //Second case: step too small
  //Corrections bug rotation + réflexion
  if (aTrack.GetStepLength()<=kCarTolerance)
    { 
      theStatus = StepTooSmallSurf;       
	
      WorkFunctionTable::iterator postStepWF;
      postStepWF = tableWF.find(pPostStepPoint->GetMaterial()->GetName());
      WorkFunctionTable::iterator preStepWF;
      preStepWF = tableWF.find(pPreStepPoint->GetMaterial()->GetName());
	
      if (postStepWF == tableWF.end()) {
	G4String str = "Material ";
	str += pPostStepPoint->GetMaterial()->GetName() + " not found!";
	G4Exception("G4Surface::G4Surface", "em0002", FatalException, str);
	return 0;
      }
	
      else if (preStepWF == tableWF.end()) {
	G4String str = "Material ";
	str += pPreStepPoint->GetMaterial()->GetName() + " not found!";
	G4Exception("G4Surface::G4Surface", "em0002", FatalException, str);
	return 0;
      }
      
      if (pPreStepPoint->GetMaterial() != pPostStepPoint->GetMaterial()) {
	
	flag_franchissement_surface = false;

	if (flag_reflexion == true && flag_normal == true) {
	  aParticleChange.ProposeMomentumDirection(-Reflexion(aStep.GetPostStepPoint()));
	  flag_reflexion = false;
	  flag_normal = false;
	}
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

  flag_normal = (theGlobalNormal.x() == 0.0 && theGlobalNormal.y() == 0.0);
  
  G4LogicalSurface* Surface = nullptr;
  
  Surface = G4LogicalBorderSurface::GetSurface
    (pPreStepPoint ->GetPhysicalVolume(),
     pPostStepPoint->GetPhysicalVolume());
  
  if (Surface == nullptr)
    {
      G4bool enteredDaughter=(pPostStepPoint->GetPhysicalVolume()
			      ->GetMotherLogical() ==
			      pPreStepPoint->GetPhysicalVolume()
			      ->GetLogicalVolume());
      if(enteredDaughter)
	{
	  Surface = G4LogicalSkinSurface::GetSurface
	    (pPostStepPoint->GetPhysicalVolume()->
	     GetLogicalVolume());
	  
	  if(Surface == nullptr)
	    Surface = G4LogicalSkinSurface::GetSurface
	      (pPreStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	}
      else 
	{
	  Surface = G4LogicalSkinSurface::GetSurface
	    (pPreStepPoint->GetPhysicalVolume()->
	     GetLogicalVolume());
	  
	  if(Surface == nullptr)
	    Surface = G4LogicalSkinSurface::GetSurface
	      (pPostStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	}
    }
    
  G4VPhysicalVolume* thePrePV = pPreStepPoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume();
  
  if (thePostPV)
    {
      WorkFunctionTable::iterator postStepWF;
      postStepWF = tableWF.find(thePostPV->GetLogicalVolume()->GetMaterial()->GetName());
      WorkFunctionTable::iterator preStepWF;
      preStepWF = tableWF.find(thePrePV->GetLogicalVolume()->GetMaterial()->GetName());
      
      if (postStepWF == tableWF.end()) {
	G4String str = "Material ";
	str += thePostPV->GetLogicalVolume()->GetMaterial()->GetName() + " not found!";
	G4Exception("G4Surface::G4Surface", "em0002", FatalException, str);
	return 0;
      }
      else if (preStepWF == tableWF.end()) {
	G4String str = "Material ";
	str += thePrePV->GetLogicalVolume()->GetMaterial()->GetName() + " not found!";
	G4Exception("G4Surface::G4Surface", "em0002", FatalException, str);
	return 0;
      }
      else
	{
	  G4double thresholdNew = postStepWF->second;
	  G4double thresholdOld = preStepWF->second;
	  energyThreshold = thresholdNew - thresholdOld;
	}
    }
  
  ekint = pPreStepPoint->GetKineticEnergy();
  thetat= GetIncidentAngle(); //angle d'incidence
  G4double ekinNormalt=ekint*std::cos(thetat)*std::cos(thetat); 
  
  G4double atet = std::sqrt(ekint/(ekint+energyThreshold))*std::sin(thetat);
  thetaft = (atet > 1.0) ? pi*0.5 : std::asin(atet);//Angle de réfraction
  
  G4double aleat=G4UniformRand(); 	
  
  const G4double waveVectort=std::sqrt(2*9.1093826E-31*1.602176487E-19)/(6.6260755E-34/(2.0*pi));
  
  //Parameter for an exponential barrier of potential (Thèse P68)
  const G4double at=0.5E-10;
  
  //G4double modif already declared in .hh
  crossingProbability=0;
  
  G4double kft=waveVectort*std::sqrt(ekint+energyThreshold)*std::cos(thetaft);
  G4double kit=waveVectort*std::sqrt(ekinNormalt); 
  
  G4double yy = std::sinh(pi*at*(kit-kft))/std::sinh(pi*at*(kit+kft));
  crossingProbability = 1 - yy*yy;
  
  //First case: the electron crosses the surface
  if((aleat<=crossingProbability)&&(ekint>std::abs(energyThreshold)))
    {
      if (pPreStepPoint->GetMaterial() != pPostStepPoint->GetMaterial()) {
	flag_franchissement_surface = true;
      }
      
      thetaft=std::abs(thetaft-thetat);
      
      G4ThreeVector zVerst = aStep.GetPostStepPoint()->GetMomentumDirection();
      G4ThreeVector xVerst = zVerst.orthogonal();
      G4ThreeVector yVerst = zVerst.cross(xVerst);

      G4double cost  = std::cos(thetaft); 
      G4double xDirt = std::sqrt(1. - cost*cost);
      G4double yDirt = xDirt;
      
      G4ThreeVector zPrimeVerst = xDirt*xVerst + yDirt*yVerst + cost*zVerst;
		
      aParticleChange.ProposeMomentumDirection(zPrimeVerst.unit());
    }
  else if ((aleat > crossingProbability) && (ekint>std::abs(energyThreshold)))
    {
      flag_reflexion = true;
      if (flag_normal) { aParticleChange.ProposeMomentumDirection(-oldMomentum.unit()); }
      else { aParticleChange.ProposeMomentumDirection(Reflexion(aStep.GetPostStepPoint())); }
      
    }
  else {
    
    if (flag_normal) { aParticleChange.ProposeMomentumDirection(-oldMomentum.unit()); }
    else { aParticleChange.ProposeMomentumDirection(Reflexion(aStep.GetPostStepPoint())); }
    flag_reflexion = true;
  }
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecSurface::GetMeanFreePath(const G4Track&, G4double,
                                                 G4ForceCondition* condition)
{ 
  *condition = Forced;
  return  DBL_MAX; 	 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecSurface::GetIncidentAngle() 
{
  theFacetNormal=theGlobalNormal; 
  
  G4double PdotN = oldMomentum * theFacetNormal;
  G4double magP= oldMomentum.mag();
  G4double magN= theFacetNormal.mag();
  
  G4double incidentangle = pi - std::acos(PdotN/(magP*magN));
  
  return incidentangle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4MicroElecSurface::Reflexion(const G4StepPoint* PostStepPoint)
{
  //Normale
  G4double Nx = theGlobalNormal.x();
  G4double Ny = theGlobalNormal.y();
  G4double Nz = theGlobalNormal.z();
  
  //PostStepPoint
  G4double PSx = PostStepPoint->GetPosition().x();
  G4double PSy = PostStepPoint->GetPosition().y();
  G4double PSz = PostStepPoint->GetPosition().z();
  
  //P(alpha,beta,gamma) - PostStep avec translation momentum
  G4double alpha = PSx + oldMomentum.x();
  G4double beta = PSy + oldMomentum.y();
  G4double gamma = PSz + oldMomentum.z();
  G4double r = theGlobalNormal.mag();
  G4double x, y, z, d, A, B, PM2x, PM2y, PM2z;
  d = -(Nx*PSx + Ny*PSy + Nz*PSz);
  
  if (Ny == 0 && Nx == 0) {
    gamma = -gamma;
  }  
  else {
    if (Ny == 0) {
      
      A = (Nz*Nz*alpha) + (Nx*Nx*PSx) + (Nx*Nz*(PSz - gamma));
      B = r*r;
      
      //M(x,y,z) - Projection de P sur la surface
      x = A / B;
      y = beta;
      z = (x - alpha)*(Nz / Nx) + gamma;
    }
    
    else {
      
      A = (r*r) / Ny;
      B = (beta / Ny)*(Nx*Nx + Nz*Nz) - (Nx*alpha + Nz*gamma + d);
      
      //M(x,y,z) - Projection de P sur la surface
      y = B / A;
      x = (y - beta)*(Nx / Ny) + alpha;
      z = (y - beta)*(Nz / Ny) + gamma;
    }
    
    //Vecteur 2*PM
    PM2x = 2 * (x - alpha);	 PM2y = 2 * (y - beta);	 PM2z = 2 * (z - gamma);
    
    //Nouveau point P
    alpha += PM2x; beta += PM2y; gamma += PM2z;
    
  }
  
  G4ThreeVector newMomentum = G4ThreeVector(alpha-PSx,beta-PSy,gamma-PSz);  
  return newMomentum.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecSurfaceStatus G4MicroElecSurface::GetStatus() const 
{ 
  return theStatus; 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  
void G4MicroElecSurface::SetFlagFranchissement() 
{ 
  flag_franchissement_surface = false; 
} 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

