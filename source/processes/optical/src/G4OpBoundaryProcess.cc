////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.cc
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    1998-05-25 - Correct parallel component of polarization
//                           (thanks to: Stefano Magni + Giovanni Pieri)
//              1998-05-28 - NULL Rindex pointer before reuse
//                           (thanks to: Stefano Magni)
//              1998-06-11 - delete *sint1 in oblique reflection
//                           (thanks to: Giovanni Pieri)
//              1998-06-19 - move from GetLocalExitNormal() to the new 
//                           method: GetLocalExitNormal(&valid) to get
//                           the surface normal in all cases
//              1998-11-07 - NULL OpticalSurface pointer before use
//                           comparison not sharp for: abs(cost1) < 1.0
//                           remove sin1, sin2 in lines 556,567
//                           (thanks to Stefano Magni)
//
// Author:      Peter Gumplinger
// 		adopted from work by Werner Keil - April 2/96
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpBoundaryProcess.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4OpBoundaryProcess::operator=(const G4OpBoundaryProcess &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4OpBoundaryProcess::G4OpBoundaryProcess(const G4String& processName)
             : G4VDiscreteProcess(processName)
{
        if ( verboseLevel > 0) {
           G4cout << GetProcessName() << " is created " << endl;
        }

	theStatus = Undefined;
	theModel = glisur;
	theFinish = polished;
}

// G4OpBoundaryProcess::G4OpBoundaryProcess(const G4OpBoundaryProcess &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4OpBoundaryProcess::~G4OpBoundaryProcess(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// ------------
//
G4VParticleChange*
G4OpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        aParticleChange.Initialize(aTrack);

        G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
        G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

        if (pPostStepPoint->GetStepStatus() != fGeomBoundary)
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	if (aTrack.GetStepLength()<=kCarTolerance/2)
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	Material1 = pPreStepPoint ->GetPhysicalVolume()->
				    GetLogicalVolume()->GetMaterial();
	Material2 = pPostStepPoint->GetPhysicalVolume()->
				    GetLogicalVolume()->GetMaterial();

	if (Material1 == Material2)
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

	thePhotonMomentum = aParticle->GetTotalMomentum();
        OldMomentum       = aParticle->GetMomentumDirection();
	OldPolarization   = aParticle->GetPolarization();

        if ( verboseLevel > 0 ) {
		G4cout << " Photon at Boundary! " << endl;
		G4cout << " Old Momentum Direction: " << OldMomentum     << endl;
		G4cout << " Old Polarization:       " << OldPolarization << endl;
        }

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
        G4MaterialPropertyVector* Rindex;

	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	}
	else {
		aParticleChange.SetStatusChange(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        if (Rindex) {
		Rindex1 = Rindex->GetProperty(thePhotonMomentum);
	}
	else {
		aParticleChange.SetStatusChange(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        Rindex = NULL;
        OpticalSurface = NULL;

	aMaterialPropertiesTable = Material2->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) 
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

        G4LogicalSurface* Surface = G4LogicalBorderSurface::GetSurface
				    (pPreStepPoint ->GetPhysicalVolume(),
				     pPostStepPoint->GetPhysicalVolume());

        if (Surface == NULL) Surface = G4LogicalSkinSurface::GetSurface
				       (pPreStepPoint->GetPhysicalVolume()->
						GetLogicalVolume());

	if (Surface != NULL) OpticalSurface = Surface->GetOpticalSurface();

	theModel = glisur;
	theFinish = polished;

	G4OpticalSurfaceType type;
	if (Rindex) {
	   type = dielectric_dielectric;
//	   if (OpticalSurface) type = OpticalSurface->GetType();
	   Rindex2 = Rindex->GetProperty(thePhotonMomentum);
	}
	else if (OpticalSurface) {
	   type = OpticalSurface->GetType();
	}
	else {
	   aParticleChange.SetStatusChange(fStopAndKill);
	   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	if (OpticalSurface) {
	   theModel  = OpticalSurface->GetModel();
	   theFinish = OpticalSurface->GetFinish();

	   aMaterialPropertiesTable = OpticalSurface->
					GetMaterialPropertiesTable();

           if (aMaterialPropertiesTable) {
	      G4MaterialPropertyVector* PropertyPointer;

	      if(!Rindex) {
	         PropertyPointer = 
	         aMaterialPropertiesTable->GetProperty("RINDEX");
	         if (PropertyPointer) Rindex2 = 
			 PropertyPointer->GetProperty(thePhotonMomentum);
	      }

	      PropertyPointer = 
	      aMaterialPropertiesTable->GetProperty("REFLECTIVITY");
	      if (PropertyPointer) theReflectivity = 
		      PropertyPointer->GetProperty(thePhotonMomentum);

	      PropertyPointer = 
	      aMaterialPropertiesTable->GetProperty("EFFICIENCY");
	      if (PropertyPointer) theEfficiency = 
		      PropertyPointer->GetProperty(thePhotonMomentum);

	      if ( theModel == unified ) {
	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("SPECULARLOBECONSTANT");
	        if (PropertyPointer) prob_sl = 
			 PropertyPointer->GetProperty(thePhotonMomentum);

	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("SPECULARSPIKECONSTANT");
	        if (PropertyPointer) prob_ss = 
			 PropertyPointer->GetProperty(thePhotonMomentum);

	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("BACKSCATTERCONSTANT");
	        if (PropertyPointer) prob_bs = 
			 PropertyPointer->GetProperty(thePhotonMomentum);
	      }
	   }
	}

	G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();

        G4Navigator* theNavigator =
                     G4TransportationManager::GetTransportationManager()->
                                              GetNavigatorForTracking();

	G4ThreeVector theLocalPoint = theNavigator->
        			      GetGlobalToLocalTransform().
				      TransformPoint(theGlobalPoint);

	G4ThreeVector theLocalNormal;	// Normal points back into volume

	G4bool valid;
	theLocalNormal = theNavigator->GetLocalExitNormal(&valid);

	if (valid) {
	  theLocalNormal = -theLocalNormal;
	}
	else {
	  G4cerr << " G4OpBoundaryProcess/PostStepDoIt(): "
	       << " The Navigator reports that it returned an invalid normal" 
	       << endl;
	}

	theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
	                                TransformAxis(theLocalNormal);

	theStatus = Undefined;

	if (type == dielectric_metal) {

	  DielectricMetal();

	}
	else if (type == dielectric_dielectric) {

	  if ( theFinish == polishedfrontpainted ||
	       theFinish == groundfrontpainted ) {
	          if( !G4BooleanRand(theReflectivity) ) {
		    DoAbsorption();
		  }
	          else {
		    if ( theFinish == groundfrontpainted )
					theStatus = LambertianReflection;
		    DoReflection();
		  }
	  }
	  else {
		  DielectricDielectric();
	  }
	}
	else {

	  G4cout << " Error: G4BoundaryProcess: illegal boundary type " << endl;
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	}

        NewMomentum = NewMomentum.unit();
        NewPolarization = NewPolarization.unit();

        if ( verboseLevel > 0) {
		G4cout << " New Momentum Direction: " << NewMomentum     << endl;
		G4cout << " New Polarization:       " << NewPolarization << endl;
		if ( theStatus == Undefined )
			G4cout << " *** Undefined *** " << endl;
		if ( theStatus == FresnelRefraction )
			G4cout << " *** FresnelRefraction *** " << endl;
		if ( theStatus == FresnelReflection )
			G4cout << " *** FresnelReflection *** " << endl;
		if ( theStatus == TotalInternalReflection )
			G4cout << " *** TotalInternalReflection *** " << endl;
		if ( theStatus == LambertianReflection )
			G4cout << " *** LambertianReflection *** " << endl;
		if ( theStatus == LobeReflection ) 
			G4cout << " *** LobeReflection *** " << endl;
		if ( theStatus == SpikeReflection )
			G4cout << " *** SpikeReflection *** " << endl;
		if ( theStatus == BackScattering )
			G4cout << " *** BackScattering *** " << endl;
		if ( theStatus == Absorption )
			G4cout << " *** Absorption *** " << endl;
		if ( theStatus == Detection )
			G4cout << " *** Detection *** " << endl;
        }

	aParticleChange.SetMomentumChange(NewMomentum);
	aParticleChange.SetPolarizationChange(NewPolarization);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4ThreeVector 
G4OpBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
        G4ThreeVector FacetNormal;

	if (theModel == unified) {

	/* This function code alpha to a random value taken from the
           distribution p(alpha) = g(alpha; 0, sigma_alpha)*sin(alpha), 
           for alpha > 0 and alpha < 90, where g(alpha; 0, sigma_alpha) 
           is a gaussian distribution with mean 0 and standard deviation 
           sigma_alpha.  */

	   G4double alpha;

	   G4double sigma_alpha = 0.0;
	   if (OpticalSurface) sigma_alpha = OpticalSurface->GetSigmaAlpha();

	   G4double f_max = min(1.0,4.*sigma_alpha);

	   do {
	      do {
	         alpha = RandGauss::shoot(0.0,sigma_alpha);
	      } while (G4UniformRand()*f_max > sin(alpha) || alpha >= halfpi );

	      G4double phi = G4UniformRand()*twopi;

	      G4double SinAlpha = sin(alpha);
	      G4double CosAlpha = cos(alpha);
              G4double SinPhi = sin(phi);
              G4double CosPhi = cos(phi);

              G4double unit_x = SinAlpha * CosPhi;
              G4double unit_y = SinAlpha * SinPhi;
              G4double unit_z = CosAlpha;

	      FacetNormal.setX(unit_x);
	      FacetNormal.setY(unit_y);
	      FacetNormal.setZ(unit_z);

	      G4ThreeVector tmpNormal = Normal;

	      FacetNormal.rotateUz(tmpNormal);
	   } while (Momentum * FacetNormal >= 0.0);
	}
	else {

	   G4double  polish = 1.0;
	   if (OpticalSurface) polish = OpticalSurface->GetPolish();

           if (polish < 1.0) {
              do {
                 G4ThreeVector smear;
                 do {
                    smear.setX(2.*G4UniformRand()-1.0);
                    smear.setY(2.*G4UniformRand()-1.0);
                    smear.setZ(2.*G4UniformRand()-1.0);
                 } while (smear.mag()>1.0);
                 smear = (1.-polish) * smear;
                 FacetNormal = Normal + smear;
              } while (Momentum * FacetNormal >= 0.0);
              FacetNormal = FacetNormal.unit();
	   }
           else {
              FacetNormal = Normal;
           }
	}
	return FacetNormal;
}

void G4OpBoundaryProcess::DielectricMetal()
{
	do {
           if( !G4BooleanRand(theReflectivity) ) {

	     DoAbsorption();

           }
           else {

	     DoReflection();

	   }

	   OldMomentum = NewMomentum;
	   OldPolarization = NewPolarization;

	} while (NewMomentum * theGlobalNormal < 0.0);
}

void G4OpBoundaryProcess::DielectricDielectric()
{
	G4bool Inside = false;
	G4bool Swap = false;

	leap:

        G4bool Through = false;
	G4bool Done = false;

	do {

	   if (Through) {
	      Swap = !Swap;
	      Through = false;
	      theGlobalNormal = -theGlobalNormal;
	      G4Swap(Material1,Material2);
	      G4Swap(&Rindex1,&Rindex2);
	   }

	   if ( theFinish == ground || theFinish == groundbackpainted ) {
		theFacetNormal = 
			     GetFacetNormal(OldMomentum,theGlobalNormal);
	   }
	   else {
		theFacetNormal = theGlobalNormal;
	   }

	   G4double PdotN = OldMomentum * theFacetNormal;
	   G4double EdotN = OldPolarization * theFacetNormal;

	   cost1 = - PdotN;
	   if (abs(cost1) < 1.0-kCarTolerance){
	      sint1 = sqrt(1-cost1*cost1);
	      sint2 = sint1*Rindex1/Rindex2;     // *** Snell's Law ***
	   }
	   else {
	      sint1 = 0.0;
	      sint2 = 0.0;
	   }

	   if (sint2 >= 1.0) {

	      // Simulate total internal reflection

	      if (Swap) Swap = !Swap;

              theStatus = TotalInternalReflection;

	      if ( theModel == unified && theFinish != polished )
						     ChooseReflection();

	      if ( theStatus == LambertianReflection ) {
		 DoReflection();
	      }
	      else if ( theStatus == BackScattering ) {
		 NewMomentum = -OldMomentum;
		 NewPolarization = -OldPolarization;
	      }
	      else {

                 PdotN = OldMomentum * theFacetNormal;
		 NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
		 EdotN = OldPolarization * theFacetNormal;
		 NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;

	      }
	   }
	   else if (sint2 < 1.0) {

	      // Calculate amplitude for transmission (Q = P x N)

	      if (cost1 > 0.0) {
	         cost2 =  sqrt(1-sint2*sint2);
	      }
	      else {
	         cost2 = -sqrt(1-sint2*sint2);
	      }

	      G4ThreeVector A_trans, Atrans, E1pp, E1pl;
	      G4double E1_perp, E1_parl;

	      if (sint1 > 0.0) {
	         A_trans = OldMomentum.cross(theFacetNormal);
                 Atrans  = A_trans.unit();
	         E1_perp = OldPolarization * Atrans;
                 E1pp    = E1_perp * Atrans;
                 E1pl    = OldPolarization - E1pp;
                 E1_parl = E1pl.mag();
	      }
	      else {
	         A_trans  = OldPolarization;
	         // Here we Follow Jackson's conventions and we set the
	         // parallel component = 1 in case of a ray perpendicular
	         // to the surface
	         E1_perp  = 0.0;
	         E1_parl  = 1.0;
	      }

	      G4double E2_perp, E2_parl, E2_total, TransCoeff;

	      if (cost1 != 0.0) {
	         G4double s1 = Rindex1*cost1;
	         E2_perp = 2.*s1*E1_perp/(Rindex1*cost1+Rindex2*cost2);
	         E2_parl = 2.*s1*E1_parl/(Rindex2*cost1+Rindex1*cost2);
	         E2_total = E2_perp*E2_perp + E2_parl*E2_parl;
	         G4double s2 = Rindex2*cost2*E2_total;
	         TransCoeff = s2/s1;
	      }
	      else {
	         TransCoeff = 0.0;
	      }

	      G4ThreeVector Refracted, Deflected;
	      G4double E2_abs, C_parl, C_perp;

	      if ( !G4BooleanRand(TransCoeff) ) {

	         // Simulate reflection

                 if (Swap) Swap = !Swap;

		 theStatus = FresnelReflection;
		 if ( theModel == unified && theFinish != polished )
						     ChooseReflection();

		 if ( theStatus == LambertianReflection ) {
		    DoReflection();
		 }
		 else if ( theStatus == BackScattering ) {
		    NewMomentum = -OldMomentum;
		    NewPolarization = -OldPolarization;
		 }
		 else {

                    PdotN = OldMomentum * theFacetNormal;
	            NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

	            if (sint1 > 0.0) {              // incident ray oblique

		       E2_parl   = Rindex2*E2_parl/Rindex1 - E1_parl;
		       E2_perp   = E2_perp - E1_perp;
		       E2_total  = E2_perp*E2_perp + E2_parl*E2_parl;
		       Refracted = theFacetNormal + PdotN * NewMomentum;
		       E2_abs    = sqrt(E2_total);
		       C_parl    = E2_parl/E2_abs;
		       C_perp    = E2_perp/E2_abs;

		       NewPolarization = C_parl*Refracted - C_perp*A_trans;

	             }

	             else if (Rindex2 > Rindex1) { // incident ray perpendicular

		       NewPolarization = - OldPolarization;

	             }
		 }
	      }
	      else { // photon gets transmitted

	         // Simulate transmission/refraction

		 Inside = !Inside;
		 Through = true;
		 theStatus = FresnelRefraction;

	         if (sint1 > 0.0) {      // incident ray oblique

		    G4double alpha = cost1 - cost2*(Rindex2/Rindex1);
		    Deflected  = OldMomentum + alpha*theFacetNormal;
		    NewMomentum = Deflected.unit();
		    PdotN = -cost2;
		    Refracted  = theFacetNormal - PdotN*NewMomentum;
		    E2_abs     = sqrt(E2_total);
		    C_parl     = E2_parl/E2_abs;
		    C_perp     = E2_perp/E2_abs;
		    NewPolarization = C_parl*Refracted + C_perp*A_trans;

	         }
	         else {                  // incident ray perpendicular

		    NewMomentum = OldMomentum;
		    NewPolarization = OldPolarization;

	         }
	      }
	   }

	   OldMomentum = NewMomentum;
	   OldPolarization = NewPolarization;

	   if (theStatus == FresnelRefraction) {
	      Done = (NewMomentum * theGlobalNormal <= 0.0);
	   } 
	   else {
	      Done = (NewMomentum * theGlobalNormal >= 0.0);
	   }

	} while (!Done);

	if (Inside && !Swap) {
          if( theFinish == polishedbackpainted ||
              theFinish == groundbackpainted ) {
	      if( !G4BooleanRand(theReflectivity) ) {
		DoAbsorption();
              }
	      else {
		if (theStatus != FresnelRefraction ) {
		   theGlobalNormal = -theGlobalNormal;
	        }
	        else {
		   Swap = !Swap;
		   G4Swap(Material1,Material2);
		   G4Swap(&Rindex1,&Rindex2);
	        }
		if ( theFinish == groundbackpainted )
					theStatus = LambertianReflection;

	        DoReflection();

	        theGlobalNormal = -theGlobalNormal;
		OldMomentum = NewMomentum;

	        goto leap;
	      }
	  }
	}
}

// GetMeanFreePath
// ---------------
//
G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track& ,
                                              G4double ,
                                              G4ForceCondition* condition)
{
	*condition = Forced;

	return DBL_MAX;
}
