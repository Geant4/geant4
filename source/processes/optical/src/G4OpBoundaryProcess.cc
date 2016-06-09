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
//                           comparison not sharp for: std::abs(cost1) < 1.0
//                           remove sin1, sin2 in lines 556,567
//                           (thanks to Stefano Magni)
//              1999-10-10 - Accommodate changes done in DoAbsorption by
//                           changing logic in DielectricMetal
//              2001-10-18 - avoid Linux (gcc-2.95.2) warning about variables
//                           might be used uninitialized in this function
//                           moved E2_perp, E2_parl and E2_total out of 'if'
//              2003-11-27 - Modified line 168-9 to reflect changes made to
//                           G4OpticalSurface class ( by Fan Lei)
//              2004-02-02 - Set theStatus = Undefined at start of DoIt
//              2005-07-28 - add G4ProcessType to constructor
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

G4OpBoundaryProcess::G4OpBoundaryProcess(const G4String& processName,
                                               G4ProcessType type)
             : G4VDiscreteProcess(processName, type)
{
        if ( verboseLevel > 0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

	theStatus = Undefined;
	theModel = glisur;
	theFinish = polished;
        theReflectivity = 1.;
        theEfficiency   = 0.;
                                                                                
        prob_sl = 0.;
        prob_ss = 0.;
        prob_bs = 0.;

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
        theStatus = Undefined;

        aParticleChange.Initialize(aTrack);

        G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
        G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

        if (pPostStepPoint->GetStepStatus() != fGeomBoundary){
	        theStatus = NotAtBoundary;
	        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	if (aTrack.GetStepLength()<=kCarTolerance/2){
	        theStatus = StepTooSmall;
	        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	Material1 = pPreStepPoint  -> GetMaterial();
	Material2 = pPostStepPoint -> GetMaterial();

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

	thePhotonMomentum = aParticle->GetTotalMomentum();
        OldMomentum       = aParticle->GetMomentumDirection();
	OldPolarization   = aParticle->GetPolarization();

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
        G4MaterialPropertyVector* Rindex;

	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	}
	else {
	        theStatus = NoRINDEX;
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        if (Rindex) {
		Rindex1 = Rindex->GetProperty(thePhotonMomentum);
	}
	else {
	        theStatus = NoRINDEX;
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        theModel = glisur;
        theFinish = polished;

        G4SurfaceType type = dielectric_dielectric;

        Rindex = NULL;
        OpticalSurface = NULL;

        G4LogicalSurface* Surface = G4LogicalBorderSurface::GetSurface
				    (pPreStepPoint ->GetPhysicalVolume(),
				     pPostStepPoint->GetPhysicalVolume());

        if (Surface == NULL){
	  G4bool enteredDaughter=(pPostStepPoint->GetPhysicalVolume()
				  ->GetMotherLogical() == 
				  pPreStepPoint->GetPhysicalVolume()
				  ->GetLogicalVolume());
	  if(enteredDaughter){
	    Surface = G4LogicalSkinSurface::GetSurface
	      (pPostStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	    if(Surface == NULL)
	      Surface = G4LogicalSkinSurface::GetSurface
	      (pPreStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	  }
	  else{
	    Surface = G4LogicalSkinSurface::GetSurface
	      (pPreStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	    if(Surface == NULL)
	      Surface = G4LogicalSkinSurface::GetSurface
	      (pPostStepPoint->GetPhysicalVolume()->
	       GetLogicalVolume());
	  }
	}

	//	if (Surface) OpticalSurface = dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty());
	if (Surface) OpticalSurface = (G4OpticalSurface*) Surface->GetSurfaceProperty();

	if (OpticalSurface) {

           type      = OpticalSurface->GetType();
	   theModel  = OpticalSurface->GetModel();
	   theFinish = OpticalSurface->GetFinish();

	   aMaterialPropertiesTable = OpticalSurface->
					GetMaterialPropertiesTable();

           if (aMaterialPropertiesTable) {

              if (theFinish == polishedbackpainted ||
                  theFinish == groundbackpainted ) {
                  Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	          if (Rindex) {
                     Rindex2 = Rindex->GetProperty(thePhotonMomentum);
                  }
                  else {
		     theStatus = NoRINDEX;
                     aParticleChange.ProposeTrackStatus(fStopAndKill);
                     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
                  }
              }

              G4MaterialPropertyVector* PropertyPointer;

	      PropertyPointer = 
	      aMaterialPropertiesTable->GetProperty("REFLECTIVITY");
	      if (PropertyPointer) { 
                      theReflectivity =
		      PropertyPointer->GetProperty(thePhotonMomentum);
              } else {
                      theReflectivity = 1.0;
              }

	      PropertyPointer = 
	      aMaterialPropertiesTable->GetProperty("EFFICIENCY");
	      if (PropertyPointer) {
                      theEfficiency =
		      PropertyPointer->GetProperty(thePhotonMomentum);
              } else {
                      theEfficiency = 0.0;
              }

	      if ( theModel == unified ) {
	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("SPECULARLOBECONSTANT");
	        if (PropertyPointer) {
                         prob_sl =
			 PropertyPointer->GetProperty(thePhotonMomentum);
                } else {
                         prob_sl = 0.0;
                }

	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("SPECULARSPIKECONSTANT");
	        if (PropertyPointer) {
                         prob_ss =
			 PropertyPointer->GetProperty(thePhotonMomentum);
                } else {
                         prob_ss = 0.0;
                }

	        PropertyPointer = 
		aMaterialPropertiesTable->GetProperty("BACKSCATTERCONSTANT");
	        if (PropertyPointer) {
                         prob_bs =
			 PropertyPointer->GetProperty(thePhotonMomentum);
                } else {
                         prob_bs = 0.0;
                }
	      }
	   }
           else if (theFinish == polishedbackpainted ||
                    theFinish == groundbackpainted ) {
                      aParticleChange.ProposeTrackStatus(fStopAndKill);
                      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
           }
        }

        if (type == dielectric_dielectric ) {
           if (theFinish == polished || theFinish == ground ) {

	      if (Material1 == Material2){
		 theStatus = SameMaterial;
		 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	      }
              aMaterialPropertiesTable = 
                     Material2->GetMaterialPropertiesTable();
              if (aMaterialPropertiesTable)
                 Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
              if (Rindex) {
                 Rindex2 = Rindex->GetProperty(thePhotonMomentum);
              }
              else {
		 theStatus = NoRINDEX;
                 aParticleChange.ProposeTrackStatus(fStopAndKill);
                 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	      }
           }
        }

        if ( verboseLevel > 0 ) {
                G4cout << " Photon at Boundary! " << G4endl;
                G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
                G4cout << " Old Polarization:       " << OldPolarization << G4endl;
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
	       << G4endl;
	}

	theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
	                                TransformAxis(theLocalNormal);
        if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4DEBUG_OPTICAL
           G4cerr << " G4OpBoundaryProcess/PostStepDoIt(): "
                  << " theGlobalNormal points the wrong direction "
                  << G4endl;
#endif
           theGlobalNormal = -theGlobalNormal;
        }
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

	  G4cerr << " Error: G4BoundaryProcess: illegal boundary type " << G4endl;
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	}

        NewMomentum = NewMomentum.unit();
        NewPolarization = NewPolarization.unit();

        if ( verboseLevel > 0) {
		G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
		G4cout << " New Polarization:       " << NewPolarization << G4endl;
		if ( theStatus == Undefined )
			G4cout << " *** Undefined *** " << G4endl;
		if ( theStatus == FresnelRefraction )
			G4cout << " *** FresnelRefraction *** " << G4endl;
		if ( theStatus == FresnelReflection )
			G4cout << " *** FresnelReflection *** " << G4endl;
		if ( theStatus == TotalInternalReflection )
			G4cout << " *** TotalInternalReflection *** " << G4endl;
		if ( theStatus == LambertianReflection )
			G4cout << " *** LambertianReflection *** " << G4endl;
		if ( theStatus == LobeReflection ) 
			G4cout << " *** LobeReflection *** " << G4endl;
		if ( theStatus == SpikeReflection )
			G4cout << " *** SpikeReflection *** " << G4endl;
		if ( theStatus == BackScattering )
			G4cout << " *** BackScattering *** " << G4endl;
		if ( theStatus == Absorption )
			G4cout << " *** Absorption *** " << G4endl;
		if ( theStatus == Detection )
			G4cout << " *** Detection *** " << G4endl;
                if ( theStatus == NotAtBoundary )
                        G4cout << " *** NotAtBoundary *** " << G4endl;
                if ( theStatus == SameMaterial )
                        G4cout << " *** SameMaterial *** " << G4endl;
                if ( theStatus == StepTooSmall )
                        G4cout << " *** StepTooSmall *** " << G4endl;
                if ( theStatus == NoRINDEX )
                        G4cout << " *** NoRINDEX *** " << G4endl;
        }

	aParticleChange.ProposeMomentumDirection(NewMomentum);
	aParticleChange.ProposePolarization(NewPolarization);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4ThreeVector 
G4OpBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
        G4ThreeVector FacetNormal;

	if (theModel == unified) {

	/* This function code alpha to a random value taken from the
           distribution p(alpha) = g(alpha; 0, sigma_alpha)*std::sin(alpha), 
           for alpha > 0 and alpha < 90, where g(alpha; 0, sigma_alpha) 
           is a gaussian distribution with mean 0 and standard deviation 
           sigma_alpha.  */

	   G4double alpha;

	   G4double sigma_alpha = 0.0;
	   if (OpticalSurface) sigma_alpha = OpticalSurface->GetSigmaAlpha();

	   G4double f_max = std::min(1.0,4.*sigma_alpha);

	   do {
	      do {
	         alpha = G4RandGauss::shoot(0.0,sigma_alpha);
	      } while (G4UniformRand()*f_max > std::sin(alpha) || alpha >= halfpi );

	      G4double phi = G4UniformRand()*twopi;

	      G4double SinAlpha = std::sin(alpha);
	      G4double CosAlpha = std::cos(alpha);
              G4double SinPhi = std::sin(phi);
              G4double CosPhi = std::cos(phi);

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
        G4int n = 0;

	do {

           n++;

           if( !G4BooleanRand(theReflectivity) && n == 1 ) {

             DoAbsorption();
             break;

           }
           else {

             if ( theModel == glisur || theFinish == polished ) {

                DoReflection();

             } else {

                if ( n == 1 ) ChooseReflection();
                                                                                
                if ( theStatus == LambertianReflection ) {
                   DoReflection();
                }
                else if ( theStatus == BackScattering ) {
                   NewMomentum = -OldMomentum;
                   NewPolarization = -OldPolarization;
                }
                else {

                   if(theStatus==LobeReflection)theFacetNormal = 
                             GetFacetNormal(OldMomentum,theGlobalNormal);

                   G4double PdotN = OldMomentum * theFacetNormal;
                   NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
                   G4double EdotN = OldPolarization * theFacetNormal;
                   NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;                                                                                
                }

             }

             OldMomentum = NewMomentum;
             OldPolarization = NewPolarization;

	   }

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
	   if (std::abs(cost1) < 1.0-kCarTolerance){
	      sint1 = std::sqrt(1.-cost1*cost1);
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
	         cost2 =  std::sqrt(1.-sint2*sint2);
	      }
	      else {
	         cost2 = -std::sqrt(1.-sint2*sint2);
	      }

	      G4ThreeVector A_trans, A_paral, E1pp, E1pl;
	      G4double E1_perp, E1_parl;

	      if (sint1 > 0.0) {
	         A_trans = OldMomentum.cross(theFacetNormal);
                 A_trans = A_trans.unit();
	         E1_perp = OldPolarization * A_trans;
                 E1pp    = E1_perp * A_trans;
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

              G4double s1 = Rindex1*cost1;
              G4double E2_perp = 2.*s1*E1_perp/(Rindex1*cost1+Rindex2*cost2);
              G4double E2_parl = 2.*s1*E1_parl/(Rindex2*cost1+Rindex1*cost2);
              G4double E2_total = E2_perp*E2_perp + E2_parl*E2_parl;
              G4double s2 = Rindex2*cost2*E2_total;

              G4double TransCoeff;

	      if (cost1 != 0.0) {
	         TransCoeff = s2/s1;
	      }
	      else {
	         TransCoeff = 0.0;
	      }

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

	            if (sint1 > 0.0) {   // incident ray oblique

		       E2_parl   = Rindex2*E2_parl/Rindex1 - E1_parl;
		       E2_perp   = E2_perp - E1_perp;
		       E2_total  = E2_perp*E2_perp + E2_parl*E2_parl;
                       A_paral   = NewMomentum.cross(A_trans);
                       A_paral   = A_paral.unit();
		       E2_abs    = std::sqrt(E2_total);
		       C_parl    = E2_parl/E2_abs;
		       C_perp    = E2_perp/E2_abs;

                       NewPolarization = C_parl*A_paral + C_perp*A_trans;

	            }

	            else {               // incident ray perpendicular

	               if (Rindex2 > Rindex1) {
		          NewPolarization = - OldPolarization;
	               }
	               else {
	                  NewPolarization =   OldPolarization;
	               }

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
		    NewMomentum = OldMomentum + alpha*theFacetNormal;
		    NewMomentum = NewMomentum.unit();
		    PdotN = -cost2;
                    A_paral = NewMomentum.cross(A_trans);
                    A_paral = A_paral.unit();
		    E2_abs     = std::sqrt(E2_total);
		    C_parl     = E2_parl/E2_abs;
		    C_perp     = E2_perp/E2_abs;

                    NewPolarization = C_parl*A_paral + C_perp*A_trans;

	         }
	         else {                  // incident ray perpendicular

		    NewMomentum = OldMomentum;
		    NewPolarization = OldPolarization;

	         }
	      }
	   }

	   OldMomentum = NewMomentum.unit();
	   OldPolarization = NewPolarization.unit();

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

