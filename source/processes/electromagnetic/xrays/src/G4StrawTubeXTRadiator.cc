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
// $Id: G4StrawTubeXTRadiator.cc,v 1.1 2005/04/22 09:44:18 grichine Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include <complex>

#include "G4StrawTubeXTRadiator.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4StrawTubeXTRadiator::G4StrawTubeXTRadiator(G4LogicalVolume *anEnvelope,
					 G4Material* foilMat,G4Material* gasMat, 
                                         G4double a, G4double b, G4Material* mediumMat,
                                         G4bool unishut,
                                         const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,1,processName)
{
  G4cout<<"Straw tube X-ray TR  radiator EM process is called"<<G4endl;

  if( unishut )
  {
    fAlphaPlate = 25;
    fAlphaGas   = 25;
    G4cout<<"straw uniform shuting: "<<"fAlphaPlate = "
          <<fAlphaPlate<<" ; fAlphaGas = "<<fAlphaGas<<G4endl;

  }
  else
  {
    fAlphaPlate = 9;
    fAlphaGas   = 9;
    G4cout<<"straw isotropical shuting: "<<"fAlphaPlate = "
          <<fAlphaPlate<<" ; fAlphaGas = "<<fAlphaGas<<G4endl;


  }
  // index of medium material

  fMatIndex3 = mediumMat->GetIndex();
  G4cout<<"medium material = "<<mediumMat->GetName()<<G4endl;

 // plasma energy squared for plate material

  fSigma3 = fPlasmaCof*mediumMat->GetElectronDensity();
  G4cout<<"medium plasma energy = "<<sqrt(fSigma3)/eV<<" eV"<<G4endl;

// Compute cofs for preparation of linear photo absorption in external medium

  ComputeMediumPhotoAbsCof();




  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable();
}

///////////////////////////////////////////////////////////////////////////

G4StrawTubeXTRadiator::~G4StrawTubeXTRadiator()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Approximation for radiator interference factor for the case of
// straw tube radiator. The plate (window, straw wall) and gas (inside straw) 
// gap thicknesses are  gamma distributed.
// The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zone.

G4double 
G4StrawTubeXTRadiator::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{


  G4double result, L2, L3, M2, M3;
  
  L2 = GetPlateFormationZone(energy,gamma,varAngle);
  L3 = GetGasFormationZone(energy,gamma,varAngle);

  M2 = GetPlateLinearPhotoAbs(energy);
  M3 = GetGasLinearPhotoAbs(energy);


  G4complex C2(1.0 + 0.5*fPlateThick*M2/fAlphaPlate, fPlateThick/L2/fAlphaPlate); 
  G4complex C3(1.0 + 0.5*fGasThick*M3/fAlphaGas, fGasThick/L3/fAlphaGas); 

  G4complex H2 = pow(C2,-fAlphaPlate);  
  G4complex H3 = pow(C3,-fAlphaGas);
  G4complex H  = H2*H3;

  G4complex Z1 = GetMediumComplexFZ(energy,gamma,varAngle);
  G4complex Z2 = GetPlateComplexFZ(energy,gamma,varAngle);
  G4complex Z3 = GetGasComplexFZ(energy,gamma,varAngle);


  G4complex R  =    ( Z1 - Z2 )*( Z1 - Z2 )*( 1. - H2*H ) +
                    ( Z2 - Z3 )*( Z2 - Z3 )*( 1. - H3 )   + 
                 2.*( Z1 - Z2 )*( Z2 - Z3 )*H2*( 1. - H3 ) ;

  result       = 2.0*real(R)*(varAngle*energy/hbarc/hbarc);
  
  return      result;

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for external medium. Omega is energy !!!

G4double G4StrawTubeXTRadiator::GetMediumFormationZone( G4double omega ,
                                                G4double gamma ,
                                                G4double varAngle    ) 
{
  G4double cof, lambda;
  lambda = 1.0/gamma/gamma + varAngle + fSigma3/omega/omega;
  cof = 2.0*hbarc/omega/lambda ;
  return cof ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates complex formation zone for external medium. Omega is energy !!!

G4complex G4StrawTubeXTRadiator::GetMediumComplexFZ( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle    ) 
{
  G4double cof, length,delta, real, image;

  length = 0.5*GetMediumFormationZone(omega,gamma,varAngle);
  delta  = length*GetMediumLinearPhotoAbs(omega);
  cof    = 1.0/(1.0 + delta*delta);

  real   = length*cof;
  image  = real*delta;

  G4complex zone(real,image); 
  return zone;
}

////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// medium material

void G4StrawTubeXTRadiator::ComputeMediumPhotoAbsCof() 
{
   G4int i, j, numberOfElements;
   static const G4MaterialTable* 
   theMaterialTable = G4Material::GetMaterialTable();

   G4SandiaTable thisMaterialSandiaTable(fMatIndex3);
   numberOfElements = (*theMaterialTable)[fMatIndex3]->GetNumberOfElements();
   G4int* thisMaterialZ = new G4int[numberOfElements];

   for(i=0;i<numberOfElements;i++)
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[fMatIndex3]->
                                      GetElement(i)->GetZ() ;
   }
   fMediumIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
   fMediumIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                           (*theMaterialTable)[fMatIndex3]->GetFractionVector() ,
                             numberOfElements,fMediumIntervalNumber);
   
   fMediumPhotoAbsCof = new G4double*[fMediumIntervalNumber];

   for(i=0;i<fMediumIntervalNumber;i++)
   {
     fMediumPhotoAbsCof[i] = new G4double[5];
   }
   for(i=0;i<fMediumIntervalNumber;i++)
   {
      fMediumPhotoAbsCof[i][0] = thisMaterialSandiaTable.
                                GetPhotoAbsorpCof(i+1,0); 
                              
      for(j=1;j<5;j++)
      {
           fMediumPhotoAbsCof[i][j] = thisMaterialSandiaTable.
                                     GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex3]->GetDensity();
      }
   }
   delete[] thisMaterialZ;
   return;
}

//////////////////////////////////////////////////////////////////////
//
// Returns the value of linear photo absorption coefficient (in reciprocal 
// length) for medium for given energy of X-ray photon omega

G4double G4StrawTubeXTRadiator::GetMediumLinearPhotoAbs(G4double omega) 
{
  G4int i ;
  G4double omega2, omega3, omega4; 

  omega2 = omega*omega;
  omega3 = omega2*omega;
  omega4 = omega2*omega2;

  for(i=0;i<fMediumIntervalNumber;i++)
  {
    if( omega < fMediumPhotoAbsCof[i][0] ) break;
  }
  if( i == 0 )
  { 
    G4Exception("Invalid (<I1) energy in G4VXTRenergyLoss::GetMediumLinearPhotoAbs");
  }
  else i-- ;
  
  return fMediumPhotoAbsCof[i][1]/omega  + fMediumPhotoAbsCof[i][2]/omega2 + 
         fMediumPhotoAbsCof[i][3]/omega3 + fMediumPhotoAbsCof[i][4]/omega4  ;
}




//
//
////////////////////////////////////////////////////////////////////////////








