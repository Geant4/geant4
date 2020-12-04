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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11
//
// Based on purging magnet advanced example.
//

#include "EMField.hh"
#include "G4Exp.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


EMField::EMField() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EMField::GetFieldValue(const double point[4], double *Bfield ) const
{ 
  // Magnetic field
  Bfield[0] = 0;
  Bfield[1] = 0;
  Bfield[2] = 0;
  
  // Electric field
  Bfield[3] = 0;
  Bfield[4] = 0;
  Bfield[5] = 0;

  G4double Bx = 0;
  G4double By = 0;
  G4double Bz = 0;
   
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2];

// ***********************
// AIFIRA SWITCHING MAGNET
// ***********************
  
  // MAGNETIC FIELD VALUE FOR 3 MeV ALPHAS
  //  G4double switchingField = 0.0589768635 * tesla ;
  G4double switchingField =   0.0590201 * tesla ;
  
  // BEAM START
  G4double beamStart = -10*m;

  // RADIUS
  G4double Rp = 0.698*m;

  // ENTRANCE POSITION AFTER ANALYSIS MAGNET
  G4double zS = 975*mm;
  
  // POLE GAP
  G4double D = 31.8*mm;
  
  // FRINGING FIELD

  G4double fieldBoundary, wc0, wc1, wc2, wc3, limitMinEntrance, limitMaxEntrance, limitMinExit, limitMaxExit;

  limitMinEntrance = beamStart+zS-4*D;
  limitMaxEntrance = beamStart+zS+4*D;
  limitMinExit =Rp-4*D;
  limitMaxExit =Rp+4*D;  
    
  wc0 = 0.3835;
  wc1 = 2.388;
  wc2 = -0.8171;
  wc3 = 0.200;

  fieldBoundary=0.62;

  G4double ws, largeS, h, dhdlargeS, dhds, dlargeSds, dsdz, dsdx, zs0, Rs0, xcenter, zcenter;
  
// - ENTRANCE OF SWITCHING MAGNET

if ( (z >= limitMinEntrance) && (z < limitMaxEntrance) ) 
{
  zs0 = fieldBoundary*D;
  ws = (-z+beamStart+zS-zs0)/D;
  dsdz = -1/D;
  dsdx = 0;

  largeS = wc0 + wc1*ws + wc2*ws*ws + wc3*ws*ws*ws;
  h = 1./(1.+G4Exp(largeS));
  dhdlargeS = -G4Exp(largeS)*h*h;  
  dlargeSds = wc1+ 2*wc2*ws + 3*wc3*ws*ws;
  dhds = dhdlargeS * dlargeSds;
      
  By = switchingField * h ;
  Bx = y*switchingField*dhds*dsdx;
  Bz = y*switchingField*dhds*dsdz;

}

// - HEART OF SWITCHING MAGNET 	  
		
 if ( 
          (z >= limitMaxEntrance)  
     &&   (( x*x + (z -(beamStart+zS))*(z -(beamStart+zS)) < limitMinExit*limitMinExit)) 
    )  	
{
   Bx=0; 
   By = switchingField; 
   Bz=0;
}			                    
	
// - EXIT OF SWITCHING MAGNET

if ( 
        (z >= limitMaxEntrance)  
     && (( x*x + (z -(beamStart+zS))*(z -(beamStart+zS))) >= limitMinExit*limitMinExit) 
     && (( x*x + (z -(beamStart+zS))*(z -(beamStart+zS))) < limitMaxExit*limitMaxExit)

   )  	
{

  xcenter = 0;
  zcenter =  beamStart+zS;
  
  Rs0 = Rp + D*fieldBoundary;
  ws = (std::sqrt((z-zcenter)*(z-zcenter)+(x-xcenter)*(x-xcenter)) - Rs0)/D;
  	
  dsdz = (1/D)*(z-zcenter)/std::sqrt((z-zcenter)*(z-zcenter)+(x-xcenter)*(x-xcenter));
  dsdx = (1/D)*(x-xcenter)/std::sqrt((z-zcenter)*(z-zcenter)+(x-xcenter)*(x-xcenter));

  largeS = wc0 + wc1*ws + wc2*ws*ws + wc3*ws*ws*ws;
  h = 1./(1.+G4Exp(largeS));
  dhdlargeS = -G4Exp(largeS)*h*h;  
  dlargeSds = wc1+ 2*wc2*ws + 3*wc3*ws*ws;
  dhds = dhdlargeS * dlargeSds;
      
  By = switchingField * h ;
  Bx = y*switchingField*dhds*dsdx;
  Bz = y*switchingField*dhds*dsdz;

}

// **************************
// MICROBEAM LINE QUADRUPOLES
// **************************
 
  // MICROBEAM LINE ANGLE
  G4double lineAngle = -10*deg;
  
  // X POSITION OF FIRST QUADRUPOLE
  G4double lineX = -1295.59*mm;

  // Z POSITION OF FIRST QUADRUPOLE
  G4double lineZ = -1327*mm;

  // Adjust magnetic zone absolute position
  lineX = lineX + 5.24*micrometer*std::cos(-lineAngle); // 5.24 = 1.3 + 3.94 micrometer (cf. DetectorConstruction)
  lineZ = lineZ + 5.24*micrometer*std::sin(-lineAngle);
       
  // QUADRUPOLE HALF LENGTH
  G4double quadHalfLength = 75*mm;
  
  // QUADRUPOLE SPACING
  G4double quadSpacing = 40*mm;
  
  // QUADRUPOLE CENTER COORDINATES
  G4double xoprime, zoprime;
  
if (z>=-1400*mm && z <-200*mm)
{
  Bx=0; By=0; Bz=0;
  
  // FRINGING FILED CONSTANTS
  G4double c0[4], c1[4], c2[4], z1[4], z2[4], a0[4], gradient[4];
  
  // QUADRUPOLE 1
  c0[0] = -5.;
  c1[0] = 2.5;
  c2[0] = -0.1;
  z1[0] = 60*mm;
  z2[0] = 130*mm;
  a0[0] = 10*mm;
  gradient[0] = 3.406526 *tesla/m;

  // QUADRUPOLE 2
  c0[1] = -5.;
  c1[1] = 2.5;
  c2[1] = -0.1;
  z1[1] = 60*mm;
  z2[1] = 130*mm;
  a0[1] = 10*mm;
  gradient[1] = -8.505263 *tesla/m;

  // QUADRUPOLE 3
  c0[2] = -5.;
  c1[2] = 2.5;
  c2[2] = -0.1;
  z1[2] = 60*mm;
  z2[2] = 130*mm;
  a0[2] = 10*mm;
  gradient[2] = 8.505263 *tesla/m;

  // QUADRUPOLE 4
  c0[3] = -5.;
  c1[3] = 2.5;
  c2[3] = -0.1;
  z1[3] = 60*mm;
  z2[3] = 130*mm;
  a0[3] = 10*mm;
  gradient[3] = -3.406526*tesla/m;

  // FIELD CREATED BY A QUADRUPOLE IN ITS LOCAL FRAME
  G4double Bx_local,By_local,Bz_local;
  Bx_local = 0; By_local = 0; Bz_local = 0;
  
  // FIELD CREATED BY A QUADRUPOOLE IN WORLD FRAME
  G4double Bx_quad,By_quad,Bz_quad;
  Bx_quad = 0; By_quad=0; Bz_quad=0;
  
  // QUADRUPOLE FRAME
  G4double x_local,y_local,z_local;
  x_local= 0; y_local=0; z_local=0;

  G4double vars = 0;
  G4double G0, G1, G2, G3;
  G4double K1, K2, K3;
  G4double P0, P1, P2,     cte;

  K1=0;
  K2=0;
  K3=0;
  P0=0;
  P1=0;
  P2=0;
  G0=0;
  G1=0;
  G2=0;
  G3=0;
  cte=0;

  G4bool largeScattering=false;
  
  for (G4int i=0;i<4; i++) 
  {
 
	 if (i==0) 
	 	{	xoprime = lineX + quadHalfLength*std::sin(lineAngle);
			zoprime = lineZ + quadHalfLength*std::cos(lineAngle);

			x_local = (x - xoprime) * std::cos (lineAngle) - (z - zoprime) * std::sin (lineAngle); 
			y_local = y; 
			z_local = (z - zoprime) * std::cos (lineAngle) + (x - xoprime) * std::sin (lineAngle); 
			if (std::sqrt(x_local*x_local+y_local*y_local)>a0[i]) largeScattering=true;

		}
		 
	 if (i==1) 
	 	{	xoprime = lineX + (3*quadHalfLength+quadSpacing)*std::sin(lineAngle);
			zoprime = lineZ + (3*quadHalfLength+quadSpacing)*std::cos(lineAngle);

			x_local = (x - xoprime) * std::cos (lineAngle) - (z - zoprime) * std::sin (lineAngle); 
			y_local = y; 
			z_local = (z - zoprime) * std::cos (lineAngle) + (x - xoprime) * std::sin (lineAngle); 
			if (std::sqrt(x_local*x_local+y_local*y_local)>a0[i]) largeScattering=true;
		}

	 if (i==2) 
	 	{	xoprime = lineX + (5*quadHalfLength+2*quadSpacing)*std::sin(lineAngle);
			zoprime = lineZ + (5*quadHalfLength+2*quadSpacing)*std::cos(lineAngle);

			x_local = (x - xoprime) * std::cos (lineAngle) - (z - zoprime) * std::sin (lineAngle); 
			y_local = y; 
			z_local = (z - zoprime) * std::cos (lineAngle) + (x - xoprime) * std::sin (lineAngle); 
			if (std::sqrt(x_local*x_local+y_local*y_local)>a0[i]) largeScattering=true;
		}
	 
	 if (i==3) 
	 	{	xoprime = lineX + (7*quadHalfLength+3*quadSpacing)*std::sin(lineAngle);
			zoprime = lineZ + (7*quadHalfLength+3*quadSpacing)*std::cos(lineAngle);

			x_local = (x - xoprime) * std::cos (lineAngle) - (z - zoprime) * std::sin (lineAngle); 
			y_local = y; 
			z_local = (z - zoprime) * std::cos (lineAngle) + (x - xoprime) * std::sin (lineAngle); 
			if (std::sqrt(x_local*x_local+y_local*y_local)>a0[i]) largeScattering=true;
		}

	 
	 if ( z_local < -z2[i] )
	 {
	  G0=0;
	  G1=0;
	  G2=0;
	  G3=0;
	 }
	 
	 if ( z_local > z2[i] )
	 {
	  G0=0;
	  G1=0;
	  G2=0;
	  G3=0;
	 }

	 if ( (z_local>=-z1[i]) & (z_local<=z1[i]) ) 
	 {
	  G0=gradient[i];
	  G1=0;
	  G2=0;
	  G3=0;
	 }
	 
	 if ( ((z_local>=-z2[i]) & (z_local<-z1[i])) ||  ((z_local>z1[i]) & (z_local<=z2[i])) ) 
	 {

	  vars = ( z_local - z1[i]) / a0[i] ;
  	  if (z_local<-z1[i]) vars = ( - z_local - z1[i]) / a0[i] ;


	  P0 = c0[i]+c1[i]*vars+c2[i]*vars*vars;

	  P1 = c1[i]/a0[i]+2*c2[i]*(z_local-z1[i])/a0[i]/a0[i];
	  if (z_local<-z1[i])  P1 = -c1[i]/a0[i]+2*c2[i]*(z_local+z1[i])/a0[i]/a0[i];

	  P2 = 2*c2[i]/a0[i]/a0[i];

	  cte = 1 + G4Exp(c0[i]);

	  K1 = -cte*P1*G4Exp(P0)/( (1+G4Exp(P0))*(1+G4Exp(P0)) );

	  K2 = -cte*G4Exp(P0)*(
	   P2/( (1+G4Exp(P0))*(1+G4Exp(P0)) )
	  +2*P1*K1/(1+G4Exp(P0))/cte
	  +P1*P1/(1+G4Exp(P0))/(1+G4Exp(P0))
	  );
 
	  K3 = -cte*G4Exp(P0)*(
	  (3*P2*P1+P1*P1*P1)/(1+G4Exp(P0))/(1+G4Exp(P0))
	  +4*K1*(P1*P1+P2)/(1+G4Exp(P0))/cte
	  +2*P1*(K1*K1/cte/cte+K2/(1+G4Exp(P0))/cte)
	   );
	  
	  G0 = gradient[i]*cte/(1+G4Exp(P0));
	  G1 = gradient[i]*K1;
	  G2 = gradient[i]*K2;
	  G3 = gradient[i]*K3;

	 }
	  
	 // PROTECTION AGAINST LARGE SCATTERING

	 if ( largeScattering ) 
	 {
	  G0=0;
	  G1=0;
	  G2=0;
	  G3=0;
	 }

	 // MAGNETIC FIELD COMPUTATION FOR EACH QUADRUPOLE
	 
	 Bx_local = y_local*(G0-(1./12)*(3*x_local*x_local+y_local*y_local)*G2);
	 By_local = x_local*(G0-(1./12)*(3*y_local*y_local+x_local*x_local)*G2);
	 Bz_local = x_local*y_local*(G1-(1./12)*(x_local*x_local+y_local*y_local)*G3);

	 Bx_quad = Bz_local*std::sin(lineAngle)+Bx_local*std::cos(lineAngle);
	 By_quad = By_local;
	 Bz_quad = Bz_local*std::cos(lineAngle)-Bx_local*std::sin(lineAngle);

	 // TOTAL MAGNETIC FIELD
	 
	 Bx = Bx + Bx_quad ;
	 By = By + By_quad ;
	 Bz = Bz + Bz_quad ;

  } // LOOP ON QUADRUPOLES

      
} // END OF QUADRUPLET

  Bfield[0] = Bx;
  Bfield[1] = By;
  Bfield[2] = Bz;

// *****************************************
// ELECTRIC FIELD CREATED BY SCANNING PLATES
// *****************************************

  Bfield[3] = 0;
  Bfield[4] = 0;
  Bfield[5] = 0;

  // POSITION OF EXIT OF LAST QUAD WHERE THE SCANNING PLATES START

  G4double electricPlateWidth1 = 5 * mm;
  G4double electricPlateWidth2 = 5 * mm;
  G4double electricPlateLength1 = 36 * mm;
  G4double electricPlateLength2 = 34 * mm;
  G4double electricPlateGap = 5 * mm;
  G4double electricPlateSpacing1 = 3 * mm;
  G4double electricPlateSpacing2 = 4 * mm;

  // APPLY VOLTAGE HERE IN VOLTS (no electrostatic deflection here)
  G4double electricPlateVoltage1 = 0 * volt;
  G4double electricPlateVoltage2 = 0 * volt;

  G4double electricFieldPlate1 = electricPlateVoltage1 / electricPlateSpacing1 ;
  G4double electricFieldPlate2 = electricPlateVoltage2 / electricPlateSpacing2 ;

  G4double  beginFirstZoneX = lineX + (8*quadHalfLength+3*quadSpacing)*std::sin(lineAngle);
  G4double  beginFirstZoneZ = lineZ + (8*quadHalfLength+3*quadSpacing)*std::cos(lineAngle);

  G4double  beginSecondZoneX = lineX + (8*quadHalfLength+3*quadSpacing+electricPlateLength1+electricPlateGap)*std::sin(lineAngle);
  G4double  beginSecondZoneZ = lineZ + (8*quadHalfLength+3*quadSpacing+electricPlateLength1+electricPlateGap)*std::cos(lineAngle);

  G4double xA, zA, xB, zB, xC, zC, xD, zD;
  G4double slope1, cte1, slope2, cte2, slope3, cte3, slope4, cte4;
 
  // WARNING : lineAngle < 0

  // FIRST PLATES
  
  xA = beginFirstZoneX + std::cos(lineAngle)*electricPlateSpacing1/2;
  zA = beginFirstZoneZ - std::sin(lineAngle)*electricPlateSpacing1/2;

  xB = xA + std::sin(lineAngle)*electricPlateLength1; 
  zB = zA + std::cos(lineAngle)*electricPlateLength1;
  
  xC = xB - std::cos(lineAngle)*electricPlateSpacing1;
  zC = zB + std::sin(lineAngle)*electricPlateSpacing1;

  xD = xC - std::sin(lineAngle)*electricPlateLength1; 
  zD = zC - std::cos(lineAngle)*electricPlateLength1;
  
  slope1 = (xB-xA)/(zB-zA);
  cte1 = xA - slope1 * zA;
  
  slope2 = (xC-xB)/(zC-zB);
  cte2 = xB - slope2 * zB;
  
  slope3 = (xD-xC)/(zD-zC);
  cte3 = xC - slope3 * zC;
  
  slope4 = (xA-xD)/(zA-zD);
  cte4 = xD - slope4 * zD;
  
   
  if 
  (
       x <= slope1 * z + cte1
    && x >= slope3 * z + cte3
    && x <= slope4 * z + cte4
    && x >= slope2 * z + cte2    
    && std::abs(y)<=electricPlateWidth1/2
  )  

  {
      Bfield[3] = electricFieldPlate1*std::cos(lineAngle);
      Bfield[4] = 0;
      Bfield[5] = -electricFieldPlate1*std::sin(lineAngle);
 
  }
      
  // SECOND PLATES
      
  xA = beginSecondZoneX + std::cos(lineAngle)*electricPlateWidth2/2;
  zA = beginSecondZoneZ - std::sin(lineAngle)*electricPlateWidth2/2;

  xB = xA + std::sin(lineAngle)*electricPlateLength2; 
  zB = zA + std::cos(lineAngle)*electricPlateLength2;
  
  xC = xB - std::cos(lineAngle)*electricPlateWidth2;
  zC = zB + std::sin(lineAngle)*electricPlateWidth2;

  xD = xC - std::sin(lineAngle)*electricPlateLength2; 
  zD = zC - std::cos(lineAngle)*electricPlateLength2;
  
  slope1 = (xB-xA)/(zB-zA);
  cte1 = xA - slope1 * zA;
  
  slope2 = (xC-xB)/(zC-zB);
  cte2 = xB - slope2 * zB;
  
  slope3 = (xD-xC)/(zD-zC);
  cte3 = xC - slope3 * zC;
  
  slope4 = (xA-xD)/(zA-zD);
  cte4 = xD - slope4 * zD;

  if 
  (     
       x <= slope1 * z + cte1
    && x >= slope3 * z + cte3
    && x <= slope4 * z + cte4
    && x >= slope2 * z + cte2    
    && std::abs(y)<=electricPlateSpacing2/2
  )

  {  
      Bfield[3] = 0;
      Bfield[4] = electricFieldPlate2;
      Bfield[5] = 0;
  }

//

}
