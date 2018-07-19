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
// Please cite the following papers if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007
// IEEE TNS 51, 4:1395-1401, 2004
//
// Based on purging magnet advanced example
//

#include "TabulatedField3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex myTabulatedField3DLock = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TabulatedField3D::TabulatedField3D(G4float gr1, G4float gr2, G4float gr3, G4float gr4, G4int choiceModel) 
{    

  G4cout << " ********************** " << G4endl;
  G4cout << " **** CONFIGURATION *** " << G4endl;
  G4cout << " ********************** " << G4endl;

  G4cout<< G4endl; 
  G4cout << "=====> You have selected :" << G4endl;
  if (choiceModel==1) G4cout<< "-> Square quadrupole field"<<G4endl;
  if (choiceModel==2) G4cout<< "-> 3D quadrupole field"<<G4endl;
  if (choiceModel==3) G4cout<< "-> Enge quadrupole field"<<G4endl;
  G4cout << "   G1 (T/m) = "<< gr1 << G4endl;
  G4cout << "   G2 (T/m) = "<< gr2 << G4endl;
  G4cout << "   G3 (T/m) = "<< gr3 << G4endl;
  G4cout << "   G4 (T/m) = "<< gr4 << G4endl;
  
  fGradient1 = gr1;
  fGradient2 = gr2;
  fGradient3 = gr3;
  fGradient4 = gr4;
  fModel = choiceModel;
  
  if (fModel==2)
  {
  //
  //This is a thread-local class and we have to avoid that all workers open the 
  //file at the same time
  G4AutoLock lock(&myTabulatedField3DLock);
  //

  const char * filename ="OM50.grid";
  
  double lenUnit= mm;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      3D Magnetic field from OPERA software "
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 
  ifstream file( filename ); // Open the file for reading.
  
  // Read table dimensions 
  file >> fNx >> fNy >> fNz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: " 
	 << fNx << " " << fNy << " " << fNz << " ] "
	 << endl;

  // Set up storage space for table
  fXField.resize( fNx );
  fYField.resize( fNx );
  fZField.resize( fNx );
  int ix, iy, iz;
  for (ix=0; ix<fNx; ix++) 
  {
    fXField[ix].resize(fNy);
    fYField[ix].resize(fNy);
    fZField[ix].resize(fNy);
    for (iy=0; iy<fNy; iy++) 
    {
      fXField[ix][iy].resize(fNz);
      fYField[ix][iy].resize(fNz);
      fZField[ix][iy].resize(fNz);
    }
  }
  
  // Read in the data
  double xval,yval,zval,bx,by,bz;
  double permeability; // Not used in this example.
  for (ix=0; ix<fNx; ix++) 
  {
    for (iy=0; iy<fNy; iy++) 
    {
      for (iz=0; iz<fNz; iz++) 
      {
        file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
        if ( ix==0 && iy==0 && iz==0 ) 
	{
          fMinix = xval * lenUnit;
          fMiniy = yval * lenUnit;
          fMiniz = zval * lenUnit;
        }
        fXField[ix][iy][iz] = bx ;
        fYField[ix][iy][iz] = by ;
        fZField[ix][iy][iz] = bz ;
      }
    }
  }
  file.close();

  //
  lock.unlock();
  //

  fMaxix = xval * lenUnit;
  fMaxiy = yval * lenUnit;
  fMaxiz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl; 

  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << fMinix/cm << " " << fMiniy/cm << " " << fMiniz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << fMaxix/cm << " " << fMaxiy/cm << " " << fMaxiz/cm << " cm " << endl;

  fDx = fMaxix - fMinix;
  fDy = fMaxiy - fMiniy;
  fDz = fMaxiz - fMiniz;

  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << fDx/cm << " " << fDy/cm << " " << fDz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;

  
  // Table normalization

  for (ix=0; ix<fNx; ix++) 
  {
    for (iy=0; iy<fNy; iy++) 
    {
      for (iz=0; iz<fNz; iz++) 
      {

	fXField[ix][iy][iz] = (fXField[ix][iy][iz]/197.736);
        fYField[ix][iy][iz] = (fYField[ix][iy][iz]/197.736);
        fZField[ix][iy][iz] = (fZField[ix][iy][iz]/197.736);

      }
    }
  }

  } // fModel==2

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void TabulatedField3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{ 
  //G4cout << fGradient1 << G4endl;
  //G4cout << fGradient2 << G4endl;
  //G4cout << fGradient3 << G4endl;
  //G4cout << fGradient4 << G4endl;
  //G4cout << "---------" << G4endl;

  G4double coef, G0;
  G0 = 0;
  
  coef=1; // for protons
  //coef=2; // for alphas

//******************************************************************

// MAP

if (fModel==2)
{
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;

  double x = point[0];
  double y = point[1];
  double z = point[2]; 

  G4int quad;
  G4double gradient[5];

  gradient[0]=fGradient1*(tesla/m)/coef;
  gradient[1]=fGradient2*(tesla/m)/coef;
  gradient[2]=fGradient3*(tesla/m)/coef; 
  gradient[3]=fGradient4*(tesla/m)/coef;
  gradient[4]=-fGradient3*(tesla/m)/coef;

  for (quad=0; quad<=4; quad++)
  {
  if ((quad+1)==1) {z = point[2] + 3720 * mm;}
  if ((quad+1)==2) {z = point[2] + 3580 * mm;}
  if ((quad+1)==3) {z = point[2] + 330  * mm;}
  if ((quad+1)==4) {z = point[2] + 190  * mm;}
  if ((quad+1)==5) {z = point[2] + 50   * mm;}

  // Check that the point is within the defined region 
       
  if 
  (
    x>=fMinix && x<=fMaxix &&
    y>=fMiniy && y<=fMaxiy &&
    z>=fMiniz && z<=fMaxiz 
  ) 
  {
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - fMinix) / fDx;
    double yfraction = (y - fMiniy) / fDy; 
    double zfraction = (z - fMiniz) / fDz;

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(fNx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(fNy-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(fNz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    
    //int xindex = static_cast<int>(xdindex);
    //int yindex = static_cast<int>(ydindex);
    //int zindex = static_cast<int>(zdindex);

    // SI 15/12/2016: modified according to bugzilla 1879
    int xindex = static_cast<int>(std::floor(xdindex));
    int yindex = static_cast<int>(std::floor(ydindex));
    int zindex = static_cast<int>(std::floor(zdindex));
    
    // Interpolated field
    Bfield[0] =
     (fXField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      fXField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      fXField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      fXField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      fXField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      fXField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      fXField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      fXField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad]
      + Bfield[0];
      
    Bfield[1] =
     (fYField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      fYField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      fYField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      fYField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      fYField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      fYField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      fYField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      fYField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad] 
      + Bfield[1];

    Bfield[2] =
     (fZField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      fZField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      fZField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      fZField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      fZField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      fZField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      fZField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      fZField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad]
      + Bfield[2];

     } 

} // loop on quads

} //end  MAP


//******************************************************************
// SQUARE

if (fModel==1)
{
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;

  // Field components 
  G4double Bx = 0;
  G4double By = 0;
  G4double Bz = 0;
   
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2];

  if (z>=-3770*mm && z<=-3670*mm)  G0 = (fGradient1/coef)* tesla/m;
  if (z>=-3630*mm && z<=-3530*mm)  G0 = (fGradient2/coef)* tesla/m;
  
  if (z>=-380*mm  && z<=-280*mm)   G0 = (fGradient3/coef)* tesla/m;
  if (z>=-240*mm  && z<=-140*mm)   G0 = (fGradient4/coef)* tesla/m;
  if (z>=-100*mm  && z<=0*mm)      G0 = (-fGradient3/coef)* tesla/m;

  Bx = y*G0;
  By = x*G0;
  Bz = 0;

  Bfield[0] = Bx;
  Bfield[1] = By;
  Bfield[2] = Bz;

}

// end SQUARE

//******************************************************************
// ENGE

if (fModel==3)
{

  // X POSITION OF FIRST QUADRUPOLE
  // G4double lineX = 0*mm;

  // Z POSITION OF FIRST QUADRUPOLE
  G4double lineZ = -3720*mm;

  // QUADRUPOLE HALF LENGTH
  // G4double quadHalfLength = 50*mm;
  
  // QUADRUPOLE CENTER COORDINATES
  G4double zoprime;
  
  G4double Grad1, Grad2, Grad3, Grad4, Grad5;
  Grad1=fGradient1;
  Grad2=fGradient2;
  Grad3=fGradient3;
  Grad4=fGradient4;
  Grad5=-Grad3;  

  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;

  double x = point[0];
  double y = point[1];
  double z = point[2]; 

  if ( (z>=-3900*mm && z<-3470*mm)  || (z>=-490*mm && z<100*mm) )
  {
  G4double Bx=0;
  G4double By=0; 
  G4double Bz=0;
  
  // FRINGING FILED CONSTANTS
  G4double c0[5], c1[5], c2[5], z1[5], z2[5], a0[5], gradient[5];
  
  // DOUBLET***************
  
  // QUADRUPOLE 1
  c0[0] = -10.;			// Ci are constants in Pn(z)=C0+C1*s+C2*s^2
  c1[0] = 3.08874;
  c2[0] = -0.00618654;
  z1[0] = 28.6834*mm;		// Fringing field lower limit
  z2[0] = z1[0]+50*mm;		// Fringing field upper limit	
  a0[0] = 7.5*mm;             	// Bore Radius
  gradient[0] =Grad1*(tesla/m)/coef;           

  // QUADRUPOLE 2
  c0[1] = -10.;			// Ci are constants in Pn(z)=C0+C1*s+C2*s^2
  c1[1] = 3.08874;
  c2[1] = -0.00618654;
  z1[1] = 28.6834*mm; 		// Fringing field lower limit
  z2[1] = z1[1]+50*mm;		// Fringing field upper limit
  a0[1] = 7.5*mm;             	// Bore Radius
  gradient[1] =Grad2*(tesla/m)/coef;
    
  // TRIPLET**********

  // QUADRUPOLE 3
  c0[2] = -10.;			// Ci are constants in Pn(z)=C0+C1*s+C2*s^2
  c1[2] = 3.08874;
  c2[2] = -0.00618654;
  z1[2] = 28.6834*mm;		// Fringing field lower limit
  z2[2] = z1[2]+50*mm;		// Fringing field upper limit
  a0[2] = 7.5*mm;             	// Bore Radius
  gradient[2] = Grad3*(tesla/m)/coef;

  // QUADRUPOLE 4
  c0[3] = -10.;			// Ci are constants in Pn(z)=C0+C1*s+C2*s^2
  c1[3] = 3.08874;
  c2[3] = -0.00618654;
  z1[3] = 28.6834*mm;		// Fringing field lower limit
  z2[3] = z1[3]+50*mm;		// Fringing field upper limit
  a0[3] = 7.5*mm;             	// Bore Radius
  gradient[3] = Grad4*(tesla/m)/coef;
  
   // QUADRUPOLE 5
  c0[4] = -10.;			// Ci are constants in Pn(z)=C0+C1*s+C2*s^2
  c1[4] = 3.08874;
  c2[4] = -0.00618654;
  z1[4] = 28.6834*mm;		// Fringing field lower limit
  z2[4] = z1[4]+50*mm;		// Fringing field upper limit
  a0[4] = 7.5*mm;             	// Bore Radius
  gradient[4] = Grad5*(tesla/m)/coef;  

  // FIELD CREATED BY A QUADRUPOLE IN ITS LOCAL FRAME
  G4double Bx_local,By_local,Bz_local;
  Bx_local = 0; By_local = 0; Bz_local = 0;
  
  // QUADRUPOLE FRAME
  G4double x_local,y_local,z_local;
  x_local= 0; y_local=0; z_local=0;

  G4double myVars = 0;          // For Enge formula
  G4double G1, G2, G3;	        // For Enge formula
  G4double K1, K2, K3; 		// For Enge formula
  G4double P0, P1, P2, cte;	// For Enge formula

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
  
  for (G4int i=0;i<5; i++) // LOOP ON MAGNETS
  {
 
	 if (i<2) // (if Doublet)
	 {	
	   zoprime = lineZ + i*140*mm; // centre of magnet nbr i 
	   x_local = x; 
	   y_local = y; 
	   z_local = (z - zoprime);
	 }
	 else    // else the current magnet is in the triplet
	 {
	   zoprime = lineZ + i*140*mm +(3150-40)*mm;

	   x_local = x; 
	   y_local = y; 
	   z_local = (z - zoprime);
	
	 }					
	 
	 if ( z_local < -z2[i] || z_local > z2[i])  // Outside the fringing field
	 {
	  G0=0;
	  G1=0;
	  G2=0;
	  G3=0;
	 }
	 
	 if ( (z_local>=-z1[i]) && (z_local<=z1[i]) ) // inside the quadrupole but outside the fringefield
	 {
	  G0=gradient[i];
	  G1=0;
	  G2=0;
	  G3=0;
	 }
	 
	 if ( ((z_local>=-z2[i]) && (z_local<-z1[i])) ||  ((z_local>z1[i]) && (z_local<=z2[i])) ) // inside the fringefield
	 {

          myVars = ( z_local - z1[i]) / a0[i];     // se (8) p1397 TNS 51
          if (z_local<-z1[i])  myVars = ( - z_local - z1[i]) / a0[i];  // see (9) p1397 TNS 51


	  P0 = c0[i]+c1[i]*myVars+c2[i]*myVars*myVars;

	  P1 = c1[i]/a0[i]+2*c2[i]*(z_local-z1[i])/a0[i]/a0[i]; // dP/fDz
	  if (z_local<-z1[i])  P1 = -c1[i]/a0[i]+2*c2[i]*(z_local+z1[i])/a0[i]/a0[i];

	  P2 = 2*c2[i]/a0[i]/a0[i]; 	// d2P/fDz2

	  cte = 1 + G4Exp(c0[i]);   // (1+e^c0)

	  K1 = -cte*P1*G4Exp(P0)/( (1+G4Exp(P0))*(1+G4Exp(P0)) );  // see (11) p1397 TNS 51

	  K2 = -cte*G4Exp(P0)*(					// see (12) p1397 TNS 51
	   P2/( (1+G4Exp(P0))*(1+G4Exp(P0)) )
	   +2*P1*K1/(1+G4Exp(P0))/cte
	   +P1*P1/(1+G4Exp(P0))/(1+G4Exp(P0))
	   );                                                            
 
	  K3 = -cte*G4Exp(P0)*(				// see (13) p1397 TNS 51	
	   (3*P2*P1+P1*P1*P1)/(1+G4Exp(P0))/(1+G4Exp(P0))
	   +4*K1*(P1*P1+P2)/(1+G4Exp(P0))/cte
	   +2*P1*(K1*K1/cte/cte+K2/(1+G4Exp(P0))/cte)
	   );
	  
	  G0 = gradient[i]*cte/(1+G4Exp(P0));    	// G = G0*K(z) see (7) p1397 TNS 51
	  G1 = gradient[i]*K1;				// dG/fDz
	  G2 = gradient[i]*K2;				// d2G/fDz2
	  G3 = gradient[i]*K3;				// d3G/fDz3

	 }
	  
	 Bx_local = y_local*(G0-(1./12)*(3*x_local*x_local+y_local*y_local)*G2); 	// see (4) p1396 TNS 51
	 By_local = x_local*(G0-(1./12)*(3*y_local*y_local+x_local*x_local)*G2);	// see (5) p1396 TNS 51
	 Bz_local = x_local*y_local*(G1-(1./12)*(x_local*x_local+y_local*y_local)*G3);	// see (6) p1396 TNS 51

	 // TOTAL MAGNETIC FIELD
	 
	 Bx = Bx + Bx_local ;
	 By = By + By_local ;
	 Bz = Bz + Bz_local ;


  } // LOOP ON QUADRUPOLES 
  
  Bfield[0] = Bx;
  Bfield[1] = By;
  Bfield[2] = Bz;
  }
  
		        
} // end ENGE

}
