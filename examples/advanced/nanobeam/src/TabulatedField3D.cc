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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "TabulatedField3D.hh"
#include "G4SystemOfUnits.hh"

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
  
  gradient1 = gr1;
  gradient2 = gr2;
  gradient3 = gr3;
  gradient4 = gr4;
  model = choiceModel;
  
  if (model==2)
  {
  const char * filename ="OM50.grid";
  
  double lenUnit= mm;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      3D Magnetic field from OPERA software "
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 
  ifstream file( filename ); // Open the file for reading.
  
  // Read table dimensions 
  file >> nx >> ny >> nz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: " 
	 << nx << " " << ny << " " << nz << " ] "
	 << endl;

  // Set up storage space for table
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }
  
  // Read in the data
  double xval,yval,zval,bx,by,bz;
  double permeability; // Not used in this example.
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
        file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx ;
        yField[ix][iy][iz] = by ;
        zField[ix][iy][iz] = bz ;
      }
    }
  }
  file.close();

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl; 
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm " << endl;

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;

  
  // Table normalization
  for (ix=0; ix<nx; ix++) 
  {
    for (iy=0; iy<ny; iy++) 
    {
      for (iz=0; iz<nz; iz++) 
      {

	xField[ix][iy][iz] = (xField[ix][iy][iz]/197.736);
        yField[ix][iy][iz] = (yField[ix][iy][iz]/197.736);
        zField[ix][iy][iz] = (zField[ix][iy][iz]/197.736);

	}
    }
  }

  } // model==2

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void TabulatedField3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{ 

  G4double coef, G0;
  G0 = 0;
  
  coef=1; //protons
  //coef=2; // alphas

//******************************************************************

// MAP
if (model==2)
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

  gradient[0]=gradient1*(tesla/m)/coef;
  gradient[1]=gradient2*(tesla/m)/coef;
  gradient[2]=gradient3*(tesla/m)/coef; 
  gradient[3]=gradient4*(tesla/m)/coef;
  gradient[4]=-gradient3*(tesla/m)/coef;

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
    x>=minx && x<=maxx &&
    y>=miny && y<=maxy &&
    z>=minz && z<=maxz 
  ) 
  {
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    
    // Interpolated field
    Bfield[0] =
     (xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad]
      + Bfield[0];
      
    Bfield[1] =
     (yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad] 
      + Bfield[1];

    Bfield[2] =
     (zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)*gradient[quad]
      + Bfield[2];

     } 

} // loop on quads

} //end  MAP


//******************************************************************
// SQUARE

if (model==1)
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

  if (z>=-3770*mm && z<=-3670*mm)  G0 = (gradient1/coef)* tesla/m;
  if (z>=-3630*mm && z<=-3530*mm)  G0 = (gradient2/coef)* tesla/m;
  
  if (z>=-380*mm  && z<=-280*mm)   G0 = (gradient3/coef)* tesla/m;
  if (z>=-240*mm  && z<=-140*mm)   G0 = (gradient4/coef)* tesla/m;
  if (z>=-100*mm  && z<=0*mm)      G0 = (-gradient3/coef)* tesla/m;

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

if (model==3)
{

  // X POSITION OF FIRST QUADRUPOLE
  //G4double lineX = 0*mm;

  // Z POSITION OF FIRST QUADRUPOLE
  G4double lineZ = -3720*mm;

  // QUADRUPOLE HALF LENGTH
  //G4double quadHalfLength = 50*mm;
  
  // QUADRUPOLE CENTER COORDINATES
  G4double zoprime;
  
  G4double Grad1, Grad2, Grad3, Grad4, Grad5;
  Grad1=gradient1;
  Grad2=gradient2;
  Grad3=gradient3;
  Grad4=gradient4;
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
  c0[0] = -10.;			// Ci are constants in the; Pn(z)=C0+C1*s+C2*s^2
  c1[0] = 3.08874;
  c2[0] = -0.00618654;
  z1[0] = 28.6834*mm;		//Fringing field lower limit
  z2[0] = z1[0]+50*mm;		//Fringing field upper limit	
  a0[0] = 7.5*mm;             	//Bore Radius
  gradient[0] =Grad1*(tesla/m)/coef;           

  // QUADRUPOLE 2
  c0[1] = -10.;			// Ci are constants in the; Pn(z)=C0+C1*s+C2*s^2
  c1[1] = 3.08874;
  c2[1] = -0.00618654;
  z1[1] = 28.6834*mm; 		//Fringing field lower limit
  z2[1] = z1[1]+50*mm;		//Fringing field upper limit
  a0[1] = 7.5*mm;             	//Bore Radius
  gradient[1] =Grad2*(tesla/m)/coef;
    
  // TRIPLET**********
  // QUADRUPOLE 3
  c0[2] = -10.;			// Ci are constants in the; Pn(z)=C0+C1*s+C2*s^2
  c1[2] = 3.08874;
  c2[2] = -0.00618654;
  z1[2] = 28.6834*mm;		//Fringing field lower limit
  z2[2] = z1[2]+50*mm;		//Fringing field upper limit
  a0[2] = 7.5*mm;             	//Bore Radius
  gradient[2] = Grad3*(tesla/m)/coef;

  // QUADRUPOLE 4
  c0[3] = -10.;			// Ci are constants in the; Pn(z)=C0+C1*s+C2*s^2
  c1[3] = 3.08874;
  c2[3] = -0.00618654;
  z1[3] = 28.6834*mm;		//Fringing field lower limit
  z2[3] = z1[3]+50*mm;		//Fringing field upper limit
  a0[3] = 7.5*mm;             	//Bore Radius
  gradient[3] = Grad4*(tesla/m)/coef;
  
   // QUADRUPOLE 5
  c0[4] = -10.;			// Ci are constants in the; Pn(z)=C0+C1*s+C2*s^2
  c1[4] = 3.08874;
  c2[4] = -0.00618654;
  z1[4] = 28.6834*mm;		//Fringing field lower limit
  z2[4] = z1[4]+50*mm;		//Fringing field upper limit
  a0[4] = 7.5*mm;             	//Bore Radius
  gradient[4] = Grad5*(tesla/m)/coef;  

  // FIELD CREATED BY A QUADRUPOLE IN ITS LOCAL FRAME
  G4double Bx_local,By_local,Bz_local;
  Bx_local = 0; By_local = 0; Bz_local = 0;
  
  // FIELD CREATED BY A QUADRUPOOLE IN WORLD FRAME
  //unused G4double Bx_quad,By_quad,Bz_quad;
  //unsued Bx_quad = 0; By_quad=0; Bz_quad=0;
  
  // QUADRUPOLE FRAME
  G4double x_local,y_local,z_local;
  x_local= 0; y_local=0; z_local=0;

  //G4double vars = 0;               // For Enges formula
  G4double G1, G2, G3;	        // For Enges formula
  //G4double K0, K1, K2, K3;	// For Enges formula
  //G4double P0, P1, P2, P3, cte;	// For Enges formula
  G4double     K1, K2, K3;	// For Enges formula
  G4double P0, P1, P2,     cte;	// For Enges formula

  //K0=0;
  K1=0;
  K2=0;
  K3=0;
  P0=0;
  P1=0;
  P2=0;
  //P3=0;
  G0=0;
  G1=0;
  G2=0;
  G3=0;
  cte=0;
  
  for (G4int i=0;i<5; i++) // LOOP ON MAGNETS
  {
 
	if (i<2) // (if Doublet)
	 	{	
			zoprime = lineZ + i*140*mm; //centre of magnet nbr i 
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
	 

	 if ( (z_local>=-z1[i]) & (z_local<=z1[i]) ) // inside the quadrupole but outside the fringefield
	 {
	  G0=gradient[i];
	  G1=0;
	  G2=0;
	  G3=0;
	 }
	 
	 if ( ((z_local>=-z2[i]) & (z_local<-z1[i])) ||  ((z_local>z1[i]) & (z_local<=z2[i])) ) // Inside the fringefield
	 {

     //vars = ( z_local - z1[i]) / a0[i];     // se (8)
     //if (z_local<-z1[i])  vars = ( - z_local - z1[i]) / a0[i];  // se (9)  p1397 Incerti et.al.


	  P0 = c0[i]+c1[i]*s+c2[i]*s*s;

	  P1 = c1[i]/a0[i]+2*c2[i]*(z_local-z1[i])/a0[i]/a0[i];       //dP/dz
	  if (z_local<-z1[i])  P1 = -c1[i]/a0[i]+2*c2[i]*(z_local+z1[i])/a0[i]/a0[i];  // --"--

	  P2 = 2*c2[i]/a0[i]/a0[i]; 	//   d2P/dz2

	  //P3 = 0;    			//  d3P/dw3 ??

	  cte = 1 + std::exp(c0[i]);   // (1+e^c0)

	  K1 = -cte*P1*std::exp(P0)/( (1+std::exp(P0))*(1+std::exp(P0)) );  // se (11) p1397 Incerti et.al.

	 K2 = -cte*std::exp(P0)*(					// se (12) p1397 Incerti et.al.
	  P2/( (1+std::exp(P0))*(1+std::exp(P0)) )
	 +2*P1*K1/(1+std::exp(P0))/cte
	 +P1*P1/(1+std::exp(P0))/(1+std::exp(P0))
	 );                                                            
 
	 K3 = -cte*std::exp(P0)*(				// se (13) p1397 Incerti et.al	
	 (3*P2*P1+P1*P1*P1)/(1+std::exp(P0))/(1+std::exp(P0))
	 +4*K1*(P1*P1+P2)/(1+std::exp(P0))/cte
	 +2*P1*(K1*K1/cte/cte+K2/(1+std::exp(P0))/cte)
	  );
	  
	 G0 = gradient[i]*cte/(1+std::exp(P0));    	// G = G0*K(z) , se (7) p1397 Incerti et.al
	 G1 = gradient[i]*K1;				// dG/dz
	 G2 = gradient[i]*K2;				// d2G/dz2
	 G3 = gradient[i]*K3;				// d3G/dz3

	 }
	  
	 Bx_local = y_local*(G0-(1./12)*(3*x_local*x_local+y_local*y_local)*G2); 	// se (4) p1396 Incerti et.al
	 By_local = x_local*(G0-(1./12)*(3*y_local*y_local+x_local*x_local)*G2);	// se (5) p1396 Incerti et.al
	 Bz_local = x_local*y_local*(G1-(1./12)*(x_local*x_local+y_local*y_local)*G3);	// se (6) p1396 Incerti et.al

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
