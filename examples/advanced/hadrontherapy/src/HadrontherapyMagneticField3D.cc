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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyMagneticField3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

namespace{  G4Mutex MyHadrontherapyLock=G4MUTEX_INITIALIZER;  }

HadrontherapyMagneticField3D::HadrontherapyMagneticField3D( const char* filename, double xOffset )
  :fXoffset(xOffset),invertX(false),invertY(false),invertZ(false)
{
   //The format file is: X Y Z Ex Ey Ez

  double lenUnit= meter;
  double fieldUnit= tesla;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";


  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl;
  G4AutoLock lock(&MyHadrontherapyLock);

  ifstream file( filename ); // Open the file for reading.

  // Ignore first blank line
  char buffer[256];
  file.getline(buffer,256);

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
  G4double xval=0.;
  G4double yval=0.;
  G4double zval=0.;
  G4double bx=0.;
  G4double by=0.;
  G4double bz=0.;
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
        file >> xval >> yval >> zval >> bx >> by >> bz ;
        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
      }
    }
  }
  file.close();

  lock.unlock();

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl;
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: "
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: "
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by " << xOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;}
  if (maxy < miny) {swap(maxy,miny); invertY = true;}
  if (maxz < minz) {swap(maxz,minz); invertZ = true;}
  G4cout << "\nAfter reordering if neccesary"
	 << "\n ---> Min values x,y,z: "
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: "
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): "
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;
}

void HadrontherapyMagneticField3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{
    double x = point[0]+ fXoffset;
    double y = point[1];
    double z = point[2];

    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy;
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

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
    int xindex = static_cast<int>(std::floor(xdindex));
    int yindex = static_cast<int>(std::floor(ydindex));
    int zindex = static_cast<int>(std::floor(zdindex));

      // Check that the point is within the defined region
    if ((xindex < 0) || (xindex >= nx - 1) ||
        (yindex < 0) || (yindex >= ny - 1) ||
        (zindex < 0) || (zindex >= nz - 1))
    {
        Bfield[0] = 0.0;
        Bfield[1] = 0.0;
        Bfield[2] = 0.0;
    }
    else
    {

#ifdef DEBUG_INTERPOLATING_FIELD
        G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
        G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
        double valx0z0, mulx0z0, valx1z0, mulx1z0;
        double valx0z1, mulx0z1, valx1z1, mulx1z1;
        valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
        valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
        valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
        valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
        Bfield[0] =
          xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
          xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
          xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
          xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
          xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
          xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
          xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
          xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

        Bfield[1] =
          yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
          yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
          yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
          yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
          yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
          yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
          yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
          yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

        Bfield[2] =
          zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
          zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
          zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
          zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
          zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
          zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
          zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
          zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    }

//In order to obtain the output file with the magnetic components read from a particle passing in the magnetic field
/*	std::ofstream MagneticField("MagneticField.out", std::ios::app);
	   MagneticField<<   Bfield[0] << '\t' << "   "
			<<   Bfield[1] << '\t' << "    "
			<<   Bfield[2] << '\t' << "   "
			<<   point[0] << '\t' << "   "
			<<   point[1] << '\t' << "    "
			<<   point[2] << '\t' << "   "
			<< G4endl;*/

}
