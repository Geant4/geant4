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

#include "HadrontherapyElectricTabulatedField3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

namespace{  G4Mutex MyHadrontherapyLockEField=G4MUTEX_INITIALIZER;  }

HadrontherapyElectricTabulatedField3D::HadrontherapyElectricTabulatedField3D( const char* filename, G4double exOffset, G4double eyOffset, G4double ezOffset)
  :feXoffset(exOffset),feYoffset(eyOffset),feZoffset(ezOffset),einvertX(false),einvertY(false),einvertZ(false)
{
   //The format file is: X Y Z Ex Ey Ez

  G4double ElenUnit= cm;
  G4double EfieldUnit= volt/m;
  G4cout << "\n-----------------------------------------------------------"
	     << "\n      Electric field"
         << "\n-----------------------------------------------------------";

  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl;
  G4AutoLock lock(&MyHadrontherapyLockEField);

  ifstream file( filename ); // Open the file for reading.

  // Ignore first blank line
  char ebuffer[256];
  file.getline(ebuffer,256);

  // Read table dimensions
  file >> Enx >> Eny >> Enz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: "
	 << Enx << " " << Eny << " " << Enz << " ] "
	 << endl;

  // Set up storage space for table
  xEField.resize( Enx );
  yEField.resize( Enx );
  zEField.resize( Enx );
  G4int ix, iy, iz;
  for (ix=0; ix<Enx; ix++) {
    xEField[ix].resize(Eny);
    yEField[ix].resize(Eny);
    zEField[ix].resize(Eny);
  for (iy=0; iy<Eny; iy++) {
      xEField[ix][iy].resize(Enz);
      yEField[ix][iy].resize(Enz);
      zEField[ix][iy].resize(Enz);
    }
  }

  // Read in the data
  G4double Exval=0.;
  G4double Eyval=0.;
  G4double Ezval=0.;
  G4double Ex=0.;
  G4double Ey=0.;
  G4double Ez=0.;
       for (iz=0; iz<Enz; iz++) {
         for (iy=0; iy<Eny; iy++) {
            for (ix=0; ix<Enx; ix++) {
        file >> Exval >> Eyval >> Ezval >> Ex >> Ey >> Ez;

        if ( ix==0 && iy==0 && iz==0 ) {
          Eminx = Exval * ElenUnit;
          Eminy = Eyval * ElenUnit;
          Eminz = Ezval * ElenUnit;
        }
        xEField[ix][iy][iz] = Ex * EfieldUnit;
        yEField[ix][iy][iz] = Ey * EfieldUnit;
        zEField[ix][iy][iz] = Ez * EfieldUnit;
      }
    }
  }
  file.close();
  lock.unlock();

  Emaxx = Exval * ElenUnit;
  Emaxy = Eyval * ElenUnit;
  Emaxz = Ezval * ElenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl;
  G4cout << " ---> assumed the order:  x, y, z, Ex, Ey, Ez "
	 << "\n ---> Min values x,y,z: "
	 << Eminx/cm << " " << Eminy/cm << " " << Eminz/cm << " cm "
	 << "\n ---> Max values x,y,z: "
	 << Emaxx/cm << " " << Emaxy/cm << " " << Emaxz/cm << " cm "
	 << "\n ---> The field will be offset in x by " << exOffset/cm << " cm "
         << "\n ---> The field will be offset in y by " << eyOffset/cm << " cm "
         << "\n ---> The field will be offset in z by " << ezOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (Emaxx < Eminx) {swap(Emaxx,Eminx); einvertX = true;}
  if (Emaxy < Eminy) {swap(Emaxy,Eminy); einvertY = true;}
  if (Emaxz < Eminz) {swap(Emaxz,Eminz); einvertZ = true;}
  G4cout << "\nAfter reordering if neccesary"
	 << "\n ---> Min values x,y,z: "
	 << Eminx/cm << " " << Eminy/cm << " " << Eminz/cm << " cm "
	 << " \n ---> Max values x,y,z: "
	 << Emaxx/cm << " " << Emaxy/cm << " " << Emaxz/cm << " cm ";

  dx1 = Emaxx - Eminx;
  dy1 = Emaxy - Eminy;
  dz1 = Emaxz - Eminz;
  G4cout << "\n ---> Dif values x,y,z (range): "
	 << dx1/cm << " " << dy1/cm << " " << dz1/cm << " cm  "
	 << "\n-----------------------------------------------------------" << endl;
}

void HadrontherapyElectricTabulatedField3D::GetFieldValue(const G4double Epoint[4],
				      G4double *Efield ) const
{
    G4double x1 = Epoint[0] + feXoffset;
    G4double y1 = Epoint[1] + feYoffset;
    G4double z1 = Epoint[2] + feZoffset;

    // Position of given point within region, normalized to the range
    // [0,1]
    G4double Exfraction = (x1 - Eminx) / dx1;
    G4double Eyfraction = (y1 - Eminy) / dy1;
    G4double Ezfraction = (z1 - Eminz) / dz1;

    if (einvertX) { Exfraction = 1 - Exfraction;}
    if (einvertY) { Eyfraction = 1 - Eyfraction;}
    if (einvertZ) { Ezfraction = 1 - Ezfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    G4double exdindex, eydindex, ezdindex;

    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    G4double exlocal = ( std::modf(Exfraction*(Enx-1), &exdindex));
    G4double eylocal = ( std::modf(Eyfraction*(Eny-1), &eydindex));
    G4double ezlocal = ( std::modf(Ezfraction*(Enz-1), &ezdindex));

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    G4int exindex = static_cast<G4int>(std::floor(exdindex));
    G4int eyindex = static_cast<G4int>(std::floor(eydindex));
    G4int ezindex = static_cast<G4int>(std::floor(ezdindex));

    if ((exindex < 0) || (exindex >= Enx - 1) ||
        (eyindex < 0) || (eyindex >= Eny - 1) ||
        (ezindex < 0) || (ezindex >= Enz - 1))
    {
        Efield[0] = 0.0;
        Efield[1] = 0.0;
        Efield[2] = 0.0;
        Efield[3] = 0.0;
        Efield[4] = 0.0;
        Efield[5] = 0.0;
    }
    else
    {

/*
#ifdef DEBUG_G4intERPOLATING_FIELD
    G4cout << "Local x,y,z: " << exlocal << " " << eylocal << " " << ezlocal << endl;
    G4cout << "Index x,y,z: " << exindex << " " << eyindex << " " << ezindex << endl;
    G4double valx0z0, mulx0z0, valx1z0, mulx1z0;
    G4double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[exindex  ][0][ezindex];  mulx0z0=  (1-exlocal) * (1-ezlocal);
    valx1z0= table[exindex+1][0][ezindex];  mulx1z0=   exlocal    * (1-ezlocal);
    valx0z1= table[exindex  ][0][ezindex+1]; mulx0z1= (1-exlocal) * ezlocal;
    valx1z1= table[exindex+1][0][ezindex+1]; mulx1z1=  exlocal    * ezlocal;
#endif
*/
        // Full 3-dimensional version

        Efield[0] = 0.0;
        Efield[1] = 0.0;
        Efield[2] = 0.0;

        Efield[3] =
          xEField[exindex  ][eyindex  ][ezindex  ] * (1-exlocal) * (1-eylocal) * (1-ezlocal) +
          xEField[exindex  ][eyindex  ][ezindex+1] * (1-exlocal) * (1-eylocal) *    ezlocal  +
          xEField[exindex  ][eyindex+1][ezindex  ] * (1-exlocal) *    eylocal  * (1-ezlocal) +
          xEField[exindex  ][eyindex+1][ezindex+1] * (1-exlocal) *    eylocal  *    ezlocal  +
          xEField[exindex+1][eyindex  ][ezindex  ] *    exlocal  * (1-eylocal) * (1-ezlocal) +
          xEField[exindex+1][eyindex  ][ezindex+1] *    exlocal  * (1-eylocal) *    ezlocal  +
          xEField[exindex+1][eyindex+1][ezindex  ] *    exlocal  *    eylocal  * (1-ezlocal) +
          xEField[exindex+1][eyindex+1][ezindex+1] *    exlocal  *    eylocal  *    ezlocal ;
        Efield[4] =
          yEField[exindex  ][eyindex  ][ezindex  ] * (1-exlocal) * (1-eylocal) * (1-ezlocal) +
          yEField[exindex  ][eyindex  ][ezindex+1] * (1-exlocal) * (1-eylocal) *    ezlocal  +
          yEField[exindex  ][eyindex+1][ezindex  ] * (1-exlocal) *    eylocal  * (1-ezlocal) +
          yEField[exindex  ][eyindex+1][ezindex+1] * (1-exlocal) *    eylocal  *    ezlocal  +
          yEField[exindex+1][eyindex  ][ezindex  ] *    exlocal  * (1-eylocal) * (1-ezlocal) +
          yEField[exindex+1][eyindex  ][ezindex+1] *    exlocal  * (1-eylocal) *    ezlocal  +
          yEField[exindex+1][eyindex+1][ezindex  ] *    exlocal  *    eylocal  * (1-ezlocal) +
          yEField[exindex+1][eyindex+1][ezindex+1] *    exlocal  *    eylocal  *    ezlocal ;
        Efield[5] =
          zEField[exindex  ][eyindex  ][ezindex  ] * (1-exlocal) * (1-eylocal) * (1-ezlocal) +
          zEField[exindex  ][eyindex  ][ezindex+1] * (1-exlocal) * (1-eylocal) *    ezlocal  +
          zEField[exindex  ][eyindex+1][ezindex  ] * (1-exlocal) *    eylocal  * (1-ezlocal) +
          zEField[exindex  ][eyindex+1][ezindex+1] * (1-exlocal) *    eylocal  *    ezlocal  +
          zEField[exindex+1][eyindex  ][ezindex  ] *    exlocal  * (1-eylocal) * (1-ezlocal) +
          zEField[exindex+1][eyindex  ][ezindex+1] *    exlocal  * (1-eylocal) *    ezlocal  +
          zEField[exindex+1][eyindex+1][ezindex  ] *    exlocal  *    eylocal  * (1-ezlocal) +
          zEField[exindex+1][eyindex+1][ezindex+1] *    exlocal  *    eylocal  *    ezlocal ;
  }
//G4cout << "Getting electric field " << Efield[3]/(volt/m) << " " << Efield[4]/(volt/m) << " " << Efield[5]/(volt/m) << endl;
//G4cout << "For coordinates: " << Epoint[0] << " " << Epoint[1] << " " << Epoint[2] << G4endl;

/*std::ofstream WriteDataIn("ElectricFieldFC.out", std::ios::app);
       WriteDataIn	<<   Epoint[0]             << '\t' << "   "
			<<   Epoint[1]           << '\t' << "   "
			<<   Epoint[2]            << '\t' << "   "
			<<   Efield[3]/(volt/m)            << '\t' << "   "
			<<   Efield[4]/(volt/m)           << '\t' << "   "
			<<   Efield[5]/(volt/m)           << '\t' << "   "
			<< G4endl;    */
}
