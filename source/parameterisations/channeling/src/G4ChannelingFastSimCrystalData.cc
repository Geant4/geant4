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

#include "G4ChannelingFastSimCrystalData.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4ChannelingFastSimCrystalData::G4ChannelingFastSimCrystalData()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingFastSimCrystalData::SetMaterialProperties(const G4Material *crystal,
							   const G4String &lattice)
{
    G4String filename=crystal->GetName(); //input file
    filename.erase(0,3);

    if (fVerbosity)
      {
	G4cout << 
	  "======================================================================="  
	       << G4endl;
	G4cout << 
	  "======                 Crystal lattice data                    ========"  
	       << G4endl;
	G4cout << 
	  "======================================================================="  
	       << G4endl;    
	G4cout << "Crystal material: " << filename << G4endl;
      }

    //choice between planes (1D model) and axes (2D model)
    if (lattice.compare(0,1,"(")==0)
    {
      iModel=1; //planes
      filename = filename + "_planes_"; //temporary name
      if (fVerbosity)
	G4cout << "Crystal planes: " << lattice << G4endl;
    }
    else if (lattice.compare(0,1,"<")==0)
    {
      iModel=2; //axes
      filename = filename + "_axes_"; //temporary name
      if (fVerbosity)
	G4cout << "Crystal axes: " << lattice << G4endl;
    }

    //input file:
    filename = filename + lattice.substr(1,(lattice.length())-2) + ".dat";

    fNelements=(G4int)crystal->GetNumberOfElements();
    for(G4int i=0; i<fNelements; i++)
    {
      fZ1.push_back(crystal->GetElement(i)->GetZ());
      fAN.push_back(crystal->GetElement(i)->GetAtomicMassAmu());
      fI0.push_back(crystal->GetElement(i)->GetIonisation()->GetMeanExcitationEnergy());
    }

    G4double var;//just variable
    G4double unitIF;//unit of interpolation function

    std::ifstream vfilein;
    vfilein.open(filename);

    //read nuclear concentration
    for(G4int i=0; i<fNelements; i++)
    {
      vfilein >> var;
      fN0.push_back(var/cm3);
    }

    //read amplitude of thermal oscillations
    for(G4int i=0; i<fNelements; i++)
    {
      vfilein >> var;
      fU1.push_back(var*cm);
    }

        if (iModel==1)
    {
      //  read channel dimensions
      vfilein >> fDx;
      fDx*=cm;
      //  read interpolation step size
      vfilein >> fNpointsx;

      fDy = fDx;
      fNpointsy = 0;
    }
    else if (iModel==2)
    {
      //  read channel dimensions
      vfilein >> fDx >> fDy;
      fDx*=cm;
      fDy*=cm;
      //  read the number of nodes of interpolation
      vfilein >> fNpointsx >> fNpointsy;
    }

    //read the height of the potential well, necessary only for step length calculation
    vfilein >> fVmax;
    fVmax*=eV;
    fVmax2=2.*fVmax;

    //read the on-zero minimal potential inside the crystal,
    //necessary for angle recalculation for entrance/exit through
    //the crystal lateral surface
    vfilein >> fVMinCrystal;
    fVMinCrystal*=eV;

    // to create the class of interpolation for any function
    fElectricFieldX =
            new G4ChannelingFastSimInterpolation(fDx,fDy,fNpointsx,fNpointsy,iModel);
    if(iModel==2) {fElectricFieldY =
            new G4ChannelingFastSimInterpolation(fDx,fDy,fNpointsx,fNpointsy,iModel);}
    fElectronDensity =
            new G4ChannelingFastSimInterpolation(fDx,fDy,fNpointsx,fNpointsy,iModel);
    fMinIonizationEnergy =
            new G4ChannelingFastSimInterpolation(fDx,fDy,fNpointsx,fNpointsy,iModel);

    // do it for any element of crystal material
    for(G4int i=0; i<fNelements; i++)
    {
        fNucleiDensity.push_back(
            new G4ChannelingFastSimInterpolation(fDx,fDy,fNpointsx,fNpointsy,iModel));
    }

    if (iModel==1)
    {
      G4double ai, bi, ci, di;
      for(G4int i=0; i<fNpointsx; i++)
      {
        //reading the coefficients of cubic spline
        vfilein >> ai >> bi >> ci >> di;
        //setting spline coefficients for electric field
        unitIF=eV/cm;
        fElectricFieldX->SetCoefficients1D(ai*unitIF, bi*unitIF,
                                           ci*unitIF, di*unitIF, i);

        //reading the coefficients of cubic spline
        vfilein >> ai >> bi >> ci >> di;
        //setting spline coefficients for nuclear density (first element)
        unitIF=1.;
        fNucleiDensity[0]->SetCoefficients1D(ai*unitIF, bi*unitIF,
                                             ci*unitIF, di*unitIF, i);

        //reading the coefficients of cubic spline
        vfilein >> ai >> bi >> ci >> di;
        //setting spline coefficients for electron density
        unitIF=1./cm3;
        fElectronDensity->SetCoefficients1D(ai*unitIF, bi*unitIF,
                                            ci*unitIF, di*unitIF, i);

        //reading the coefficients of cubic spline
        vfilein >> ai >> bi >> ci >> di;
        //setting spline coefficients for minimal ionization energy
        unitIF=eV;
        fMinIonizationEnergy->SetCoefficients1D(ai*unitIF, bi*unitIF,
                                                ci*unitIF, di*unitIF, i);

        for(G4int ii=1; ii<fNelements; ii++)
        {
            //reading the coefficients of cubic spline
            vfilein >> ai >> bi >> ci >> di;
            //setting spline coefficients for nuclear density (other elements if any)
            unitIF=1.;
            fNucleiDensity[ii]->SetCoefficients1D(ai*unitIF, bi*unitIF,
                                                  ci*unitIF, di*unitIF, i);
        }
      }
    }
    else if (iModel==2)
    {
      G4double ai3D, bi3D, ci3D;
      for(G4int j=0; j<fNpointsy; j++)
      {
        for(G4int i=0; i<fNpointsx+1; i++)
        {
          for(G4int k=0; k<2; k++)
          {
            //reading the coefficients of cubic spline
            vfilein >> ai3D >> bi3D >> ci3D;
            unitIF=eV;
            //setting spline coefficients for minimal ionization energy
            fMinIonizationEnergy->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF, ci3D*unitIF,
                                                    i, j, k);

            //reading the coefficients of cubic spline
            vfilein >> ai3D >> bi3D >> ci3D;
            //setting spline coefficients for horizontal electric field
            unitIF=eV/cm;
            fElectricFieldX->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF, ci3D*unitIF,
                                               i, j, k);

            //reading the coefficients of cubic spline
            vfilein >> ai3D >> bi3D >> ci3D;
            //setting spline coefficients for vertical electric field
            unitIF=eV/cm;
            fElectricFieldY->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF, ci3D*unitIF,
                                               i, j, k);

            //reading the coefficients of cubic spline
            vfilein >> ai3D >> bi3D >> ci3D;
            //setting spline coefficients for nuclear density (first element)
            unitIF=1.;
            fNucleiDensity[0]->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF, ci3D*unitIF,
                                                 i, j, k);

            //reading the coefficients of cubic spline
            vfilein >> ai3D >> bi3D >> ci3D;
            //setting spline coefficients for electron density
            unitIF=1./cm3;
            fElectronDensity->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF, ci3D*unitIF,
                                                i, j, k);

            for(G4int ii=1; ii<fNelements; ii++)
            {
                //reading the coefficients of cubic spline
                vfilein >> ai3D >> bi3D >> ci3D;
                //setting spline coefficients for nuclear density (other elements if any)
                unitIF=1.;
                fNucleiDensity[ii]->SetCoefficients2D(ai3D*unitIF, bi3D*unitIF,
                                                      ci3D*unitIF,
                                                      i, j, k);
            }

          }
        }
      }
    }

    vfilein.close();

    //set special values and coefficients
    G4double alphahbarc2=std::pow(CLHEP::fine_structure_const*CLHEP::hbarc ,2.);
    fK30=2.*CLHEP::pi*alphahbarc2/CLHEP::electron_mass_c2;

    for(G4int i=0; i<fNelements; i++)
    {
      fRF.push_back((std::pow(9*CLHEP::pi*CLHEP::pi/128/fZ1[i],1/3.))
                    *0.5291772109217*angstrom);//Thomas-Fermi screening radius

      fTetamax0.push_back(CLHEP::hbarc/(fR0*std::pow(fAN[i],1./3.)));
      fTeta10.push_back(CLHEP::hbarc/fRF[i]);
      fPu11.push_back(std::pow(fU1[i]/CLHEP::hbarc,2.));

      fK20.push_back(alphahbarc2*4*CLHEP::pi*fN0[i]*fZ1[i]*fZ1[i]);

      fK40.push_back(3.76*std::pow(CLHEP::fine_structure_const*fZ1[i],2.));

      fKD.push_back(fK30*fZ1[i]*fN0[i]);
    }

    fBB.resize(fNelements);
    fE1XBbb.resize(fNelements);
    fBBDEXP.resize(fNelements);
    fPzu11.resize(fNelements);
    fTeta12.resize(fNelements);
    fTetamax2.resize(fNelements);
    fTetamax12.resize(fNelements);
    fK2.resize(fNelements);

    fChangeStep = CLHEP::pi*std::min(fDx,fDy)/fNsteps;//necessary to define simulation
                                         //step = fChannelingStep =
                                         // = fChannelingStep0*sqrt(pv)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4ChannelingFastSimCrystalData::CoordinatesFromBoxToLattice
                                              (const G4ThreeVector &pos0)
{
   G4double x0=pos0.x(),y0=pos0.y(),z0=pos0.z();
   z0+=fHalfDimBoundingBox.z();
   G4double x,y,z;

   if (fBent)
   {
       // for bent crystal
       G4double rsqrt = std::sqrt(fBendingRsquare -
                                  fBending2R*(x0*fCosMiscutAngle - z0*fSinMiscutAngle) +
                                  x0*x0 + z0*z0);
       //transform to co-rotating reference system connected with "central plane/axis"
       x = fBendingR - rsqrt;
       y = y0;
       z = fBendingR*std::asin((z0*fCosMiscutAngle + x0*fSinMiscutAngle)/rsqrt);
   }
   else
   {
       //for straight crystal
       x = x0*fCosMiscutAngle - z0*fSinMiscutAngle;
       y = y0;
       z = x0*fSinMiscutAngle + z0*fCosMiscutAngle;

       //for crystalline undulator
       if(fCU){x-=GetCUx(z);}
   }

   //calculation of coordinates within a channel (periodic cell)
   fNChannelx=std::floor(x/fDx); //remember the horizontal channel number
                                 //to track the particle
   x-=fNChannelx*fDx;
   fNChannely=std::floor(y/fDy);//remember the vertical channel number
                                //to track the particle (=0 for planar case)
   y-=fNChannely*fDy;
   //correction of the longitudinal coordinate
   if (fBent) {fCorrectionZ = fBendingR/(fBendingR-fNChannelx*fDx);}

   return G4ThreeVector(x,y,z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4ChannelingFastSimCrystalData::CoordinatesFromLatticeToBox(
        const G4ThreeVector &pos)
{
   G4double x=pos.x(),y=pos.y(),z=pos.z();

   //transform to co-rotating reference system connected with "central plane/axis"
   x+=fNChannelx*fDx;
   y+=fNChannely*fDy;

   G4double x0,y0,z0;

   if (fBent)
   {
       // for bent crystal
       G4double rcos = (fBendingR - x)*(1. - std::cos(z/fBendingR));
       G4double a = x + rcos;
       G4double b = std::sqrt(x*x + fBending2R*rcos - a*a);

       //transform to Box coordinates
       x0 = a*fCosMiscutAngle + b*fSinMiscutAngle;
       y0 = y;
       z0 = b*fCosMiscutAngle - a*fSinMiscutAngle;
   }
   else
   {
       //for crystalline undulator
       if(fCU){x+=GetCUx(z);}

       //for straight crystal
       x0 = x*fCosMiscutAngle + z*fSinMiscutAngle;
       y0 = y;
       z0 =-x*fSinMiscutAngle + z*fCosMiscutAngle;
   }

   return G4ThreeVector(x0,y0,z0-fHalfDimBoundingBox.z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4ChannelingFastSimCrystalData::ChannelChange(G4double& x,
							    G4double& y,
							    G4double& z)
{

    //test of enter in other channel
    if (x<0)
    {
        fNChannelx-=1;
        x+=fDx; //enter in other channel
        //correction of the longitudinal coordinate
        if (fBent) {fCorrectionZ = fBendingR/(fBendingR-fNChannelx*fDx);}
    }
    else if (x>=fDx)
    {
        fNChannelx+=1;
        x-=fDx; //enter in other channel
        //correction of the longitudinal coordinate
        if (fBent) {fCorrectionZ = fBendingR/(fBendingR-fNChannelx*fDx);}
    }

    //test of enter in other channel
    if (y<0)
    {
        fNChannely-=1;
        y+=fDy; //enter in other channel
    }
    else if (y>=fDy)
    {
        fNChannely+=1;
        y-=fDy; //enter in other channel
    }

    return G4ThreeVector(x,y,z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
