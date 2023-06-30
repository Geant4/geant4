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

/// \file G4ChannelingFastSimInterpolation.cc
/// \brief Implementation of the G4ChannelingFastSimInterpolation class

#include "G4ChannelingFastSimInterpolation.hh"
#include "G4SystemOfUnits.hh"

G4ChannelingFastSimInterpolation::G4ChannelingFastSimInterpolation(G4double dx0,
                                                               G4double dy0,
                                                               G4int nPointsx0,
                                                               G4int nPointsy0,
                                                               G4int iModel0)
{
    fDx = dx0;
    fDy = dy0;
    nPointsx = nPointsx0;
    nPointsy = nPointsy0;
    fStepi = fDx/nPointsx;
    fStepi2= fStepi*fStepi;
    iModel = iModel0;

    //resize vectors for interpolation coefficients
    if(iModel==1)
    {
        fAI.resize(nPointsx+1);
        fBI.resize(nPointsx+1);
        fCI.resize(nPointsx+1);
        fDI.resize(nPointsx+1);
    }
    else if(iModel==2)
    {
        fStepj = fDy/nPointsy;
        fAI3D.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));
        fBI3D.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));
        fCI3D.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));
        fAI3D3.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));
        fBI3D3.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));
        fCI3D3.resize(nPointsx+1, std::vector<G4double>(nPointsy+1));

        for(G4int i=0; i<=nPointsx; i++)
        {
            fCI3D[i][nPointsy]=0.;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ChannelingFastSimInterpolation::GetIF(G4double xx, G4double yy){
    G4double SplineF=0.;
    if(iModel==1)
    {
      SplineF = Spline1D(xx);
    }
    else if(iModel==2)
    {
      SplineF = Spline2D(xx, yy);
    }
    return SplineF;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingFastSimInterpolation::SetCoefficients1D(G4double AI0,
                                                       G4double BI0,
                                                       G4double CI0,
                                                       G4double DI0,
                                                       G4int i){
    fAI[i] = AI0;
    fBI[i] = BI0/cm;
    fCI[i] = CI0/cm2;
    fDI[i] = DI0/cm3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingFastSimInterpolation::SetCoefficients2D(G4double AI3D0,
                                                       G4double BI3D0,
                                                       G4double CI3D0,
                                                       G4int i,
                                                       G4int j,
                                                       G4int k){
    if (k==0)
    {
        fAI3D[i][j] = AI3D0/fStepj/fStepi/6.;
        fBI3D[i][j] = BI3D0/fStepj/fStepi/6.;
        fCI3D[i][j] = CI3D0/fStepj/fStepi/6./cm2;
    }
    else if (k==1)
    {
        fAI3D3[i][j] = AI3D0/fStepj/fStepi/6./cm2;
        fBI3D3[i][j] = BI3D0/fStepj/fStepi/6./cm2;
        fCI3D3[i][j] = CI3D0/fStepj/fStepi/6./cm2/cm2;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ChannelingFastSimInterpolation::Spline1D(G4double xx) //cubic spline of
                                                               //1-variable function
{
     G4double x1 = xx;

     //if a particle escapes the interpolation area
     if (x1<0.)
     {
         x1 += fDx;
     }
     else if (x1>=fDx)
     {
         x1 -= fDx;
     }

//   calculation of interpolation function
     G4int mmx = std::floor(x1/fStepi);
     x1 -= (mmx+1)*fStepi;
     G4double Spline3 = fAI[mmx]+x1*(fBI[mmx]+x1*(fCI[mmx]+fDI[mmx]*x1));
     return Spline3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   G4double G4ChannelingFastSimInterpolation::Spline2D(G4double xx, G4double yy)//cubic
                                                      //spline of 2-variable function
   {
     G4double x1 = xx;
     G4double y1 = yy;

     //if a particle escapes the interpolation area
     if (x1<0.)
     {
         x1 += fDx;
     }
     else if (x1>=fDx)
     {
         x1 -= fDx;
     }
     if (y1<0.)
     {
         y1 += fDy;
     }
     else if (y1>=fDy)
     {
         y1 -= fDy;
     }

     //  calculation of interpolation function
     G4int mmx = std::floor(x1/fStepi);
     G4int mmy = std::floor(y1/fStepj);

     G4double tt1 = x1-mmx*fStepi;
     G4double tt2 = fStepi-tt1;
     G4double tt13 = tt1*tt1*tt1;
     G4double tt23 = tt2*tt2*tt2;

     G4double tt1y = y1-mmy*fStepj;
     G4double tt2y = fStepj-tt1y;
     G4double tt1y3 = tt1y*tt1y*tt1y;
     G4double tt2y3 = tt2y*tt2y*tt2y;

     G4double spl3dxx1 = fCI3D3[mmx  ][mmy]*tt2y3 + fCI3D3[mmx  ][mmy+1]*tt1y3 +
                         fAI3D3[mmx  ][mmy]*tt2y  + fBI3D3[mmx  ][mmy]*tt1y;
     G4double spl3dxx2 = fCI3D3[mmx+1][mmy]*tt2y3 + fCI3D3[mmx+1][mmy+1]*tt1y3 +
                         fAI3D3[mmx+1][mmy]*tt2y  + fBI3D3[mmx+1][mmy]*tt1y;
     G4double spl3d1 =   fCI3D[ mmx  ][mmy]*tt2y3 +  fCI3D[mmx  ][mmy+1]*tt1y3 +
                         fAI3D[ mmx  ][mmy]*tt2y  +  fBI3D[mmx  ][mmy]*tt1y;
     G4double spl3d2 =   fCI3D[ mmx+1][mmy]*tt2y3 +  fCI3D[mmx+1][mmy+1]*tt1y3 +
                         fAI3D[ mmx+1][mmy]*tt2y  +  fBI3D[mmx+1][mmy]*tt1y;

     G4double spl3d = spl3dxx1*tt23 + spl3dxx2*tt13 +
             (spl3d1*6.-spl3dxx1*fStepi2)*tt2 + (spl3d2*6.-spl3dxx2*fStepi2)*tt1;

     return spl3d;
   }


