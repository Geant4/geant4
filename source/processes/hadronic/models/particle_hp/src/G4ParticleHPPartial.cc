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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 13-Jan-06 fix in Sample by T. Koi
//
#include "G4ParticleHPPartial.hh"  
#include "G4ParticleHPInterpolator.hh"
#include "Randomize.hh"

G4ParticleHPVector * G4ParticleHPPartial::GetY(G4double e1)
  {
    G4ParticleHPVector * aBuffer = new G4ParticleHPVector();
    G4int i;
    if(nData==1)
    {
      for(i=0; i<data[0].GetVectorLength(); i++)
      {
        aBuffer->SetInterpolationManager(data[0].GetInterpolationManager());
        aBuffer->SetData(i , data[0].GetX(i), data[0].GetY(i));
      }
      return aBuffer;
    }
    for (i=0; i<nData; i++)
    {
       if(X[i]>e1) break;
    }
    if(i==nData) i--;
    if(0==i) i=1;
    G4double x1,x2,y1,y2,y;
    G4int i1=0, ib=0;
    G4double E1 = X[i-1];
    G4double E2 = X[i];
    for(G4int ii=0; ii<data[i].GetVectorLength(); ii++)
    {
      x1 = data[i-1].GetX(std::min(i1, data[i-1].GetVectorLength()-1));
      x2 = data[i].GetX(ii);
      if(x1<x2&&i1<data[i-1].GetVectorLength())
      {
        y1 = data[i-1].GetY(i1);
        y2 = data[i].GetY(x1);
        if(E2-E1!=0)
        {
          y = theInt.Interpolate(theManager.GetScheme(i), e1, E1, E2, y1, y2);
        }
        else
        {
          y = 0.5*(y1+y2);
        }
        aBuffer->SetData(ib, x1, y);
        aBuffer->SetScheme(ib++, data[i-1].GetScheme(i1));
        i1++;
        if(x2-x1>0.001*x1)
        {
          ii--;
        }
      }
      else
      {
        y1 = data[i-1].GetY(x2);
        y2 = data[i].GetY(ii);
        if(E2-E1!=0)
        {
          y = theInt.Interpolate(theManager.GetScheme(i), e1, E1, E2, y1, y2);
        }
        else
        {
          y = 0.5*(y1+y2);
        }       
        aBuffer->SetData(ib, x2, y);
        aBuffer->SetScheme(ib++, data[i].GetScheme(ii));
        if(x1-x2<0.001*x2) i1++;
      }
    }
    return aBuffer;
  }
  
  G4double G4ParticleHPPartial::Sample(G4double x)
  {
    G4int i;
    for (i=0; i<nData; i++)
    {
      if(x<X[i]) break;
    }
    G4ParticleHPVector theBuff;
    if(i==0) 
    {
      theBuff.SetInterpolationManager(data[0].GetInterpolationManager());
      for(G4int ii=0;ii<GetNEntries(0);ii++) 
      {
        theBuff.SetX(ii, GetX(0,ii));
        theBuff.SetY(ii, GetY(0,ii));
      }
    }
    //else if(i==nData-1) this line will be delete 
    else if ( i == nData )  
    {
      for(i=0;i<GetNEntries(nData-1);i++) 
      {
        theBuff.SetX(i, GetX(nData-1,i));
        theBuff.SetY(i, GetY(nData-1,i));      
        theBuff.SetInterpolationManager(data[nData-1].GetInterpolationManager());
      }
    }
    else
    {
      G4int low = i-1;
      G4int high = low+1;
      G4double x1,x2,y1,y2;
      G4int i1=0, i2=0, ii=0;
      x1 = X[low];
      x2 = X[high];
      while(i1<GetNEntries(low)||i2<GetNEntries(high))
      {
	if(   (GetX(low,i1)<GetX(high,i2) && i1<GetNEntries(low))
            ||(i2==GetNEntries(high)) )
	{
          theBuff.SetX(ii, GetX(low,i1));
          y1 = GetY(low,i1);
          y2 = GetY(high, GetX(low,i1));  //prob at ident theta
          theBuff.SetY(ii, theInt.Interpolate(theManager.GetScheme(high),
                                              x, x1, x2, y1, y2)); //energy interpol
          theBuff.SetScheme(ii, data[low].GetScheme(i1));
          if(std::abs(GetX(low,i1)-GetX(high,i2))<0.001) i2++;
          i1++;
          ii++;
	}
	else
	{
          theBuff.SetX(ii, GetX(high,i2));
	  //*******************************************************************************
	  //Change by E.Mendoza and D.Cano (CIEMAT):
          //y1 = GetY(high,i2);
          //y2 = GetY(low, GetX(high,i2));  //prob at ident theta
          y2 = GetY(high,i2);
          y1 = GetY(low, GetX(high,i2));  //prob at ident theta
	  //*******************************************************************************
          theBuff.SetY(ii, theInt.Interpolate(theManager.GetScheme(high),
                                              x, x1, x2, y1, y2)); //energy interpol
          theBuff.SetScheme(ii, data[high].GetScheme(i2));
          if(std::abs(GetX(low,i1)-GetX(high,i2))<0.001) i1++;
          i2++;
          ii++;
	} 
      }
    }
    //buff is full, now sample.
    return theBuff.Sample();
  }
