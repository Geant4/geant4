// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPPartial.hh"  
#include "G4NeutronHPInterpolator.hh"
#include "Randomize.hh"

G4NeutronHPVector * G4NeutronHPPartial::GetY(G4double e1)
  {
    G4NeutronHPVector * aBuffer = new G4NeutronHPVector();
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
    G4double x1,x2,y1,y2,y, off, slope;
    G4int i1=0, ib=0;
    G4double E1 = X[i-1];
    G4double E2 = X[i];
    for(G4int ii=0; ii<data[i].GetVectorLength(); ii++)
    {
      x1 = data[i-1].GetX(min(i1, data[i-1].GetVectorLength()-1));
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
  
  G4double G4NeutronHPPartial::Sample(G4double x)
  {
    G4double result=0;
    G4int i;
    for (i=0; i<nData; i++)
    {
      if(x<X[i]) break;
    }
    G4NeutronHPVector theBuff;
    if(i==0) 
    {
      theBuff.SetInterpolationManager(data[0].GetInterpolationManager());
      for(G4int ii=0;ii<GetNEntries(0);i++) 
      {
        theBuff.SetX(ii, GetX(0,ii));
        theBuff.SetY(ii, GetY(0,ii));
      }
    }
    else if(i==nData-1)
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
          if(abs(GetX(low,i1)-GetX(high,i2))<0.001) i2++;
          i1++;
          ii++;
	}
	else
	{
          theBuff.SetX(ii, GetX(high,i2));
          y1 = GetY(high,i2);
          y2 = GetY(low, GetX(high,i2));  //prob at ident theta
          theBuff.SetY(ii, theInt.Interpolate(theManager.GetScheme(high),
                                              x, x1, x2, y1, y2)); //energy interpol
          theBuff.SetScheme(ii, data[high].GetScheme(i2));
          if(abs(GetX(low,i1)-GetX(high,i2))<0.001) i1++;
          i2++;
          ii++;
	} 
      }
    }
    //buff is full, now sample.
    return theBuff.Sample();
  }
