// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPVector.hh"
 
  // if the ranges do not match, constant extrapolation is used.
  G4NeutronHPVector & operator + (const G4NeutronHPVector & left, const G4NeutronHPVector & right)
  {
    G4NeutronHPVector * result = new G4NeutronHPVector;
    G4int j=0;
    G4double x;
    G4double yl, yr, y;
    G4int running = 0;
    for(G4int i=0; i<left.nEntries; i++)
    {
      while(j<right.nEntries)
      {
        if(right.GetX(j)<left.GetX(i)*1.001)
        {
          x = right.GetX(j);
          y = right.GetY(j)+left.GetY(x);
          result->SetData(running++, x, y);
          j++;
        }
        else if(abs((right.GetX(j)-left.GetX(i))/(left.GetX(i)+right.GetX(j)))>0.001)
        {
          x = left.GetX(i);
          y = left.GetY(i)+right.GetY(x);
          result->SetData(running++, x, y);
          break;
        }
        else
        {
          break;
        }
      }
      if(j==right.nEntries)
      {
        x = left.GetX(i);
        y = left.GetY(i)+right.GetY(x);
        result->SetData(running++, x, y);     
      }
    }
    result->ThinOut(0.02);
    return *result;
  }

  G4NeutronHPVector::G4NeutronHPVector()
  {
    theData = new G4NeutronHPDataPoint[100]; 
    nPoints=100;
    nEntries=0;
    Verbose=0;
    theIntegral=NULL;
    totalIntegral=-1;
  }
  
  G4NeutronHPVector::~G4NeutronHPVector()
  {
//    if(Verbose==1)G4cout <<"G4NeutronHPVector::~G4NeutronHPVector"<<endl;
    if(theData!=NULL)     
    {
      delete [] theData;
    }
//    if(Verbose==1)G4cout <<"Vector: delete theData"<<endl;
    if(theIntegral!=NULL) delete [] theIntegral;
//    if(Verbose==1)G4cout <<"Vector: delete theIntegral"<<endl;
  }
  
  G4NeutronHPVector & G4NeutronHPVector::  
  operator = (const G4NeutronHPVector & right)
  {
    if(&right == this) return *this;
    
    G4int i;
    
    totalIntegral = right.totalIntegral;
    if(right.theIntegral!=NULL) theIntegral = new G4double[right.nEntries];
    for(i=0; i<right.nEntries; i++)
    {
      SetPoint(i, right.GetPoint(i)); // copy theData
      if(right.theIntegral!=NULL) theIntegral[i] = right.theIntegral[i];
    }
    theManager = right.theManager; 
    label = right.label;
  
    Verbose = right.Verbose;
    return *this;
  }

  G4double G4NeutronHPVector::GetXsec(G4double e) const
  {
    if(nEntries <= 1) 
    {
      if(nEntries == 0) return 0;
      return theData[0].GetY();
    }
    G4int found = 0;
    G4int low   = 0;
    G4int high  = 0;
    G4double eps = 0.00001*e;  // fast fix of precision problems, needs improvement @@@
//    if(Verbose==1) G4cout <<"G4NeutronHPVector::GetXsec";
    if(e<=theData[0].GetX()) return theData[0].GetY();
    G4int i=0, ii;
    for (ii=0; ii<nEntries/10+1; ii++) // von null weg, weil sonst <10 im argen liegt.
    {
      i = ii;
      if(theData[10*i].GetX()+eps>e) break;
    }
//    if(Verbose==1) G4cout << low<<" "<<high<<" "<<i<<" "<<nEntries<<" ";
    if(i!=(nEntries/10))
    {
      i=10*i;
      for (G4int j=0; j<11; j++)
      {
        if(theData[i].GetX()<e+eps) break;
        i--;
      }
      if(i>nEntries-2) i = nEntries-2;
      low = i;
      high = i+1;
    }
    else
    {
      i=max(0,10*(i-1));
      while (i<nEntries)
      {
        if(theData[i].GetX()>e) break;
        i++;
      } 
      if(i>nEntries-1) i = nEntries-1;
      low  = i-1;
      high = i;
    }
//    if(Verbose==1) G4cout << "sss"<<low<<" "<<high<<" ";
    G4double x1, x2, y1, y2, x, y;
    while ( theData[low].GetX()-e > 0.0000001*e ) 
    {
      low--;
      if(low<0) return theData[0].GetY();
    }
    while ( theData[high].GetX()-e < -0.0000001*e && high!=nEntries-1)
    {
      high++;
    }
    while( theData[high].GetX()-theData[low].GetX()<0.0000001*e)
    {
      if(high<nEntries-1)
      {
        high++;
      }
      else
      {
        low--;
        if(low<0) return theData[0].GetY();
      }
    }
//    if(Verbose==1) G4cout << "ddd"<<low<<" "<<high<<" ";
    x = e;
    x1 = theData[low] .GetX();
    x2 = theData[high].GetX();
    y1 = theData[low] .GetY();
    y2 = theData[high].GetY();
    y = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, y1, y2);
    if(e>=theData[nEntries-1].GetX()) 
    {
      y=theData[nEntries-1].GetY();
    }
    return y;
  }

  void G4NeutronHPVector::Dump()
  {
    G4cout << nEntries<<endl;
    for(G4int i=0; i<nEntries; i++)
    {
      G4cout << theData[i].GetX()<<" ";
      G4cout << theData[i].GetY()<<" ";
//      if (i!=1&&i==5*(i/5)) G4cout << endl;
      G4cout << endl;
    }
    G4cout << endl;
  }
  
  void G4NeutronHPVector::Check(G4int i)
  {
//    G4cout << "1: i: "<<i<<" nEntries: "<<nEntries<<" nPoints: "<<nPoints<<endl;
    if(i>nEntries) G4Exception("Skipped some index numbers in G4NeutronHPVector");
    if(i==nPoints)
    {
      nPoints += 50;
//    G4cout << "2a: i: "<<i<<" nEntries: "<<nEntries<<" nPoints: "<<nPoints<<endl;
      G4NeutronHPDataPoint * buff = new G4NeutronHPDataPoint[nPoints];
//    G4cout << "2b: i: "<<i<<" nEntries: "<<nEntries<<" nPoints: "<<nPoints<<endl;
      if(nPoints!=50)
      {
//        G4cout << "copying 1: nEntries="<<nEntries<<" nPoints="<<nPoints<<endl;
        for (G4int j=0; j<nEntries; j++) buff[j] = theData[j];
//        G4cout << "copying 2"<<endl;
        delete [] theData;
//    G4cout << "3: i: "<<i<<" nEntries: "<<nEntries<<" nPoints: "<<nPoints<<endl;
      }
      theData = buff;
    }
    if(i==nEntries) nEntries=i+1;
//    G4cout << "4: i: "<<i<<" nEntries: "<<nEntries<<" nPoints: "<<nPoints<<endl;
  }

  void G4NeutronHPVector::
  Merge(G4InterpolationScheme aScheme, G4double aValue, 
        G4NeutronHPVector * active, G4NeutronHPVector * passive)
  { 
    // interpolate between labels according to aScheme, cut at aValue, 
    // continue in unknown areas by substraction of the last difference.
    
    CleanUp();
    G4int s = 0, n=0, i=0, m=0;
    G4bool flag;
    G4NeutronHPVector * tmp;
    G4int a = s, p = n, t;
    while ( a<active->GetVectorLength() )
    {
      if(active->GetEnergy(a) <= passive->GetEnergy(p))
      {
        G4double xa  = active->GetEnergy(a);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(a), passive->GetXsec(xa));
        SetData(m, xa, yy);
        theManager.AppendScheme(m, active->GetScheme(a));
        m++;
        a++;
        G4double xp = passive->GetEnergy(p);
        if( abs(abs(xp-xa)/xa)<0.0000001&&a<active->GetVectorLength() ) 
        {
          p++;
          tmp = active; t=a;
          active = passive; a=p;
          passive = tmp; p=t;
        }
      } else {
        tmp = active; t=a;
        active = passive; a=p;
        passive = tmp; p=t;
      }
    }
    
    G4double deltaX = passive->GetXsec(GetEnergy(m-1)) - GetXsec(m-1);
    while (p!=passive->GetVectorLength()&&passive->GetEnergy(p)<=aValue)
    {
      G4double anX;
      anX = passive->GetXsec(p)-deltaX;
      if(anX>0)
      {
        if(abs(GetEnergy(m-1)-passive->GetEnergy(p))/passive->GetEnergy(p)>0.0000001)
        {
          SetData(m, passive->GetEnergy(p), anX);
          theManager.AppendScheme(m++, passive->GetScheme(p));
        }
      }
      p++;
    }
  }
  
  void G4NeutronHPVector::ThinOut(G4double precision)
  {
    // make the new vector
    G4NeutronHPDataPoint * aBuff = new G4NeutronHPDataPoint[nPoints];
    G4double x, x1, x2, y, y1, y2;
    G4int count = 0, current = 2, start = 1;
    
    // First element always goes and is never tested.
    aBuff[0] = theData[0];
    
    // Find the rest
    while(current < GetVectorLength())
    {
      x1=aBuff[count].GetX();
      y1=aBuff[count].GetY();
      x2=theData[current].GetX();
      y2=theData[current].GetY();
      for(G4int j=start; j<current; j++)
      {
	x = theData[j].GetX();
	y = theInt.Lin(x, x1, x2, y1, y2);
	if (abs(y-theData[j].GetY())>precision*y)
	{
	  aBuff[++count] = theData[current-1]; // for this one, everything was fine
          start = current; // the next candidate
	  break;
	}
      }
      current++ ;
    }
    // The last one also always goes, and is never tested.
    aBuff[++count] = theData[GetVectorLength()-1];
    delete [] theData;
    theData = aBuff;
    nEntries = count+1;
  }
