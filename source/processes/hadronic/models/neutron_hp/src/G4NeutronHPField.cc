// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPField.hh"
#include "G4ios.hh"


  G4NeutronHPField::G4NeutronHPField()
  {
    theData = new G4NeutronHPFieldPoint[100]; 
    nPoints=100;
    nEntries=0;
    theData->SetData(0,0,0);
  }
  
  G4NeutronHPField::~G4NeutronHPField(){ delete [] theData;}
  
  G4double G4NeutronHPField::GetY(G4double e, G4int j)
  {
    G4int found = 0;
    G4int low   = 0;
    G4int high  = 0;
    G4int i;
    for (i=1; i<nEntries/10; i++)
    {
      if(theData[10*i].GetX()>e) break;
    }
    if(i==(nEntries/10))
    {
      i=10*i;
      while (i<nEntries)
      {
        if(theData[i++].GetX()>e) break;
      } 
      if (i==nEntries)
      {
        low  = nEntries-1;
        high = nEntries-2;
      }else{
        low = i-1;
        high = i;
      }
    }else{
      for (G4int j=0; j<10; j++)
      {
        if(theData[i].GetX()<e) break;
        i--;
      }
      low = i;
      high = i+1;
    }
    G4double x1, x2, y1, y2, x, y;
    x = e;
    x1 = theData[low] .GetX();
    x2 = theData[high].GetX();
    y1 = theData[low] .GetY(j);
    y2 = theData[high].GetY(j);
    y = x*(y2-y1)/(x2-x1);
    return y += y2-x2*(y2-y1)/(x2-x1);
  }

  void G4NeutronHPField::Dump()
  {
    G4cout << nEntries<<endl;
    for(G4int i=0; i<nEntries; i++)
    {
      G4cout << theData[i].GetX()<<" ";
      for(G4int j=0; j<theData[i].GetDepth(); j++)
      {
        G4cout << theData[i].GetY(j)<<" ";
      }
      G4cout << endl;
    }
  }
  
  void G4NeutronHPField::Check(G4int i)
  {
    if(i>nEntries) G4Exception("Skipped some index numbers in G4NeutronHPField");
    if(i==nPoints)
    {
      nPoints += 50;
      G4NeutronHPFieldPoint * buff = new G4NeutronHPFieldPoint[nPoints];
//      G4cout << "copying 1"<<endl;
      for (G4int j=0; j<nEntries; j++) 
      {
        buff[j] = theData[j];
      }
//      G4cout << "copying 2"<<endl;
      delete [] theData;
      theData = buff;
    }
    if(i==nEntries) nEntries=i+1;
  }
