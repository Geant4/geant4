// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPList.hh"

  G4NeutronHPList::G4NeutronHPList()
  {
    theData = new G4double[100]; 
    nPoints=100;
    nEntries=0;
  }
  
  G4NeutronHPList::~G4NeutronHPList()
  {
    delete [] theData;
  }

  void G4NeutronHPList::Check(G4int i)
  {
    if(i<0) 
    {
      G4int dummy; cin >> dummy;
    }
    if(i>nEntries) G4Exception("Skipped some index numbers in G4NeutronHPList");
    if(i==nPoints)
    {
      nPoints += 50;
      G4double * buff = new G4double[nPoints];
      for (G4int j=0; j<nEntries; j++) buff[j] = theData[j];
      delete [] theData;
      theData = buff;
    }
    if(i==nEntries) nEntries++;
  }
  
  void G4NeutronHPList::Init(ifstream & aDataFile, G4int nPar, G4double unit)
  {
    G4int i;
    G4double y;
    for (i=0; i<nPar; i++)
    {
      aDataFile >> y;
      SetValue(i,y*unit);
    }
  }

  void G4NeutronHPList::Init(ifstream & aDataFile, G4double unit)
  {
    G4int total, i;
    aDataFile >> total;
    G4double y;
    for (i=0;i<total;i++)
    {
      aDataFile >>y;
      SetValue(i,y*unit);
    }
  }
  
  G4double G4NeutronHPList::GetValue(G4int i) 
  { 
//    G4cout << "TestList "<<i<<" "<<nEntries<<endl;
    if(nEntries<0)
    {
//      G4cout <<nPoints<<" "<<nEntries<<" "<<theData<<endl;
//      for(G4int ii=0; ii<2; ii++) G4cout << theData[ii]<<" ";
//      G4cout << endl;
    }
    if (i<0) i=0;
    if(i>=GetListLength()) i=GetListLength()-1;
    return theData[i];
  }
  
