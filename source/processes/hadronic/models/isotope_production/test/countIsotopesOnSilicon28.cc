#include "g4std/fstream"
#include "g4std/iostream"
// #include "G4NeutronHPNames.hh"
// #include "G4SandiaTable.hh"
#include <stdlib.h>
#include "globals.hh"
// #include "G4NeutronHPDataUsed.hh"
// #include "G4NeutronHPVector.hh"

int main()
{
  int iiii=0;
  G4String theNames[11]={"isoinputs/Silicon/isoinput28_14_01.out",
                         "isoinputs/Silicon/isoinput28_14_10.out",
                         "isoinputs/Silicon/isoinput28_14_20.out",
                         "isoinputs/Silicon/isoinput28_14_30.out",  
                         "isoinputs/Silicon/isoinput28_14_40.out",
                         "isoinputs/Silicon/isoinput28_14_50.out",
                         "isoinputs/Silicon/isoinput28_14_60.out",  
                         "isoinputs/Silicon/isoinput28_14_70.out",
                         "isoinputs/Silicon/isoinput28_14_80.out",
                         "isoinputs/Silicon/isoinput28_14_90.out",
                         "isoinputs/Silicon/isoinput28_14_99.out"};

  double energies[11]={1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 99.};
  double xsec[11]={0.00, 1.913, 1.869,  1.869, 1.869, 1.869, 1.869, 1.869, 1.869, 1.869, 1.869};
  G4String * theList = new G4String[1000000];
  int ** theCounter = new int * [1000000];
  for(int i5=0; i5<1000000; i5++)
  {
    theCounter[i5] = new int[11];
  }
  int i,j,k;
  for(i=0; i<1000000; i++) for(j=0; j<11; j++){ theCounter[i][j]=0; theList[i]=""; }
  int counter=0; // number of iso types produced
  for(i=0; i<11; i++)
  {
    G4std::ifstream aDataSet(theNames[i], G4std::ios::in);
    for(j=0; j<1000000; j++) // 1000000 events each test-run 
    {
      G4String label;
      aDataSet >> label>>label;
      if(label=="UNCHANGED") continue;
      for(k=0; k<counter; k++)
      {
        if(label==theList[k]) break;
      }
      if(k==counter)
      {
        theList[counter]=label;
        counter++;
      }
      theCounter[k][i]++;
    }
  }
  cout << "ve/cr ee(11) r ";
  for(i=0; i<11; i++) cout << energies[i]<<" ";
  cout << G4endl;
  for(j=0; j<counter; j++)
  {
    cout << "ve/cr "<<theList[j]<<"(11) r ";
    for(i=0; i<11; i++)
    {
//      cout << energies[i]<<" "<<theCounter[j][i]/1000000.<<" ";
      cout <<xsec[i]*theCounter[j][i]/1000000.<<" ";
    }
    cout <<G4endl;
    cout << "ve/cr err"<<theList[j]<<"(11) r ";
    for(i=0; i<11; i++)
    {
//      cout << energies[i]<<" "<<theCounter[j][i]/1000000.<<" ";
      cout <<xsec[i]*sqrt(theCounter[j][i])/1000000.+0.1*xsec[i]*theCounter[j][i]/1000000.<<" ";
    }
    cout <<G4endl;
  }
  G4std::ifstream aDataSet("/data/NeutronProductionData/G4NDL0.2//IsotopeProduction/CrossSection/14_28_Silicon", G4std::ios::in);
  int niso;
  aDataSet >> niso;
  double ** theXsec = new double * [niso+1];
  double ** theEs = new double * [niso+1];
  int * isoSize = new int[niso+1];
  G4String * isoLabels = new G4String[niso+1];
  double dummy;
  for(i=0; i<niso; i++)
  {
    aDataSet >> isoLabels[i];
    aDataSet >> dummy >> dummy;
    aDataSet >> isoSize[i];
    theXsec[i] = new double[isoSize[i]];
    theEs[i] = new double[isoSize[i]];
    for(j=0; j<isoSize[i]; j++)
    {
      aDataSet >> theEs[i][j];
      theEs[i][j]/=1000000.;
      aDataSet >> theXsec[i][j];
    }
  }
  
  // now dump out
  cout << "opt a4"<<G4endl;
  cout << "zone 2 3"<<G4endl;
  cout << "opt pto"<<G4endl;
  cout << "opt nbox"<<G4endl;
  for(i=0; i<counter; i++)
  {
    if(theList[i]=="UNCHANGED") i++;
    int it;
    for(it=0; it<niso; it++) 
    { 
      if(isoLabels[it]==theList[i]) break;
    }
    
    cout << "ve/cr E_";
    cout <<isoLabels[it];
    cout <<"("<<isoSize[it];
    cout <<") r ";
    for(j=0; j<isoSize[it]; j++)
    {
      if(j==30*(j/30)) cout << " _ "<<G4endl;
      cout << theEs[it][j]<<" ";
    }
    cout << G4endl;
    cout << "ve/cr X_"<<isoLabels[it]<<"("<<isoSize[it]<<") r ";
    for(j=0; j<isoSize[it]; j++)
    {
      if(j==30*(j/30)) cout << " _ "<<G4endl;
      cout << theXsec[it][j]<<" ";
    }
    cout << G4endl;
  }
  
  // now make the plots
  int pictureCount = 0;
  cout << "ve/cr err(11) r 11*0.000000001"<<G4endl;
  for(i=0; i<counter; i++)
  {
    if(6*(pictureCount/6) == pictureCount) 
    {
       if(pictureCount!=0) cout <<"close 66"<<G4endl;
       cout << "fo/file 66 silicon28"<<pictureCount/6<<".ps"<<G4endl<<"meta 66 -113"<<G4endl;
    }
    if(theList[i]=="UNCHANGED") i++;
    for(j=0; j<niso; j++)
    {
      if(isoLabels[j]==theList[i]) break;
    }
    double max=0;
    for(k=0; k<isoSize[j]; k++)
    {
      if(max<theXsec[j][k]) max = theXsec[j][k];
    }
    cout << "title_gl 'isotopes produced by neutrons on 28Si'"<<G4endl;
    cout << "null -10 110 0 "<<1.5*max<<G4endl;
    cout << "atit 'E?n! \"m#MeV\"n#' '[s]?"<<isoLabels[j]<<"! \"m#barn\"n#'"<<G4endl;
    cout << "text 0 "<<0.9*1.5*max<<" '"<<isoLabels[j]<<"' 0.3"<<G4endl;
    cout << "ve/cr err_"<<isoLabels[j]<<"("<<isoSize[j]<<") r "<<isoSize[j]<<"*0.0000000001"<<G4endl;
    cout << "hpl/err E_"<<isoLabels[j]<< " X_"<<isoLabels[j]
                 <<" err_"<<isoLabels[j]<<" err_"<<isoLabels[j]<<" "
                 << isoSize[j] << " 19 0.1"<<G4endl;
    cout << "hpl/err ee "<<theList[i]<<" err err"<<theList[i]<<" 11 20 0.3"<<G4endl;
    pictureCount++;
  }
  cout << "close 66"<<G4endl;
  
  // some debugging
  
  double sumX = 0;
  double testEnergy = 1.;
  for(i=0; i<niso; i++)
  {
    double lowE, highE, lowX, highX;
    lowE = 0;
    for(j=0; j<isoSize[i]; j++)
    {
      if(theEs[i][j]>testEnergy)
      {
        highE = theEs[i][j];
        highX = theXsec[i][j];
        break;
      }
      lowE = theEs[i][j];
      lowX = theXsec[i][j];
    }
    double avX = lowX+(testEnergy-lowE)*(highX-lowX)/(highE-lowE);
    sumX += avX;
  }
  G4cerr << "the total cross-section at "<<testEnergy<<" MeV is "<<sumX<<G4endl;
}
