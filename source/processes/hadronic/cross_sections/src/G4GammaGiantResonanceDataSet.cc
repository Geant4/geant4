// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GammaGiantResonanceDataSet.cc,v 1.2 2000-09-27 12:23:27 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: GammaGiantResonanceDataSet for gamma+A cross sections
// M.V. Kossov, ITEP(Moscow), 10-SEP-00
// 

#include "G4GammaGiantResonanceDataSet.hh"

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4GammaGiantResonanceDataSet::GetCrossSection(const G4DynamicParticle* aPart,
                                                       const G4Element* anEle)
{
  const G4double kineticEnergy = aPart->GetKineticEnergy()/GeV;
  const G4double targetAtomicNumber = anEle->GetN();
  const G4double nTargetProtons = anEle->GetZ();
  const G4double nTargetNeutrons = targetAtomicNumber-nTargetProtons;
  G4double AA=targetAtomicNumber+targetAtomicNumber;
  G4int i=6;
  if     (AA<ANucl(0)+ANucl(1)) i=0;
  else if(AA<ANucl(1)+ANucl(2)) i=1;
  else if(AA<ANucl(2)+ANucl(3)) i=2;
  else if(AA<ANucl(3)+ANucl(4)) i=3;
  else if(AA<ANucl(4)+ANucl(5)) i=4;
  else if(AA<ANucl(5)+ANucl(6)) i=5;
  G4double f = AbsorbtionByNucleus(i,kineticEnergy);
  G4double sqCharge = 2*nTargetNeutrons + 3*nTargetProtons;
  return (f/SumQQ(i))*sqCharge*millibarn;        // Scaling
}

// Combined member function for different nuclei
G4double G4GammaGiantResonanceDataSet::AbsorbtionByNucleus(G4int i, G4double E)
{
  static const G4int nOfN  = 7;
  //                            -H-   -D-   -Be-   -C-   -Al-  -Cu-  -Pb-
  static G4double X0[nOfN] = {  0.1, 0.001,   0. ,  5.0,  5.0,  0.0,  0.0};
  static G4double XS[nOfN] = {100. , 1.   , 300. ,  20.,  20.,  10.,  0.0};
  static G4double XM[nOfN] = {175.8, 130.5, 100.5, 71.9, 70.6,  10.,  0.0};
  static G4double YS[nOfN] = {  1. , 1.   ,   4. ,   6.,  12., 100.,  0.0};
  static G4double YM[nOfN] = {106.3, 86.5 ,  44.4,  50., 51.6, 42.2,  0.0};
  // Cross section gamma-proton taken from PDG figure
  static const G4int nPtH  = 39;
  static G4double XH[nPtH] = {17.1,19.5,21.9,23.0,24.4,27.5,29.1,30.5,32.8,34.9,
                              36.5,38.3,40.6,42.9,45.1,47.0,49.0,51.0,52.6,54.3,
                              55.6,56.9,58.2,59.3,60.4,62.0,65.0,69.1,71.4,73.8,
                              75.9,79.8,86.0,93.6,99.3,103.9,111.1,116.9,134.5};
  static G4double YH[nPtH] = {11.1,16.6,21.9,30.5,43.5,54.0,56.5,53.1,44.0,33.0,
                              23.8,21.0,18.4,19.3,22.7,24.9,28.0,29.3,26.0,22.5,
                              21.4,21.6,22.5,23.0,20.8,17.5,15.7,16.0,15.3,15.0,
                              15.0,14.1,13.9,13.2,13.0,12.9,12.9,12.7,11.9};
  static G4double YN[nPtH] = {20.5,30.0,38.9,54.3,80.0,93.3,96.1,95.2,82.7,68.5,
                              56.6,48.6,41.5,40.8,44.6,48.3,51.6,50.1,46.9,42.5,
                              39.3,37.0,37.9,36.6,35.0,32.0,30.1,30.4,30.3,29.8,
                              28.5,27.7,27.3,25.5,25.0,24.8,24.5,24.0,22.5};

  // Cross section gamma-D for low E_gamma are taken from the figure in the BLACK BOOK
  static const G4int nPtD  = 20;
  static G4double XD[nPtD] = {16.5,17.5,18.5, 27.2, 32.6,38.7,50.6,53.9,57.0,60.2,
                              63.6,69.0,75.0, 79.1, 85.4,88.9,91.1,94.0,96.5,99.5};
  static G4double YD[nPtD] = {89.6,92.2,94.4,104.0,101.8,98.1,84.0,78.5,75.4,70.9,
                              65.0,58.4,52.0, 45.4, 37.0,34.8,34.0,35.5,37.0,38.0};
  // Cross section gamma-Be (J.Ahrens Nucl.Phys. A335(1980)67-74)"Pion threshold"
  static const G4int nPtPT = 25;
  static G4double XPT[nPtPT]={32.3,35.7,38.6,42.1,45.6,49.3,52.2,54.9,57.0,60.0,
                              63.0,66.0,68.3,70.1,73.7,77.4,80.5,83.1,86.3,89.9,
                              92.6,94.4,96.8,99.3,102.2};
  static G4double YPT[nPtPT]={ 5.2, 4.5, 4.7, 5.2, 6.0, 7.3, 8.7, 9.9,11.0,13.9,
                              15.9,17.5,19.5,21.0,23.8,26.9,29.5,31.7,33.9,36.0,
                              37.2,38.0,39.0,40.0,39.8};
  // Cross section gamma-Be (J.Ahrens Nucl.Phys. A335(1980)67-74)"Low energy part"
  static const G4int nPtB  = 19;
  static G4double XB[nPtB] = { 4.0, 5.1, 5.5, 6.6, 7.6, 8.8, 9.7,10.5,11.3,12.2,
							   13.1,14.0,15.5,17.2,20.1,23.0,26.5,29.4,32.3};
  static G4double YB[nPtB] = { 8.9, 9.8,24.7,56.0,59.5,58.6,52.9,49.0,44.2,40.1,
							   36.0,34.2,25.2,20.9,15.6,12.2, 9.5, 7.0, 5.2};
  // gamma-C
  static const G4int nPtC  = 23;
  static G4double XC[nPtC] = {46.7, 51.2, 53.0,54.5,55.6,56.7,57.8,58.7,59.7,61.9,
                              62.9, 65.5, 67.5,69.1,71.9,75.0,84.0,87.2,89.1,93.0,
                              98.0,107.0,113.5};
  static G4double YC[nPtC] = { 0.3,  2.1,  5.2, 9.6,15.1,23.9,36.3,47.8,56.5,53.5,
                              60.3, 54.0, 40.0,28.5,34.7,25.3,18.5,15.5,21.0,13.5,
							  10.5, 14.5, 10.0};
  // gamma-AL
  static const G4int nPtA  = 25;
  static G4double XA[nPtA] = {28.2,31.1,35.3, 37.4, 39.7,43.8,45.0,47.0,50.8,52.0,
                              55.5,57.1,59.0, 63.0, 66.4,69.5,72.3,76.2,78.6,82.8,
                              85.5,92.1,97.8,101.0,108.8};
  static G4double YA[nPtA] = { 1.0, 4.1, 7.9, 15.2, 21.6,27.1,36.9,42.7,44.0,52.5,
                              53.2,62.0,58.5, 54.5, 48.5,42.6,38.0,34.7,37.7,35.0,
							  33.6,29.0,22.5, 22.8, 18.5};
  // gamma-cU
  static const G4int nPtU  = 55;
  static G4double XU[nPtU] = { 9.84375,10.15630,10.46880,10.78130,11.09380,
                              11.40630,11.71880,12.03130,12.34380,12.65630,
                              12.96880,13.28130,13.59380,13.90630,14.21880,
                              14.53130,14.84380,15.15630,15.46880,15.78130,
                              16.09380,16.40630,16.71880,17.03130,17.34380,
                              17.65630,17.96880,18.28130,18.59380,18.90630,
                              19.21880,19.53130,19.84380,20.15630,20.46880,
                              20.78130,21.09380,21.40630,21.71880,22.03130,
                              22.34380,22.65630,22.96880,23.28130,23.59380,
                              23.90630,24.21880,24.53130,25.15630,25.78130,
							  26.40630,27.03130,27.65630,30.00001,33.00001};
  static G4double YU[nPtU] = { 0.6, 1.3, 2.5, 3.5, 4.6, 5.3, 6.2, 7.2, 8.3,10.1,
                              10.1,12.9,13.8,15.7,18.6,22.4,23.2,24.2,26.3,27.5,
                              30.0,29.0,30.5,28.3,27.6,25.5,25.0,25.8,26.4,24.4,
                              27.3,25.0,25.0,24.9,21.7,21.2,23.8,23.2,21.9,21.9,
                              20.2,18.7,20.8,19.2,17.5,15.2,14.5,15.5,11.9,12.0,
                               9.8, 8.8, 8.8, 5.4, 2.5};
  // gamma-Pb is made as a function
  G4double f=0.;
  //if(i>1&&E>0.3) f= ANucl(i)*HighEnergyOld(E);  // General High Energy approximation
  //else if(i>1&&E>0.106) f=LinearFit(E*335., nPtPT, XPT, YPT)*SumQQ(i)/11.1/SumQQ(2);
  if(i>1&&E>0.086) f= ANucl(i)*HighEnergyNew(E);  // General High Energy approximation
  else if(i>1&&E>PbDataEnergy(i)) f=(23.73333-106.6667*E)*CorA(E,i)*SumQQ(i)/SumQQ(6);
  else if(!i)                                // Hidgogen
  {
    if(E>ThresholdEnergy(i))
      return f=(YS[i]/YM[i])*LinearFit(log10(E/X0[i])*XM[i]/log10(XS[i]/X0[i]),nPtH,XH,YH);
    else return 0.;
  }
  else if(i==1)                                // Deuterium
  {
    if(E>ThresholdEnergy(0))
      return f=(YS[0]/YM[0])*LinearFit(log10(E/X0[0])*XM[0]/log10(XS[0]/X0[0]),nPtH,XH,YN);
    else if(E>ThresholdEnergy(i))
      return f=0.01*pow(10,(LinearFit(log10(E/X0[i])*XM[i]/log10(XS[i]/X0[i]),nPtH,XD,YD)
                         *log10(YS[i]/0.01)/YM[i]));
    else return 0.;
  }
  else if(i==2)                                // Berilium
  {
    if(E>ThresholdEnergy(i))
	  return f=(YS[i]/YM[i])*LinearFit((E*1000.-X0[i])*XM[i]/XS[i],nPtB,XB,YB)*CorA(E,i);
    else return 0.;
  }
  else if(i==3)                                // Carbon
  {
    if(E>ThresholdEnergy(i))
	  return f=(YS[i]/YM[i])*LinearFit((E*1000.-X0[i])*XM[i]/XS[i],nPtC,XC,YC)*CorA(E,i);
    else return 0.;
  }
  else if(i==4)                                // Aluminum
  {
    if(E>ThresholdEnergy(i))
	  return f=(YS[i]/YM[i])*LinearFit((E*1000.-X0[i])*XM[i]/XS[i],nPtA,XA,YA)*CorA(E,i);
    else return 0.;
  }
  else if(i==5)                                // Copper
  {
    if(E>ThresholdEnergy(i))
	  return f=(YS[i]/YM[i])*LinearFit((E*1000.-X0[i])*XM[i]/XS[i],nPtU,XU,YU)*CorA(E,i);
    else return 0.;
  }
  else if(i==6)                                // Lead
  {
    if(E>ThresholdEnergy(i))
	  return f=(640./(1.+pow((E*E-0.00018)/E/0.00405,2)))*CorA(E,i);
    else return 0.;
  }
  return f;
}

// Gamma-A in the range 0.15-5 GeV (from PR-89-002) 
G4double G4GammaGiantResonanceDataSet::HighEnergyOld(G4double E)
{
  static G4double A=-79.55;
  static G4double B=17.673;
  static G4double S=0.0125;
  static const G4int nPtA = 22; // Coveres 0.1475 - 5.299 GeV range
  static G4double X[nPtA] = { 8.7,10.7,12.7,14.7,16.1,18.0,19.7,21.5,23.0,24.4,25.6,
                             27.6,29.5,30.6,32.9,38.1,43.5,50.0,55.5,61.0,65.6,72.};
  static G4double Y[nPtA] = { 2.0, 4.9, 8.9,12.9,17.6,23.6,28.0,31.3,32.8,31.0,28.3,
                             24.5,21.5,17.9,16.0,14.5,12.0,10.5, 9.5, 8.0, 7.5, 7.};
  if(E<5.)
    return S*LinearFit(A+B*log(E*1000.), nPtA, X, Y);
  else
  {
    G4double s=0.880355+1.876545*E;
    return 0.071*pow(s,0.075)+0.12*pow(s,-0.46);
  }
}

// Gamma-A for E>0.086 GeV (from PL-127B-331-CHOLLRT-83 & PR-54-1688-BIANCHI-96) 
G4double G4GammaGiantResonanceDataSet::HighEnergyNew(G4double E)
{
  static const G4int nPtA = 73; // Coveres 0.086 - 1.163 GeV range
  static G4double X[nPtA] = {0.086,0.092,0.098,0.103,0.106,0.112,0.119,0.127,0.132,.1424,
                             .1475,.1525,.1572,.1619,.1669,.1719,.1770,.1825,.1879,.1926,
                             .1978,.2025,.2076,.2128,.2177,.2229,.2282,.2328,.2373,.2423,
                             .2473,.2522,.2572,.2624,.2677,.2729,.2783,.2825,.2867,.2926,
                             .2973,0.301,0.317,0.334,0.369,0.389,0.408,0.439,0.465,0.490,
                             0.514,0.540,0.568,0.598,0.616,0.636,0.664,0.684,0.717,0.751,
                             0.768,0.788,0.817,0.840,0.865,0.895,0.908,0.936,0.973,1.044,
                             1.081,1.119,1.163};
  static G4double Y[nPtA] = {.0579,.0603,.0603,.0579,.0555,.0579,.0579,.0579,.0724,.0688,
                             .0644,.0625,.0810,.0865,.0965,.0999,.1265,.1387,.1498,.1631,
                             .1653,.1698,.1942,.2008,.2230,.2330,.2508,.2663,.2796,.2900,
                             .3074,.3218,.3307,.3440,.3440,.3584,.3651,.3939,.3839,.3961,
                             .3950,.4015,.4195,.4212,.4071,.3836,.3547,.3178,.2920,.2516,
                             .2509,.2357,.2211,.2188,.2026,.2042,.2064,.1924,.1884,.1860,
                             .1898,.1763,.1803,.1837,.1811,.1813,.1691,.1603,.1682,.1572,
                             .1509,.1516,.1407};
  if(E<1.16)
    return LinearFit(E, nPtA, X, Y);
  else
  {
    G4double s=0.880355+1.876545*E;
    return 0.071*pow(s,0.075)+0.12*pow(s,-0.46);
  }
}

// Correction function for Be,C @@ Move to header
G4double G4GammaGiantResonanceDataSet::LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN)
{
  G4double Xj=XN[0];
  G4double Xh=XN[N-1];
  if(X<=Xj) return Xj; //-----+
  else if(X>=Xh) return Xh;//-|
  G4double Xp; //             |
  G4int j=0;   //             |
  while (X>Xj && j<N)//<------+
  {
    j++;
    Xp=Xj;
    Xj=XN[j];
  }
  return YN[j]-(Xj-X)*(YN[j]-YN[j-1])/(Xj-Xp);
}


