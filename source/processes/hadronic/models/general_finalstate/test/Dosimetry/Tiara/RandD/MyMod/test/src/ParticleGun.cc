#include "ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Interfaces/IVector.h"
#include "Interfaces/IVectorFactory.h"
#include "GunMessenger.hh"
#include "stdlib.h"
#include "AnalyzerLin.hh"

//#include "stdafx.h"
IVectorFactory* alGetVectorFactory();

struct __tagGunStruct{
  double minEnergy;
  double BinWidth;
  double b1;
  double b2;
};

static unsigned long g_ulBins;

static double* g_pBinsSurface=NULL;
static struct __tagGunStruct* g_pBinsInfo=NULL;

/*==================================================================================*/
float ran2(long* idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy=0;
  static long iv[32];
  float temp;

  if(*idum <-0){
    if(-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for(j=39;j>=0;j--){
      k = (*idum)/53668;
      *idum = 40014*(*idum-k*53668)-k*12211;
      if(*idum < 0) *idum += 2147483563;
      if(j<32) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k=(*idum)/53668;
  *idum=40014*(*idum-k*53668)-k*12211;
  if(*idum < 0) *idum += 2147483563;
  k = idum2;
  idum2 = 40692*(idum2-k*52774)-k*3791;
  if(idum2 < 0) idum2 += 2147483399;
  j = iy/(1+2147483562/32);
  iy = iv[j]-idum2;
  iv[j] = *idum;
  if(iy < 1) iy += 2147483562;
  if((temp=1./2147483563*iy) > (1.- 1.2e-7)) return 1.2e-7;
  else return temp;
}
/*================================================================================*/

ParticleGun::ParticleGun()
{
  pGun = new G4ParticleGun();
  m_pMessenger = new GunMessenger(this);
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pParticle = pTable->FindParticle("proton");
  pGun->SetParticleDefinition(pParticle);
  pGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  pGun->SetParticleEnergy(43*MeV);
  pGun->SetParticlePosition(G4ThreeVector(0.,0.,-350.*cm));
  m_bCalculateNeutrons = false;
  m_pNeutronData = NULL;
  m_dArea = m_dRArea = 0;
}
ParticleGun::~ParticleGun()
{
  delete pGun;
  delete m_pMessenger;
  if(g_pBinsSurface) delete[] g_pBinsSurface;
  if(g_pBinsInfo) delete[] g_pBinsInfo;
  g_pBinsSurface = (double*)g_pBinsInfo = NULL;
  g_ulBins = 0;
}
void ParticleGun::GeneratePrimaries(G4Event* pEvent)
{
  if(m_bCalculateNeutrons)
    pGun->SetParticleEnergy(GetEnergy());
  pGun->GeneratePrimaryVertex(pEvent);
}

void ParticleGun::CalculateNeutrons(bool bCalculate,char* filename)
{
  G4int i;
  m_dArea=0;
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  if(m_bCalculateNeutrons){
    if(m_pNeutronData)delete m_pNeutronData;
    m_pNeutronData = NULL;
    if(g_pBinsSurface) delete[] g_pBinsSurface;
    if(g_pBinsInfo) delete[] g_pBinsInfo;
    g_pBinsSurface = NULL;
    g_pBinsInfo = NULL;
    m_dArea = m_dRArea = 0;
  }
  m_bCalculateNeutrons = bCalculate;
  if(bCalculate){
    IPoint* pPoint1;
    IPoint* pPoint2;
    pGun->SetParticleDefinition(pTable->FindParticle("neutron"));
    pGun->SetParticlePosition(G4ThreeVector(0.,0.,-325.*cm));
    m_pNeutronData = alGetVectorFactory()->create();
    m_pNeutronData->fromAscii(filename);
    if(m_pNeutronData->nPoints()==0){
      delete m_pNeutronData;
      m_pNeutronData = NULL;
      m_bCalculateNeutrons = false;
      G4cout<<"Error in reading neutron spectrum. Switching to proton gun"<<G4endl;
      return;
    }
    g_ulBins = m_pNeutronData->nPoints();
    g_pBinsSurface = new double[g_ulBins+1];
    g_pBinsInfo = new (struct __tagGunStruct)[g_ulBins+1];
    g_pBinsSurface[0] = m_pNeutronData->point(0)->coordinate(0)*m_pNeutronData->point(0)->value()*0.5;
    g_pBinsInfo[0].minEnergy = 0.;
    g_pBinsInfo[0].BinWidth = m_pNeutronData->point(1)->coordinate(0);
    g_pBinsInfo[0].b1 = 0;
    g_pBinsInfo[0].b2 = m_pNeutronData->point(0)->value();
    if(g_pBinsInfo[0].minEnergy + g_pBinsInfo[0].BinWidth > 10){
      m_dRArea = (2*g_pBinsInfo[0].b2 - g_pBinsInfo[0].b1)/(2*g_pBinsInfo[0].BinWidth)*(g_pBinsInfo[0].minEnergy + g_pBinsInfo[0].BinWidth-10.);
    }
    for(i=0;i<g_ulBins-1;i++){
      pPoint1 = m_pNeutronData->point(i);
      pPoint2 = m_pNeutronData->point(i+1);
      g_pBinsSurface[i+1] = g_pBinsSurface[i]+(pPoint2->coordinate(0)-pPoint1->coordinate(0))*
	(pPoint1->value()+pPoint2->value())*0.5;
      g_pBinsInfo[i+1].minEnergy = pPoint1->coordinate(0);
      g_pBinsInfo[i+1].BinWidth = pPoint2->coordinate(0)-pPoint1->coordinate(0);
      g_pBinsInfo[i+1].b1 = pPoint1->value();
      g_pBinsInfo[i+1].b2 = pPoint2->value();
      if(g_pBinsInfo[i+1].minEnergy>10){
	m_dRArea += (pPoint2->coordinate(0)-pPoint1->coordinate(0))*(pPoint1->value()+pPoint2->value())*0.5;
      }
      else if(g_pBinsInfo[i+1].minEnergy + g_pBinsInfo[i+1].BinWidth > 10){
	m_dRArea = (2*pPoint2->value()-pPoint1->value())/(2*g_pBinsInfo[i+1].BinWidth)*(pPoint2->coordinate(0)-10);
      }
      //      m_dArea += (pPoint2->coordinate(0) - pPoint1->coordinate(0))*(pPoint1->value() + pPoint2->value())/2.;
    }
    pPoint2 = m_pNeutronData->point(g_ulBins-1);
    g_pBinsSurface[g_ulBins] = g_pBinsSurface[g_ulBins-1]+g_pBinsInfo[g_ulBins-1].BinWidth*pPoint2->value()*0.5;
    g_pBinsInfo[g_ulBins].minEnergy = pPoint2->coordinate(0);
    g_pBinsInfo[g_ulBins].BinWidth = g_pBinsInfo[g_ulBins-1].BinWidth;
    g_pBinsInfo[g_ulBins].b1 = pPoint2->value();
    g_pBinsInfo[g_ulBins].b2 = 0.;
    if(g_pBinsInfo[g_ulBins].minEnergy > 10){
      m_dRArea += g_pBinsInfo[g_ulBins].BinWidth*g_pBinsInfo[g_ulBins].b1*0.5;
    }
    else if(g_pBinsInfo[g_ulBins].minEnergy+g_pBinsInfo[g_ulBins].BinWidth > 10){
      m_dRArea = (2*g_pBinsInfo[g_ulBins].b2-g_pBinsInfo[g_ulBins].b1)/(2*g_pBinsInfo[g_ulBins].BinWidth)*(g_pBinsInfo[g_ulBins].minEnergy+g_pBinsInfo[g_ulBins].BinWidth-10);
    }
    m_dArea = g_pBinsSurface[g_ulBins];
  }
  else{
    pGun->SetParticleDefinition(pTable->FindParticle("proton"));
    pGun->SetParticlePosition(G4ThreeVector(0.,0.,-350.*cm));
  }
}


G4double ParticleGun::GetEnergy()
{
  //  IPoint* pPoint1,*pPoint2;
  G4double dTmp2,dTmp1,a,b1,b2;
  G4double randval;
  G4int i,j,k;
  long dummy;
  if(m_bCalculateNeutrons){
    do{
      randval = (double)rand()/((double)RAND_MAX+1.)*g_pBinsSurface[g_ulBins];
      //  G4double area=0;
      if(randval < g_pBinsSurface[0]){
	if(g_pBinsInfo[0].b2 != 0)
	  dTmp2 = sqrt(2*g_pBinsInfo[0].BinWidth*randval/g_pBinsInfo[0].b2);
	else
	  dTmp2 = g_pBinsInfo[0].BinWidth;
      }
      else if(randval >g_pBinsSurface[g_ulBins-1]){
	if(g_pBinsInfo[g_ulBins].b1 != 0)
	  dTmp2 = g_pBinsInfo[g_ulBins].minEnergy + sqrt(2*g_pBinsInfo[g_ulBins].BinWidth*randval/g_pBinsInfo[g_ulBins].b1);
	else
	  dTmp2 = g_pBinsInfo[g_ulBins].minEnergy;
      }
      else{
	i = 0;j = g_ulBins;
	while(i+1 < j){
	  k = (i+j)>>1;
	  if(randval < g_pBinsSurface[k]) j = k;
	  else if(randval > g_pBinsSurface[k]) i=k;
	  else{
	    dTmp2 = g_pBinsInfo[k].minEnergy;
	    goto Energy_found;
	  }
	}
	randval -= g_pBinsSurface[i];
	b1 = g_pBinsInfo[j].b1;
	b2 = g_pBinsInfo[j].b2;
	a = g_pBinsInfo[j].BinWidth;
	dTmp1 = b1*b1 + (b2-b1)*2*randval/a;
	if(dTmp1<0.)
	  dTmp2 = 0;
	else if(fabs(dTmp1)<1e-9)
	  {
	    dTmp2 = b1*a/(b1-b2);
	    if(dTmp2 < 0) dTmp2 = -dTmp2;
	    if(dTmp2 > a) dTmp2 = a;
	  }
	else{
	  dTmp2 = (b1 - sqrt(dTmp1))*a/(b1-b2);
	  dTmp1 = dTmp2 - 2*sqrt(dTmp1)*a/(b2-b1);
	  if(dTmp2 <0) dTmp2 = 0;
	  else if(dTmp2 > a) dTmp2 = a;
	  if(dTmp1 < 0) dTmp1 = 0;
	  else if(dTmp1 > a) dTmp1 = a;
	  if(dTmp2 == 0) dTmp2 = dTmp1;
	  else if(dTmp1 != 0) dTmp2 = ((dTmp2 < dTmp1) ? dTmp2 : dTmp1);
	}
	dTmp2 += g_pBinsInfo[j].minEnergy;
      }
Energy_found:
      dTmp2 +=0;
    }
    while(dTmp2 < 0.07);
    G4double alpha1 = ran2(&dummy)*0.945e-4;
    G4double beta = ran2(&dummy)*2*pi;
    pGun->SetParticleMomentumDirection(G4ThreeVector(alpha1*sin(beta),alpha1*cos(beta),
						     sqrt((1-alpha1)*(1+alpha1))));
    return dTmp2;
  }
  return pGun->GetParticleEnergy();
}

G4double ParticleGun::GetMax()
{
  G4double dTmp=0,retVal=0;
  IPoint* pPoint;
  if(!m_bCalculateNeutrons) return dTmp;
  if(m_pNeutronData==NULL) return dTmp;
  for(G4int i=0;i<m_pNeutronData->nPoints();i++){
    pPoint = m_pNeutronData->point(i);
    if(pPoint->value() >= dTmp){
      retVal = pPoint->coordinate(0);
      dTmp = pPoint->value();
    }
  }
  return retVal;
}

