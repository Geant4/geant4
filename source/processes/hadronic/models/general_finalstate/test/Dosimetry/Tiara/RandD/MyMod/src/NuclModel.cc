#include "NuclModel.h"
#include "globals.hh"
#include "Randomize.hh"

struct VECTOR{
  double x;
  double y;
  double z;
};
double NormalizeDensity(double A,unsigned nSlices)
{
  double Rad = 1.07*fermi*pow(A,1./3.)+1*fermi;
  unsigned i;
  double currSum=0;
  double currRad;
  for(i=0;i<nSlices;i++){
    currRad = (double)i*Rad/(double)nSlices;
    currSum += 1/(1+exp((currRad-Rad)/0.54));
  }
  return A/currSum;
}

double GetDensity(double Norm,double Rad,double A)
{
  double NuclRad = 1.07*fermi*pow(A,1./3.);
  return Norm/(1+exp((Rad-NuclRad)/0.54));
}

double GetPotential(double Rad,double A,double Z)
{
  double Vcou;
  double NuclRad = 1.07*fermi*pow(A,1./3.)+1*fermi;
  if(Z==0) Vcou = 0;
  else if(Rad < NuclRad){
    Vcou = Z*eplus*eplus/MeV/MeV*0.5/coulomb/coulomb/NuclRad;
    Vcou *= (3-Rad*Rad/NuclRad/NuclRad);
  }
  else{
    Vcou = Z*eplus*eplus/coulomb/coulomb/MeV/MeV/Rad;
  }
  double V0 = -52 +((Z==0) ? 33*(A-2*Z)/A : -33*(A-2*Z)/A);
  return -V0/(1+exp((Rad-NuclRad)/0.65))+Vcou;
}

/*==========================================================================

Sechenija. A sega de ...

NP stoi sa neutron-proton, ne za non-polinomial (hard problem :-)
PP - za proton-proton (i za neutron-neutron) 

============================================================================*/

#define NP_ELASTIC_SIZE 10
#define PP_ELASTIC_SIZE 13
#define NP_ANGLES_SIZE 25
#define PP_ANGLES_SIZE 9

typedef struct __tagElasticStr{
  unsigned nEntries;
  double* pEnergies;
  double* pElastic;
  double* pTotal;
}ELASTICSTRUCT;

typedef ELASTICSTRUCT * LPELASTICSTR;

typedef struct __tagDiffElStruct{
  unsigned nAngles;
  double* pAngles;
  double* pProbabilites;
}DIFFELASTICSTRUCT;
typedef DIFFELASTICSTRUCT * LPDIFFELASTICSTR;

typedef struct __tagDiffInelStruct{
  unsigned nEnergies;
  double* pEnergy;
  double* pProbabilites;
}DIFFINELASTICSTRUCT;
typedef DIFFINELASTICSTRUCT * LPDIFFINELASTICSTR;

typedef struct __tagEnergyXSec{
  DIFFELASTICSTRUCT pDiffElasticProbabilities;
  DIFFINELASTICSTRUCT pDiffInelasticProbabilities;
} ENERGYINDEXEDSTRUCT;

typedef ENERGYINDEXEDSTRUCT * LPENERGYSTR;


typedef struct __tagProbs{
  unsigned type;
  unsigned nEntries;
  double* dEnergies;
  LPELASTICSTR pXSec;
  LPENERGYSTR pData;
}PROBSSTR;
typedef PROBSSTR * LPPROBABILITIES;

static double NPEnergies[NP_ELASTIC_SIZE] = {7,19.9,29.7,44.4,66.3,97.4,144.4,190.2,279.1,364.4};
static double PPEnergies[PP_ELASTIC_SIZE] = {10,19.9,29.8,39.6,49.4,83.2,121.1,149.1,167.5,217.4,261.7,305.2,364.6};

static double NPElasticData[NP_ELASTIC_SIZE] = {361.6,214.5,130,78,53.5,41.6,35.9,33.5,25.9,18.9};
static double NPTotalData[NP_ELASTIC_SIZE] = {361.6,214.5,130,78,53.5,41.6,35.9,34.2,34.3,34.9};
static double PPElasticData[PP_ELASTIC_SIZE] = {156.8,80.6,48,36.6,31.6,25.9,24,23.1,23.1,23.5,24.9,25,22.3};
static double PPTotalData[PP_ELASTIC_SIZE] = {156.8,80.6,48,36.6,31.6,25.9,24,23.1,24,28.3,33.6,41.5,47};

static ELASTICSTRUCT ElDataPP = {PP_ELASTIC_SIZE,PPEnergies,PPElasticData,PPTotalData};
static ELASTICSTRUCT ElDataNP = {NP_ELASTIC_SIZE,NPEnergies,NPElasticData,NPTotalData};

static double NPAngles[NP_ANGLES_SIZE] ={-1,-0.99,-0.96,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.96,0.99,1};
static double NPAnglesData[NP_ELASTIC_SIZE][NP_ANGLES_SIZE]=
{
  {0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08},
  {.095,.094,.09,.086,.082,.079,.077,.075,.073,.072,.071,.071,.071,.071,.073,.075,.079,.082,.085,.088,.091,.095,.097,.099,.099},
  {.113,.112,.107,.101,.090,.084,.079,.075,.072,.069,.067,.065,.064,.064,.065,.065,.068,.072,.079,.088,.098,.111,.120,.125,.127},
  {.187,.159,.141,.123,.105,.092,.082,.075,.068,.063,.06,.057,.055,.055,.057,.060,.063,.067,.073,.082,.096,.113,.132,.154,.162},
  {.224,.216,.193,.160,.120,.096,.083,.074,.067,.060,.055,.051,.0486,.0478,.0476,.0480,.0490,.053,.062,.076,.096,.121,.139,.163,.183},
  {.279,.254,.196,.151,.121,.101,.090,.080,.071,.064,.057,.051,.0458,.0436,.0445,.0478,.052,.057,.062,.069,.080,.104,.137,.213,.260},
  {.331,.234,.182,.143,.108,.086,.072,.063,.057,.052,.0494,.048,.0488,.050,.054,.059,.066,.073,.081,.089,.099,.111,.134,.160,.179,},
  {.389,.252,.180,.134,.102,.082,.068,.058,.052,.0482,.0467,.0472,.0493,.052,.056,.062,.068,.074,.081,.091,.104,.123,.148,.180,.226},
  {.329,.235,.170,.126,.089,.067,.052,.043,.0365,.0334,.0331,.0348,.0378,.0422,.0481,.055,.063,.074,.089,.109,.137,.179,.235,.343,.426},
  {.268,.183,.138,.109,.078,.056,.0407,.0313,.0291,.0284,.0293,.0317,.0353,.0398,.0447,.052,.061,.072,.091,.116,.152,.222,.317,.465,.626}
};
static double PPAngles[PP_ANGLES_SIZE]={0,0.2,0.4,0.6,0.7,0.8,0.9,0.96,1};
static double PPAnglesData[PP_ELASTIC_SIZE][PP_ANGLES_SIZE]=
{
  {0.0795,0.0795,0.0795,0.0795,0.0795,0.0795,0.0795,0.0795,0.0795},
  {0.715,0.715,.073,.08,.08,.088,.0935,.0975,.1005},
  {.077,.077,.077,.0785,.08,.081,.0855,.092,.097},
  {.0725,.073,.074,.078,.082,.087,.0935,.101,.108},
  {.069,.069,.072,.0775,.082,.0895,.101,.1125,.123},
  {.07,.07,.07,.073,.078,.087,.150,.1275,.1495},
  {.074,.074,.075,.075,.079,.08,.093,.1155,.139},
  {.0795,.079,.078,.078,.0775,.0775,.0795,.0935,.1105},
  {.081,.081,.0785,.078,.078,.078,.0795,.084,.095},
  {.0735,.0755,.0785,.079,.0805,.0815,.0855,.096,.099},
  {.0615,.064,.0695,.0765,.0825,.091,.150,.1305,.184},
  {.041,.0455,.0575,.078,.092,.1095,.134,.159,.36},
  {.01995,.02375,.028,.0715,.0965,.1325,.1815,.222,.46}
};

static ENERGYINDEXEDSTRUCT PPangles[PP_ELASTIC_SIZE] =
{
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[0]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[1]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[2]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[3]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[4]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[5]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[6]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[7]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[8]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[9]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[10]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[11]},{0}},
  {{PP_ANGLES_SIZE,PPAngles,PPAnglesData[12]},{0}}
};

static ENERGYINDEXEDSTRUCT NPangles[NP_ELASTIC_SIZE]=
{
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[0]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[1]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[2]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[3]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[4]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[5]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[6]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[7]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[8]},{0}},
  {{NP_ANGLES_SIZE,NPAngles,NPAnglesData[9]},{0}}
};

static PROBSSTR ElData[2]=
{
  {EQUAL_TYPE,PP_ELASTIC_SIZE,PPEnergies,&ElDataPP,PPangles},
  {DIFF_TYPE,NP_ELASTIC_SIZE,NPEnergies,&ElDataNP,NPangles}
};

/*-----------------------------------------------------------------

dannite dotuk. sega coda

------------------------------------------------------------------*/

double GetElasticData(double Ecm,unsigned nType,double& Total)
{
  Total = 0;
  if(nType>DIFF_TYPE) return 0;
  if(Ecm > 800){
    Total = 35.5;
    return 0;
  }
  signed i=0,j=ElData[nType].nEntries-1,k;
  bool bExact=false;
  if(ElData[nType].dEnergies[0] > Ecm) i=-1;
  else{
    do{
      k = (j+i)/2;
      if(ElData[nType].dEnergies[k] > Ecm){
	j = k;
      }
      else if(ElData[nType].dEnergies[k] < Ecm){
	i = k;
      }
      else{
	bExact = true;
	break;
      }
    }
    while(i+1<j);
  }
  if(bExact){
    Total = ElData[nType].pXSec->pTotal[k];
    return ElData[nType].pXSec->pElastic[k];
  }
  //tuk se vzimat tezi ot i i ot j;
  if(i==-1){
    Total = ElData[nType].pXSec->pTotal[0];
    return ElData[nType].pXSec->pElastic[0];
  }
  else if(j==(signed)ElData[nType].nEntries){
    Total = ElData[nType].pXSec->pTotal[j-1];
    return ElData[nType].pXSec->pElastic[j-1];
  }
  Total = (ElData[nType].pXSec->pTotal[i]+ElData[nType].pXSec->pTotal[j])/2.;
  return ((ElData[nType].pXSec->pElastic[i]+ElData[nType].pXSec->pElastic[j])/2.);
}

static double* CompouseWorkingArray(double Ecm,unsigned nType,double* pArr)
{
  if(pArr==NULL) pArr = new double[ElData[nType].nEntries];           /* tova ne biva da se izpylnjava i oshe
									 masiva treve da e dostatychno goljam*/
  unsigned i=0,j=ElData[nType].nEntries,k;
  double* dHelp,*dHelp1/*,dAngles*/;
  bool bExact=false;
  if(Ecm > 800){
    j = ElData[nType].pData[0].pDiffElasticProbabilities.nAngles;
    dHelp =  ElData[nType].pData[0].pDiffElasticProbabilities.pProbabilites;
    for(i=0;i<j;i++) pArr[i] = dHelp[i];
    return pArr;
  }
  else if(Ecm < ElData[nType].dEnergies[0]){
    j = ElData[nType].pData[ElData[nType].nEntries].pDiffElasticProbabilities.nAngles;
    dHelp = ElData[nType].pData[ElData[nType].nEntries].pDiffElasticProbabilities.pProbabilites;
    for(i=0;i<j;i++) pArr[i] = dHelp[i];
    return pArr;
  }
  else{
    do{
      k = (i+j)/2;
      if(ElData[nType].dEnergies[k] > Ecm) j = k;
      else if(ElData[nType].dEnergies[k]<Ecm) i = k;
      else{
	bExact = true;
	break;
      }
    }
    while(i+1<j);
    if(bExact){
      j = ElData[nType].pData[k].pDiffElasticProbabilities.nAngles;
      dHelp = ElData[nType].pData[k].pDiffElasticProbabilities.pProbabilites;
      for(i=0;i<j;i++) pArr[i] = dHelp[i];
      return pArr;
    }
    else{
      //Predpolagame che i dvata seta sa s ednakva golemina
      j = ElData[nType].pData[k].pDiffElasticProbabilities.nAngles;
      dHelp = ElData[nType].pData[i].pDiffElasticProbabilities.pProbabilites;
      dHelp1 = ElData[nType].pData[j].pDiffElasticProbabilities.pProbabilites;
      for(i=0;i<j;i++) pArr[i] = (dHelp[i]+dHelp1[i])/2.;
    }
  }
  return pArr;
}

double SampleAngle(double Ecm,unsigned nType)
{
  //Zasega imame max 25 elementa, zatova pArr e s razmer 30 > 25
  double pArr[30],*pAng;
  double Int,newInt;
  unsigned nSize,i;
  nSize = ElData[nType].pData[0].pDiffElasticProbabilities.nAngles;
  pAng = ElData[nType].pData[0].pDiffElasticProbabilities.pAngles;
  CompouseWorkingArray(Ecm,nType,pArr);
  Int= -pArr[0]*pAng[0]/2.;
  for(i=0;i<nSize-1;i++){
    Int += (pArr[i]+pArr[i+1])*(pAng[i]-pAng[i-1])/2.;
  }
  //v Int - liceto na tova govedo
  double rand = G4UniformRand()*Int;
  /*double ret = pAng[0];*/
  newInt = pArr[0]*pAng[0]/2.;
  if(rand < newInt){
    return sqrt(-2*pAng[0]*rand/pArr[0]);
  }
  double Det,t1,t2;
  for(i=0;i<nSize;i++){
    newInt = -(pArr[i]+pArr[i+1])*(pAng[i+1]-pAng[i])/2.;
    if(rand < newInt){
      Det = pArr[i]*pArr[i]-2*(pArr[i+1]-pArr[i])*(newInt-rand)/(pAng[i+1]-pAng[i]);
      if(Det<0){
	return pAng[i];
      }
      else if(Det==0){
	return -pArr[i]*(pAng[i+1]-pAng[i])/(pArr[i+1]-pArr[i]);
      }
      t1 = -(-pArr[i]+sqrt(Det))*(pAng[i+1]-pAng[i])/(pArr[i+1]-pArr[i]);
      t2 = -(-pArr[i]-sqrt(Det))*(pAng[i+1]-pAng[i])/(pArr[i+1]-pArr[i]);
      if(t2 < t1) t1 = t2;
      return t1;
    }
  }
  return pAng[nSize-1];
}

double GetAngle(double E1,double M1,double E2,double M2)
{
  double v1 = sqrt(2*M1*E1);
  double v2 = sqrt(2*M2*E2);
  v2 = (M1*v1+M2*v2)*(M1*v1+M2*v2)/2/(M1+M2);
  if(M1==M2) return SampleAngle(v2,EQUAL_TYPE);
  double angle = SampleAngle(v2,DIFF_TYPE);
  if(angle<-1) angle = -1;
  else if(angle > 1) angle=1;
  if(angle+M2/(M1+M2)==0){
    return 1;
  }
  return cos(atan(sin(acos(angle))/(angle+M2/(M1+M2))));
}

void ModifyMomentum(double costheta,double phi,VECTOR& Mom)
{
  double quat[4];
  double quat1[4];
  double res[4];
  double MomSize = sqrt(Mom.x*Mom.x+Mom.y*Mom.y+Mom.z*Mom.z);
  res[0] = sin(asin(sqrt(1-costheta*costheta))/2.);
  res[1] = sqrt(1-res[0]*res[0]);
  quat[0] = Mom.y/MomSize*res[0]*((costheta<0) ? 1 : -1);
  quat[1] = Mom.x/MomSize*res[0]*((costheta<0) ? 1: -1);
  quat[2] = 0;
  quat[3] = res[1];
  quat1[0] = quat1[1] = 0;
  quat1[2] = sin(phi/2);
  quat1[3] = cos(phi/2);
  //sega gi umnojavame quat*quat1 - purvo po z na phi sled tova obratnata na tazi kojato kara p po ostta z
  res[0] = quat[0]*quat1[0] - quat[1]*quat1[1]-quat[2]*quat1[2]-quat[3]*quat1[3];
  res[1] = quat[0]*quat1[1] + quat[1]*quat1[0] + quat[2]*quat1[3]-quat[3]*quat1[2];
  res[2] = quat[0]*quat1[2] - quat[1]*quat1[3] + quat[2]*quat1[0]+ quat[3]*quat1[1];
  res[3] = quat[0]*quat1[3] + quat[1]*quat1[2] - quat[2]*quat1[1] + quat[3]*quat1[0];
  // a ot dolnite 5 reda ne trjabva da ima nujda - tova e umnojenie na edinici => trjabva da e edinica
  MomSize = sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]+res[3]*res[3]);
  res[0] /= MomSize;
  res[1] /= MomSize;
  res[2] /= MomSize;
  res[3] /= MomSize;
  //sega pravim matricata i vyrtim :-). res e normaliziran
  quat[0] = res[0]*Mom.x - res[1]*Mom.y-res[2]*Mom.z-res[3]*0;
  quat[1] = res[0]*Mom.y + res[1]*Mom.x + res[2]*0-res[3]*Mom.z;
  quat[2] = res[0]*Mom.z - res[1]*0 /*tova compiler-a she go mahne*/+ res[2]*Mom.x+ res[3]*Mom.y;
  quat[3] = res[0]*0/*tova compiler-a she go mahne*/ + res[1]*Mom.z - res[2]*Mom.y + res[3]*Mom.x;
  res[0] *=-1;
  res[1] *=-1;
  res[2] *=-1;
  Mom.x = quat[0]*res[0] - quat[1]*res[1] - quat[2]*res[2] - quat[3]*res[3];
  Mom.y = quat[0]*res[1] + quat[1]*res[0] + quat[2]*res[3] - quat[3]*res[2];
  Mom.z = quat[0]*res[2] - quat[1]*res[3] + quat[2]*res[0] + quat[3]*res[1];
  //i poslednija red treve da dava nula
}

double GetXSection(double E,unsigned nType)
{
  double XSec=0;
  if(nType==EQUAL_TYPE){
    if(E==0) return pi*1.5*1.5;
    if(E<40){
      XSec = (-1174.8/E/E+3088.5/E+5.3107);
    }
    else if(E<310){
      XSec = (93074/E/E - 11.148/E + 22.429);
    }
    else{
      XSec = (887.37/E + 0.053331*E + 3.5475);
    }
  }
  if(E<40)
    XSec = -5057.4/E/E + 9069.2/E + 6.9466;
  else if(E<400)
    XSec = 239380/E/E + 1802/E + 27.147;
  else XSec = 34.5;
  return ((XSec<0.) ? 0. : XSec);
}
