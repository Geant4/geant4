#include "AnalyzerLin.hh"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IVector.h"
#include "Interfaces/IHistoFactory.h"
#include "Interfaces/IVectorFactory.h"
#include "Interfaces/IPlotter.h"
#include "Interfaces/IAxis.h"

#include "cfortran.h"
#include "hplot.h"
#include "hbook.h"
#include "higz.h"
#include "kuip.h"


#include "strings.h"
#include <ctype.h>

//#include "stdafx.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        Types and local variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define GRID_BIT 1
#define DVY_BIT 1<<1
#define DVX_BIT 1<<2
#define STA_BIT 1<<3
#define TIC_BIT 1<<4
#define BOX_BIT 1<<5
#define LINX_BIT 1<<6
#define LINY_BIT 1<<7
#define AST_BIT 1<<8

#define DEFAULT_OPTIONS GRID_BIT|DVY_BIT|DVX_BIT|BOX_BIT|LINX_BIT|LINY_BIT

enum STATLINETYPES {STAT_INITIALIZE,STAT_SOLIDLINE,STAT_DASHLINE,STAT_DOTLINE,STAT_DASHDOTLINE,STAT_DOTDOTLINE};

unsigned int OptArray[9] ={GRID_BIT,DVY_BIT,DVX_BIT,STA_BIT,TIC_BIT,BOX_BIT,LINX_BIT,LINY_BIT,AST_BIT};

struct __tagHistoSet{
  char* String;
  float value;
  __tagHistoSet* pNext;
};

typedef __tagHistoSet HISTOSET;
typedef HISTOSET* LPHISTOSET;

struct __tagHistoOptions{
  struct{
    int AxisColor;
    int TitleColor;
    int AxisTitleColor;
    int HistoColor;
    int FuncColor;
  }Colors;
  unsigned cFlags;
  char* XAxis;
  char *YAxis;
  float minX;
  float minY;
  float maxX;
  float maxY;
  char HistoStyle[2];
  LPHISTOSET pSetOptions;
};

typedef __tagHistoOptions HISTOPT;
typedef HISTOPT* LPHISTOPT;

struct _tagHistoNode{
  G4int          ID;
  IHistogram1D*  pHisto;
  IHistogram1D** pAdditionalHistos;
  unsigned char  cAdditionalHistos;
  unsigned*       idAdditionalHistos;
  IVector*       pDataSet;
  IVector*       pBinErrors;
  IVector*       pBinEntries;
  unsigned long  pViewer;
  unsigned long  cEvents;
  LPHISTOPT      pOptions;
  float          fScale;
  _tagHistoNode* pNext;
  _tagHistoNode* pPrev;
};
typedef _tagHistoNode HISTONODE;
typedef HISTONODE* LPHISTONODE;

static LPHISTONODE HistoList = NULL;

static IHistoFactory* pHistoFactory=NULL;
static char*       szStoreName=NULL;
static IVectorFactory* pVectorFactory=NULL;

static char* g_szGlobalTitle = NULL;

static G4bool bViewerIsShown;
static G4bool bOptimizeHistos=true;

static unsigned int nZones=0;

static char OptionsArray[9][5] = {"NGRI", "DVYR", "DVXR",
				  "NSTA", "TIC ", "BOX ",
				  "LINX", "LINY", "AST "};
static char NegOptionsArray[9][5] = {"GRID","DVYI","DVXI",
				     "NSTA","NTIC","NBOX",
				     "LOGX","LOGY","NAST"};

static unsigned short int AAType=AA_NONE;

static unsigned char cHistos=0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static  G4int FindID();
static LPHISTONODE FindHisto(G4int ID,LPHISTONODE* ppPrev);
static void ScaleHisto(LPHISTONODE pCurr,G4double val);
static double alFindPointFor(G4double xCoord,IVector* pVector,G4double dMaxEnergy);
static void alShiftZones(unsigned StartFrom);
static unsigned int FindOption(G4String szKey);
static unsigned int FindColor(G4String szKey);
static unsigned int FindAxis(G4String szKey);
static LPHISTOSET FindSetNode(LPHISTOPT pOpt,G4String szKey,LPHISTOSET* pPrev);
static void CalculateZones(unsigned int *xDim,unsigned int *yDim);
static void SetOptions(LPHISTOPT pOpt);
static void AddStat(unsigned char bStatType,char* szString);
static void Plot(LPHISTONODE pCurr);
static void PlotHisto(unsigned int ID,bool bRecalculateZones);
static void alAddData(IVector* pData);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local functions implementation
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

static double TentFilter(double f1,double f2,double f3,double x1,double x2,double xCoord)
{
  double binWidth = x2-x1;
  double res=0.f;
  if(xCoord-binWidth>x1){
    if(xCoord+binWidth<x2) return f2;
    res = 1.f-(x2-xCoord-binWidth)/binWidth*0.5;
    f2 = f2*res*(x2-xCoord+binWidth);
    res = res - 1.f;
    f3 = f3*res*(binWidth-x2+xCoord);
    res = (f2+f3)/binWidth*0.5;
  }
  else if(xCoord+binWidth<x2){
    res = 1.f-(x1-xCoord-binWidth)/binWidth*0.5;
    f1 = f1*res*(x1-xCoord+binWidth);
    res = res - 1.f;
    f2 = f2*res*(binWidth-x1+xCoord);
    res = (f1+f2)/binWidth*0.5;
  }
  else{
    res = 1.f-(x1-xCoord-binWidth)/binWidth*0.5;
    f1 = f1*res*(x1-xCoord+binWidth);
    res = 1.f-(x1+x2-2*xCoord)/binWidth*0.5;
    f2 = f2*res*(x2-x1);
    res = 1.f-(x2-xCoord+binWidth)/binWidth*0.5;
    f3 = f3*res*(binWidth-x2+xCoord);
    res = (f1+f2+f3)/binWidth*0.5;
  }
  return res;
}
static double BoxFilter(double f1,double f2,double f3,double x1,double x2,double xCoord)
{
  double binWidth = x2-x1;
  double res = 0.f;
  if(xCoord-binWidth>x1){
    if(xCoord+binWidth<x2) return f2;
    f2 = f2*(res=x2-xCoord+binWidth,binWidth*=2.f);
    f3 = f3*(binWidth-res);
    res = (f2+f3)/binWidth;
  }
  else if(xCoord+binWidth<x2){
    f1 = f1*(res=x1-xCoord+binWidth,binWidth*=2.f);
    f2 = f2*(binWidth-res);
    res = (f1+f2)/binWidth;
  }
  else{
    f1 = f1*(x1-xCoord+binWidth);
    f2 = f2*(x2-x1);
    f3 = f3*(xCoord+binWidth-x2);
    res = (f1+f2+f3)/binWidth*0.5;
  }
  return res;
}

static double alFindPointFor(G4double xCoord,IVector* pVector,G4double dMaxEnergy)
{
  IPoint* pCurrPoint,*pPrevPoint,*pNextPoint;
  for(G4int i=0;i< pVector->nPoints();i++){
    if(AAType != AA_NONE) if(i != 0) pPrevPoint = pCurrPoint;
    pCurrPoint = pVector->point(i);
    if(pCurrPoint->coordinate(0) > xCoord) break;
  }
  if (i <= 1) return 0;
  if(i==pVector->nPoints()) return 1;
  pCurrPoint = pVector->point(--i);
  pPrevPoint = pVector->point(i-1);
  if(AAType == AA_NONE) return pCurrPoint->value();
  if(i==0) return 0;
  pNextPoint = pVector->point(++i);
  if(pCurrPoint->coordinate(0) >= dMaxEnergy) return pCurrPoint->value();
  if(AAType == AA_BOX)
    return BoxFilter(pPrevPoint->value(),pCurrPoint->value(),
		     pNextPoint->value(),pCurrPoint->coordinate(0),
		     pNextPoint->coordinate(0),xCoord);
  return TentFilter(pPrevPoint->value(),pCurrPoint->value(),
		    pNextPoint->value(),pCurrPoint->coordinate(0),
		    pNextPoint->coordinate(0),xCoord);
}

static G4int FindID()
{
  LPHISTONODE Curr;
  for(G4int i=1;i<100;i++){
    Curr = HistoList;
    while(Curr){
      if(Curr->ID == i) break;
      Curr = Curr->pNext;
    }
    if(Curr) continue;
    return i;
  }
  return -1;
}

static LPHISTONODE FindHisto(G4int nHisto,LPHISTONODE* ppPrev)
{
  *ppPrev = NULL;
  LPHISTONODE Curr = HistoList;
  while(Curr){
    if(Curr->ID == nHisto){
      if((bOptimizeHistos)&&(*ppPrev !=NULL)){      //if this node is not first
	if(Curr->cEvents > Curr->pPrev->cEvents){
	  if(Curr->pNext) Curr->pNext->pPrev = Curr->pPrev;
	  Curr->pPrev->pNext = Curr->pNext;
	  Curr->pNext = Curr->pPrev;
	  if(Curr->pPrev->pPrev !=NULL){
	    Curr->pPrev->pPrev->pNext = Curr;
	    LPHISTONODE pTmp = Curr->pPrev->pPrev;
	    Curr->pPrev->pPrev = Curr;
	    Curr->pPrev = pTmp;
	  }
	  else{
	    HistoList = Curr;
	    Curr->pPrev->pPrev = Curr;
	    Curr->pPrev = NULL;
	  }
	}
      }
      return Curr;
    }
    *ppPrev = Curr;
    Curr = Curr->pNext;
  }
  return Curr;
}

static void ScaleHisto(LPHISTONODE pCurr,G4double val)
{
  G4int i;
  IVector* pVect = pVectorFactory->from1D(pCurr->pHisto);
  pCurr->pHisto->reset();
  bOptimizeHistos=false;
  pCurr->fScale *= val;
  for(i=0;i<pVect->nPoints();i++){
    IPoint* pCurrPoint = pVect->point(i);
    alFill(pCurr->ID,pCurrPoint->coordinate(0),pCurrPoint->value()*val,1,false);
  }
  bOptimizeHistos=true;
  delete pVect;
}

  static unsigned int FindOption(G4String szKey)
{
  char tmp[4];
  strncpy(tmp,szKey.data(),3);
  tmp[3] = 0;
  if(strcasecmp(tmp,"sty")==0) return 10;
  if(strcasecmp(tmp,"gri")==0) return 0;
  if(strcasecmp(tmp,"dvy")==0) return 1;
  if(strcasecmp(tmp,"dvx")==0) return 2;
  if(strcasecmp(tmp,"sta")==0) return 3;
  if(strcasecmp(tmp,"tic")==0) return 4;
  if(strcasecmp(tmp,"box")==0) return 5;
  if(strcasecmp(tmp,"lin")==0){
    if(toupper(szKey.data()[3])=='X') return 6;
    if(toupper(szKey.data()[3])=='Y') return 7;
  }
  if(strcasecmp(tmp,"ast")==0) return 8;
  if(strcasecmp(tmp,"min")==0){
    if(toupper(szKey.data()[3])=='X') return 11;
    if(toupper(szKey.data()[3])=='Y') return 12;
  }
  return ((unsigned int )-1);
}

static unsigned int FindColor(G4String szKey)
{
  if(strcasecmp(szKey.data(),"AxisColor")==0) return 0;
  if(strcasecmp(szKey.data(),"TitleColor")==0) return 1;
  if(strcasecmp(szKey.data(),"AxTitColor")==0) return 2;
  if(strcasecmp(szKey.data(),"HistoColor")==0) return 3;
  if(strcasecmp(szKey.data(),"FuncColor")==0) return 4;
  return ((unsigned int)-1);
}

static unsigned int FindAxis(G4String szKey)
{
  if(strcasecmp(szKey.data(),"xTitle")==0) return 0;
  if(strcasecmp(szKey.data(),"YTitle")==0) return 1;
  if(strcasecmp(szKey.data(),"Title")==0) return 2;
  return ((unsigned int)-1);
}

static LPHISTOSET FindSetNode(LPHISTOPT pOpt,G4String szKey,LPHISTOSET* pPrev)
{
  LPHISTOSET pSet = pOpt->pSetOptions;
  if(pPrev != NULL) *pPrev = NULL;
  while(pSet){
    if(strcasecmp(pSet->String,szKey.data())==0) return pSet;
    if(pPrev != NULL) *pPrev = pSet;
    pSet = pSet->pNext;
  }
  return pSet;
}

static void CalculateZones(unsigned int *xDim,unsigned int* yDim)
{
  float temp = sqrt(nZones);
  *xDim = (unsigned)floor(temp+0.5);
  if ((*xDim)<=0){
    *xDim = 0;
    *yDim = 0;
    return;
  }
  *yDim = (nZones%(*xDim) == 0) ? nZones/(*xDim) : nZones/(*xDim)+1;
}

static char XAxisStr[5] = "XCOL";
static char YAxisStr[5] = "YCOL";
static char HistoStr[5] = "HCOL";
static char FunStr[5] = "FCOL";

static void SetOptions(LPHISTOPT pOpt)
{
  char options[10][5];
  char i;
  char res;
  float tmp;
  LPHISTOSET pCurr = pOpt->pSetOptions;
  for(i=0;i<=9;i++){
    res = (pOpt->cFlags & (1<<i));
    if(res==0)
      strncpy(options[i],NegOptionsArray[i],5);
    else
      strncpy(options[i],OptionsArray[i],5);
  }
  HPLOPT(options,9);
  for(;pCurr != NULL;pCurr = pCurr->pNext)
    HPLSET(pCurr->String,pCurr->value);
  tmp = (float)pOpt->Colors.AxisColor;
  if(tmp!=0.){
    HPLSET(XAxisStr,tmp);
    HPLSET(YAxisStr,tmp);
  }
  tmp = (float)pOpt->Colors.HistoColor;
  if(tmp!=0.)
    HPLSET(HistoStr,tmp);
  tmp = (float)pOpt->Colors.FuncColor;
  if(tmp!=0.)
    HPLSET(FunStr,tmp);
}

static void Plot(LPHISTONODE pCurr)
{
  char temp[3] = " ";
  // char temp1[10]="SHISTO";
  char temp2[6] = "L";
  //  char tmp3[10]="NDVY";
  char tmpOver[10]="S";
  char dModStr[5] = "DMOD";
  float *XArr,*YArr,*XErr,*YErr;
  float *XInt,*YInt;
  unsigned nBins;
  int statType = STAT_DASHLINE;
  int i,j,i1;
  //  float fTmp=1010.f;
  float minX,minY,maxX,maxY;
  float cHistos=1.;
  //strncpy(&(temp[2]),pCurr->pOptions->HistoStyle,2);
  if((pCurr->pOptions->cFlags&STA_BIT))
    AddStat(STAT_INITIALIZE,NULL);
  nBins = pCurr->pHisto->xAxis()->bins();
  XArr = (float*)malloc(sizeof(float)*(4*nBins+6*(nBins+2)));
  YArr = XArr+nBins;
  XErr = YArr+nBins;
  YErr = XErr+nBins;
  XInt = YErr+nBins;
  YInt = XInt+3*(nBins+2);
  HREBIN(pCurr->ID,XArr,YArr,XErr,YErr,nBins,1,nBins);
  if(pCurr->pOptions->maxX==0){
    minX = XArr[0]-XErr[0];
    maxX = XArr[nBins-1]+XErr[nBins-1];
  }
  else{
    minX = pCurr->pOptions->minX;
    maxX = pCurr->pOptions->maxX;
  }
  if(pCurr->pOptions->maxY==0){
    maxY = 0.f;
    minY = 1e+12f;
    for(i=0;i<(signed)nBins;i++){
      if(minY > YArr[i]) minY = YArr[i];
      if(maxY < YArr[i]) maxY = YArr[i];
    }
    minY *= 0.5f;
    maxY *= 1.5f;
  }
  else{
    minY = pCurr->pOptions->minY;
    maxY = pCurr->pOptions->maxY;
  }
  for(i=0;i<(signed)nBins;i++) if(YArr[i]>1e-5) break;
  for(j=nBins-1;j>0;j--) if(YArr[j]>1e-5) break;
  for(i1=i;i<=j;i++){
    XInt[3*(i-i1)] = XArr[i]-XErr[i];
    YInt[3*(i-i1)] = (YArr[i]==0) ? 1e-12 :YArr[i];
    XInt[3*(i-i1)+1] = XArr[i]+XErr[i];
    YInt[3*(i-i1)+1] = (YArr[i]==0) ? 1e-12 :YArr[i];
    XInt[3*(i-i1)+2] = XArr[i]+XErr[i];
    if(i+1<(signed)nBins) YInt[3*(i-i1)+2] = (YArr[i+1]==0) ? 1e-12:YArr[i+1];
    else YInt[3*(i-i1)+2] = (YArr[i]==0) ? 1e-12:YArr[i];
  }
  strncat(tmpOver,pCurr->pOptions->HistoStyle,8);
  HPLFRA(minX,maxX,minY,maxY,temp);
  HPLSET(dModStr,cHistos);
  //HPLOT(pCurr->ID,tmpOver,temp,0);
  HPLINE(XInt,YInt,3*(j-i1+1),temp);
  free(XArr);
  cHistos++;
  if((pCurr->pOptions->cFlags&STA_BIT))
    AddStat(STAT_SOLIDLINE,"Simulation");
  //  int nErrLen;
  //  float* dXArray,*dYArray,*dXErr,*dYErr;
  IAxis* pAxis;
  int isym = 0;
  float usiz = 0;
  for(int j1=0;j1<pCurr->cAdditionalHistos;j1++){
    pAxis = pCurr->pAdditionalHistos[j1]->xAxis();
    nBins = pAxis->bins();
    XArr = (float*)malloc(sizeof(float)*(4*nBins+6*(nBins+2)));
    YArr = XArr+nBins;
    XErr = YArr+nBins;
    YErr = XErr+nBins;
    XInt = YErr+nBins;
    YInt = XInt + 3*(nBins+2);
    HREBIN(pCurr->idAdditionalHistos[j1],XArr,YArr,XErr,YErr,nBins,1,nBins);
    for(i=0;i<(signed)nBins;i++) if(YArr[i]>1e-5) break;
    for(j=nBins-1;j>0;j--) if(YArr[j]>1e-5) break;
    for(i1=i;i<=j;i++){
      XInt[3*(i-i1)] = XArr[i]-XErr[i];
      YInt[3*(i-i1)] = (YArr[i]==0) ? 1e-12 :YArr[i];
      XInt[3*(i-i1)+1] = XArr[i]+XErr[i];
      YInt[3*(i-i1)+1] = (YArr[i]==0) ? 1e-12 :YArr[i];
      XInt[3*(i-i1)+2] = XArr[i]+XErr[i];
      if(i+1<(signed)nBins) YInt[3*(i-i1)+2] = (YArr[i+1]==0) ? 1e-12:YArr[i+1];
      else YInt[3*(i-i1)+2] = (YArr[i]==0) ? 1e-12:YArr[i];
    }
    HPLSET(dModStr,cHistos);
    HPLINE(XInt,YInt,3*(j-i1+1),temp);
    //    HPLOT(pCurr->idAdditionalHistos[j],temp1,temp,0);
    if((++cHistos)>5) cHistos -= 5;
    for(i=0;i<nBins;i++) XErr[i]=0.f; //memset(dXErr,0,sizeof(float)*nErrLen);
    HPLERR(XArr,YArr,XErr,YErr,nBins,temp2,isym,usiz);
    free(XArr);
    if((pCurr->pOptions->cFlags&STA_BIT)){
      AddStat(statType,(char*)pCurr->pAdditionalHistos[j1]->title().c_str());
      if(++statType > STAT_DOTDOTLINE) statType = STAT_DASHLINE;
    }
  }
  cHistos=1;
  HPLSET(dModStr,cHistos);
  nBins = pCurr->pBinErrors->nPoints();
  XArr = (float*)malloc(sizeof(float)*4*nBins);
  YArr = XArr+nBins;
  XErr = YArr+nBins;
  YErr = XErr+nBins;
  pAxis = pCurr->pHisto->xAxis();
  IPoint* pPoint,*pPoint1;
  for(i=0;i<nBins;i++){
    XArr[i] = pAxis->binCentre(i);
    XErr[i] = pAxis->binCentre(i)-pAxis->binLowerEdge(i);
    pPoint = pCurr->pBinErrors->point(i);
    pPoint1 = pCurr->pBinEntries->point(i);
    YArr[i] = pCurr->pHisto->binHeight(i);
    if(pPoint->value()!=0)
      YErr[i] = YArr[i]*sqrt(pPoint->value())/pPoint->value();
    else{
      if((pCurr->pOptions->cFlags & LINY_BIT))
	YErr[i] = 0.f;
      else
	YErr[i] = 1.f;
    }
  }
  HPLERR(XArr,YArr,XErr,YErr,nBins,temp,isym,usiz);
  free(XArr);
}  

static void AddStat(unsigned char bStatType,char* szString)
{
  static int statLine=0;
  static float LineHeight=0.35;
  static float StartX,StartY;
  float CurrX = StartX,CurrY=StartY;
  float X1,Y1,X2,Y2;
  static float TextSize = 0.3;
  float LineCoordsX[2];
  float LineCoordsY[2];
  char cmString[2] = "C";
  static int IOPT = -1;
  static float ANGLE = 0.f;
  static float MaxSiz = 0.f;
  char szTmp[32];
  CurrY -= LineHeight*statLine;
  switch(bStatType){
  case STAT_INITIALIZE:
    HPLGIV(X1,Y1,X2,Y2);
    statLine = 0;
    StartX = X1;
    StartX += (X2-X1)/20.f;
    StartY = Y2-LineHeight;
    return;
  case STAT_SOLIDLINE:
    CurrY -= LineHeight*statLine;
    LineCoordsX[0] = StartX+0.1;
    LineCoordsY[0] = StartY+LineHeight/2.f;
    LineCoordsX[1] = StartX+0.1f+TextSize*5.f;
    LineCoordsY[1] = StartY+LineHeight/2.f;
    HPLINE(LineCoordsX,LineCoordsY,2,cmString);
    CurrX += 0.2f+TextSize*5.f;
    break;
  case STAT_DASHLINE:
    szTmp[0] = szTmp[2] = szTmp[4] = '-';
    szTmp[1] = szTmp[3] = ' ';
    szTmp[5] = 0;
    CurrX += 0.1f;
    HPLSOF(CurrX,CurrY,szTmp,TextSize,ANGLE,MaxSiz,IOPT);
    CurrX += 0.2f+TextSize*5.f;
    break;
  case STAT_DASHDOTLINE:
    szTmp[0] = szTmp[2] = szTmp[4] = '-';
    szTmp[1] = szTmp[3] = '.';
    szTmp[5] = 0;
    CurrX += 0.1f;
    HPLSOF(CurrX,CurrY,szTmp,TextSize,ANGLE,MaxSiz,IOPT);
    CurrX += 0.02f+TextSize*5.f;
    break;
  case STAT_DOTLINE:
    szTmp[0]=szTmp[2] = szTmp[4] = '.';
    szTmp[1] = szTmp[3] = ' ';
    szTmp[5] = 0;
    CurrX += 0.1f;
    HPLSOF(CurrX,CurrY,szTmp,TextSize,ANGLE,MaxSiz,IOPT);
    CurrX += 0.2f+TextSize*5.f;
    break;
  case STAT_DOTDOTLINE:
    szTmp[0] = szTmp[1] = szTmp[2] = '.';
    szTmp[3] = 0;
    CurrX += 0.01f+TextSize;
    HPLSOF(CurrX,CurrY,szTmp,TextSize,ANGLE,MaxSiz,IOPT);
    CurrX += 0.02f+TextSize+5.f;
  }
  strncpy(szTmp,szString,31);
  szTmp[31]=0;
  HPLSOF(CurrX,CurrY,szTmp,TextSize,ANGLE,MaxSiz,IOPT);
  statLine++;
}

static void PlotHisto(unsigned int ID,bool bRecalculateZones)
{
  unsigned int x,y;
  char* xTitle,*yTitle;
  char temp[3] = " ";
  LPHISTONODE pCurr,pDummy;
  bOptimizeHistos=false;
  pCurr = FindHisto(ID,&pDummy);
  SetOptions(pCurr->pOptions);
  if(bRecalculateZones){
    if(nZones==1){
      HPLINT(1);
      SetOptions(pCurr->pOptions);
    }
    CalculateZones(&x,&y);
    HPLZON(x,y,1,temp);
  }
  bOptimizeHistos=true;
  Plot(pCurr);
  if(pCurr->pOptions->XAxis==NULL) xTitle = temp;
  else xTitle = pCurr->pOptions->XAxis;
  if(pCurr->pOptions->YAxis==NULL) yTitle = temp;
  else yTitle = pCurr->pOptions->YAxis;
  HPLLGD(xTitle,yTitle,temp,222,temp);
}

static void alAddData(IVector* pData)
{
  float* x,*y;
  unsigned Points = pData->nPoints();
  char temp[3] = " ";
  x = new float[Points];
  y = new float[Points];
  for(unsigned i=0;i<Points;i++){
    x[i] = pData->point(i)->coordinate(0);
    y[i] = pData->point(i)->value();
  }
  HPLFUN(x,y,Points,temp);
  delete[] x;
  delete[] y;
}
static void alShiftZones(unsigned StartFrom)
{
  char temp[3] = " ";
  unsigned x,y;
  LPHISTONODE pCurr = HistoList;
  while(pCurr){
    if(pCurr->pViewer >= StartFrom) pCurr->pViewer--;
    pCurr = pCurr->pNext;
  }
  nZones--;
  if(nZones> 0){
    CalculateZones(&x,&y);
    HPLZON(x,y,1,temp);
  }
  else{
    HPLEND();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  exported functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void alSetTitle(char* title)
{
  if(g_szGlobalTitle) delete g_szGlobalTitle;
  g_szGlobalTitle = new char[strlen(title)+1];
  strcpy(g_szGlobalTitle,title);
}

void alSetAA(int AA)
{
  if(AA<=AA_TENT)
    AAType = AA;
}

IVectorFactory* alGetVectorFactory()
{
  return pVectorFactory;
}

void alInit()
{
  HistoList = NULL;
  pHistoFactory = createIHistoFactory();
  if(pHistoFactory==NULL){
    G4cerr<<"Error: cannot create Histogram Factory object"<<G4endl;
    exit(1);
  }
  pVectorFactory = createIVectorFactory();
  if(pVectorFactory==NULL){
    G4cerr<<"Error: cannot create Vector Factory object"<<G4endl;
    exit(1);
  }
  szStoreName = NULL;
  bViewerIsShown = false;
  bOptimizeHistos = true;
}

G4int alAddHisto(G4String szName,G4int nBins,G4double dMin,G4double dMax)
{
  LPHISTONODE NewNode = new HISTONODE;
  G4int HistoID = FindID();
  if(NewNode==NULL) return -1;
  if(HistoID==-1){
    delete NewNode;
    return -1;
  }
  NewNode->ID = HistoID;
  NewNode->cEvents = 0;
  NewNode->pOptions = new HISTOPT;
  NewNode->pOptions->XAxis = NewNode->pOptions->YAxis = NULL;
  NewNode->pOptions->cFlags = DEFAULT_OPTIONS;
  NewNode->pOptions->Colors.AxisColor = 0;
  NewNode->pOptions->Colors.TitleColor = 0;
  NewNode->pOptions->Colors.AxisTitleColor = 0;
  NewNode->pOptions->Colors.HistoColor = 0;
  NewNode->pOptions->Colors.FuncColor = 0;
  NewNode->pOptions->minX = 0.f;
  NewNode->pOptions->minY = 0.f;
  NewNode->pOptions->maxX = 0.f;
  NewNode->pOptions->maxY = 0.f;
  NewNode->pOptions->HistoStyle[0]=' ';
  NewNode->pOptions->HistoStyle[1]=0;
  NewNode->pOptions->pSetOptions = NULL;
  NewNode->pHisto = pHistoFactory->create1D(szName.c_str(),nBins,dMin,dMax,HistoID);
  NewNode->pBinErrors = pVectorFactory->create();
  NewNode->pBinEntries = pVectorFactory->create();
  if(NewNode->pHisto==NULL){
    delete NewNode->pOptions;
    delete NewNode;
    return -1;
  }
  IAxis* pAxis = NewNode->pHisto->xAxis();
  for(int i=0;i<pAxis->bins();i++){
    IPoint* pPoint = pVectorFactory->createPoint();
    pPoint->setDimension(1);
    IScalar* pScalar = pPoint->coordScalar(0);
    pScalar->setValue(pAxis->binCentre(i));
    pPoint->vScalar()->setValue(0);
    NewNode->pBinErrors->addPoint(pPoint);
    pPoint = pVectorFactory->createPoint();
    pPoint->setDimension(1);
    pScalar = pPoint->coordScalar(0);
    pScalar->setValue(pAxis->binCentre(i));
    pPoint->vScalar()->setValue(0);
    NewNode->pBinEntries->addPoint(pPoint);
  }
  NewNode->pAdditionalHistos = NULL;
  NewNode->cAdditionalHistos = 0;
  NewNode->idAdditionalHistos = NULL;
  NewNode->fScale = 1.f;
  NewNode->pDataSet = NULL;
  NewNode->pViewer = 0;
  NewNode->pNext = HistoList;
  NewNode->pPrev = NULL;
  if(HistoList) HistoList->pPrev = NewNode;
  HistoList = NewNode;
  return HistoID;
}
  
void alDeleteHisto(G4int nHisto)
{
  LPHISTONODE pDelNode,pPrevNode;
  bOptimizeHistos = false;
  pDelNode = FindHisto(nHisto,&pPrevNode);
  if(pDelNode==NULL){
    G4cerr<<"Warning: Histogram: "<<nHisto<<" is not in use"<<G4endl;
    bOptimizeHistos = true;
    return;
  }
  if(pPrevNode==NULL){  //deleting first node
    HistoList = pDelNode->pNext;
    if(HistoList) HistoList->pPrev = NULL;
  }
  else{
    if(pDelNode->pNext) pDelNode->pNext->pPrev = pPrevNode;
    pPrevNode->pNext = pDelNode->pNext;
  }
  if(pDelNode->pViewer)  alHideHisto(pDelNode->ID);
  if(pDelNode->pDataSet) delete pDelNode->pDataSet;
  if(pDelNode->pOptions->pSetOptions){
    LPHISTOSET pCurr = pDelNode->pOptions->pSetOptions;
    while(pCurr){
      pDelNode->pOptions->pSetOptions = pCurr->pNext;
      delete pCurr->String;
      delete pCurr;
      pCurr = pDelNode->pOptions->pSetOptions;
    }
  }
  for(int i=0;i<pDelNode->cAdditionalHistos;i++){
    delete pDelNode->pAdditionalHistos[i];
  }
  if(pDelNode->cAdditionalHistos){
    free(pDelNode->pAdditionalHistos);
    free(pDelNode->idAdditionalHistos);
  }
  delete pDelNode->pOptions;
  delete pDelNode->pHisto;
  delete pDelNode->pBinErrors;
  delete pDelNode->pBinEntries;
  delete pDelNode;
  bOptimizeHistos=true;
}

void alShowHisto(G4int nHisto,G4int Style)
{
  LPHISTONODE pShowNode,pDummy;
  pShowNode = FindHisto(nHisto,&pDummy);
  if(pShowNode==NULL){
    G4cerr<<"Warning: Histogram ID: "<<nHisto<<" does not point to histogram"<<G4endl;
    return;
  }
  if(pShowNode->pViewer==0){
    if(nZones>=4){
      G4cout<<"Due to lack of room only 4 histograms might be shown at ones"<<G4endl;
      return;
    }
    pShowNode->pViewer = nZones+1;
    alRefreshHisto(true);
  }
  else
    alRefreshHisto(false);
}

void alHideHisto(G4int nHisto)
{
  LPHISTONODE pHideNode,pDummy;
  bOptimizeHistos = false;
  pHideNode = FindHisto(nHisto,&pDummy);
  if(pHideNode==NULL){
    G4cerr<<"Warning: Histogram ID: "<<nHisto<<" does not point to histogram"<<G4endl;
    bOptimizeHistos = true;
    return;
  }
  if(pHideNode->pViewer){
    alShiftZones(pHideNode->pViewer+1);
    pHideNode->pViewer = 0;
    alRefreshHisto(false);
  }
  bOptimizeHistos=true;
}

void alMakePS(G4int nHisto, G4String file)
{
  int istat;
  char temp[5] = "NEW";
  char temp2[3] = " ";
  char filename[1024];
  LPHISTONODE pPSNode,pDummy;
  pPSNode = FindHisto(nHisto,&pDummy);
  if(pPSNode==NULL){
    G4cerr<<"Warning: Histogram ID: "<<nHisto<<" does not point to histogram"<<G4endl;
    return;
  }
  strncpy(filename,file.c_str(),1023);
  filename[1023] = 0;
  if(nZones==0){ HPLINT(0);}
  else HPLZON(1,1,1,temp2);
  KUOPEN(77,filename,temp,istat);
  IGMETA(-77,-111);
  PlotHisto(pPSNode->ID,false);
  if(pPSNode->pDataSet!=NULL) alAddData(pPSNode->pDataSet);
  IGMETA(999,0);
  KUCLOS(77,temp2,1);
  if(nZones==0){HPLEND();}
  else{
  }
}

void alOpenStorage(G4String szName)
{
  pHistoFactory->selectStore(szName);
  if(szStoreName) delete szStoreName;
  if(szName.data()!=NULL){
    szStoreName = new char[strlen(szName.data())+1];
    strcpy(szStoreName,szName.data());
  }
}

void alSaveHisto(G4int nHisto)
{
  LPHISTONODE pSaveNode;
  LPHISTONODE pDummy;
  pSaveNode = FindHisto(nHisto,&pDummy);
  if(pSaveNode==NULL){
    G4cerr<<"Warning: Histogram ID: "<<nHisto<<" does not point to histogram"<<G4endl;
    return;
  }
  pHistoFactory->store1D(pSaveNode->pHisto);
}
void alRemoveFromStore(G4int nHisto)
{
  char szHistoLabel[32];
  int i=0;
  LPHISTONODE pRemoveNode;
  LPHISTONODE pDummy;
  pRemoveNode = FindHisto(nHisto,&pDummy);
  if(pRemoveNode==NULL){
    G4cerr<<"Warning: HistogramID: "<<nHisto<<" dos not point to histogram"<<G4endl;
    return;
  }
  szHistoLabel[i]=0;
  if(nHisto<10){
    szHistoLabel[0]=(char)nHisto+'0';
    szHistoLabel[1]=0;
  }
  else{
    while(nHisto/((i+1)*10)!=0){
      unsigned char c = (char)((nHisto%((i+1)*10))/(i*10));
      szHistoLabel[i]=c+'0';
      szHistoLabel[i+1]=0;
      i++;
    }
    int j=0;
    while(i>j+1){
      szHistoLabel[j] ^= szHistoLabel[i-1];
      szHistoLabel[i-1] ^= szHistoLabel[j];
      szHistoLabel[j] ^= szHistoLabel[i-1];
      j++;
      i--;
    }
  }
  pHistoFactory->scratchHisto(szHistoLabel);
}
void alAddData(G4int nHisto,IVector* pVector)
{
  LPHISTONODE pAddData;
  LPHISTONODE pDummy;
  pAddData = FindHisto(nHisto,&pDummy);
  if(pAddData){
    if(pAddData->pDataSet) delete pAddData->pDataSet;
    pAddData->pDataSet = pVector;
  }
}
//-------------------------------------------------------------------
void alAdditionalHisto(G4int nHist,const char* Storage,const char* label)
{
  LPHISTONODE pCurr,pDummy;
  IHistogram1D* pTmp=NULL;
  unsigned  HistoID;
  //  char* pLast;
  HistoID = atol(label);
  //  if(*pLast!=0) return;
  pCurr = FindHisto(nHist,&pDummy);
  if(pCurr==NULL) return;
  if(Storage!=NULL){
    pHistoFactory->selectStore(Storage);
  }
  pTmp = pHistoFactory->load1D(label);
  if(pTmp==NULL){
    if((Storage!=NULL)&&(szStoreName!=NULL))
      pHistoFactory->selectStore(szStoreName);
    return;
  }
  if(pCurr->cAdditionalHistos==0){
    pCurr->pAdditionalHistos = (IHistogram1D**)malloc(sizeof(void*)*5);
    pCurr->idAdditionalHistos = (unsigned*)malloc(sizeof(unsigned)*5);
  }
  else if((pCurr->cAdditionalHistos+1)%5==0){
    pCurr->pAdditionalHistos = (IHistogram1D**)realloc(pCurr->pAdditionalHistos,sizeof(void*)*(pCurr->cAdditionalHistos/5 + 1)*5);
    pCurr->idAdditionalHistos = (unsigned*)realloc(pCurr->idAdditionalHistos,sizeof(unsigned)*(pCurr->cAdditionalHistos/5+1)*5);
  }
  pCurr->pAdditionalHistos[pCurr->cAdditionalHistos]=pTmp;
  pCurr->idAdditionalHistos[pCurr->cAdditionalHistos]=HistoID;
  pCurr->cAdditionalHistos++;
  if(Storage!=NULL){
    if(szStoreName!=NULL)
      pHistoFactory->selectStore(szStoreName);
  }
}
//-----------------------------------------------------------------------
void alAddDataFromFile(G4int nHist,G4String szFileName)
{
  IVector* pVector = pVectorFactory->create();
  szFileName.strip(2);
  pVector->fromAscii(szFileName.data());
  alAddData(nHist,pVector);
}

void ShowViewer(G4bool bYes)
{
  bViewerIsShown = bYes;
}

void alSetProperty(G4int nHisto,G4String szKey,G4String szValue)
{
  LPHISTONODE pSetProp;
  LPHISTONODE pDummy;
  unsigned int Value;
  pSetProp = FindHisto(nHisto,&pDummy);
  if(pSetProp==NULL){
    G4cout<<"Histogram with id "<<nHisto<<" does not exist"<<G4endl;
    return;
  }
  unsigned int OptionID = FindOption(szKey);
  if(OptionID == ((unsigned int)-1)){
    OptionID = FindColor(szKey);
    if(OptionID==((unsigned int)-1)){
      OptionID = FindAxis(szKey);
      if(OptionID==((unsigned int)-1)){
	char* ptmp;
	float Val = strtod(szValue.data(),&ptmp);
	if((*ptmp!=' ')&&(*ptmp != 0)){
	  G4cout<<"Unable to determine what is: "<<szKey.data()<<G4endl;
	  return;
	}
	LPHISTOSET pNew;
	szKey.toUpper();
	if ((pNew = FindSetNode(pSetProp->pOptions,szKey,NULL)) == NULL){
	  pNew = new HISTOSET;
	  pNew->String = new char[strlen(szKey.data())+1];
	  strcpy(pNew->String,szKey.data());
	  pNew->pNext = pSetProp->pOptions->pSetOptions;
	  pSetProp->pOptions->pSetOptions = pNew;
	}
	pNew->value = Val;
      }
      else{
	char* tmp = new char[strlen(szValue.data())+1];
	strcpy(tmp,szValue.data());
	if(OptionID==0){
	  if(pSetProp->pOptions->XAxis) delete pSetProp->pOptions->XAxis;
	  pSetProp->pOptions->XAxis = tmp;
	}
	else if(OptionID==1){
	  if(pSetProp->pOptions->YAxis) delete pSetProp->pOptions->YAxis;
	  pSetProp->pOptions->YAxis = tmp;
	}
	else{
	  if(g_szGlobalTitle) delete g_szGlobalTitle;
	  g_szGlobalTitle = tmp;
	}
      }
    }
    else{
      Value = atoi(szValue);
      switch(OptionID){
      case 0: pSetProp->pOptions->Colors.AxisColor= Value; break;
      case 1: pSetProp->pOptions->Colors.TitleColor = Value; break;
      case 2: pSetProp->pOptions->Colors.AxisTitleColor = Value; break;
      case 3: pSetProp->pOptions->Colors.HistoColor = Value; break;
      case 4: pSetProp->pOptions->Colors.FuncColor = Value;
      }  
    }
  }
  else{
    if(OptionID == 10){
      char C = szValue.data()[0];
      if((C == 'H')||(C == 'L')||(C=='*')||(C=='C')||(C=='B'))
	pSetProp->pOptions->HistoStyle[0] = toupper(C);
    }
    else if(OptionID == 11){
      G4String max;
      szValue.strip(2);
      int pos = szValue.find(' ');
      max = szValue;
      max.remove(0,pos+1);
      szValue.remove(pos,szValue.length()-pos);
      float fMin = atof(szValue.data());
      float fMax = atof(max);
      pSetProp->pOptions->minX = fMin;
      pSetProp->pOptions->maxX = fMax;
    }
    else if(OptionID == 12){
      G4String max;
      szValue.strip(2);
      int pos = szValue.find(' ');
      max = szValue;
      max.remove(0,pos+1);
      szValue.remove(pos,szValue.length()-pos);
      float fMin = atof(szValue);
      float fMax = atof(max);
      pSetProp->pOptions->minY = fMin;
      pSetProp->pOptions->maxY = fMax;
    }
    else{
      Value = atoi(szValue);
      if(Value==0)
	pSetProp->pOptions->cFlags &= ~OptArray[OptionID];
      else
	pSetProp->pOptions->cFlags |= OptArray[OptionID];
    }
  }
  alRefreshHisto(false);
}

void alReset()
{
  LPHISTONODE pCurr = HistoList;
  int index,i;
  while(pCurr){
    pCurr->pHisto->reset();
    pCurr->cEvents = 0;
    index = pCurr->pBinErrors->nPoints();
    for(i=0;i<index;i++){
      pCurr->pBinErrors->point(i)->vScalar()->setValue(0);
      pCurr->pBinEntries->point(i)->vScalar()->setValue(0);
    }
    if(pCurr->cAdditionalHistos){
      for(int i=0;i<pCurr->cAdditionalHistos;i++)
	delete pCurr->pAdditionalHistos[i];
      free(pCurr->pAdditionalHistos);
      free(pCurr->idAdditionalHistos);
      pCurr->cAdditionalHistos=0;
      pCurr->idAdditionalHistos = NULL;
      pCurr->pAdditionalHistos = NULL;
    }
    pCurr = pCurr->pNext;
  }
  alRefreshHisto(false);
}
void alClear()
{
  int i;
  LPHISTONODE pCurr = HistoList;
  LPHISTONODE pTmp;
  while(pCurr){
    if(pCurr->pDataSet) delete pCurr->pDataSet;
    if(pCurr->pBinErrors) delete pCurr->pBinErrors;
    if(pCurr->pBinEntries) delete pCurr->pBinEntries;
    if(pCurr->pOptions){
      LPHISTOSET pCurrOption = pCurr->pOptions->pSetOptions;
      LPHISTOSET pNextOption;
      while(pCurrOption){
	pNextOption = pCurrOption->pNext;
	delete pCurrOption;
	pCurrOption = pNextOption;
      }
      if(pCurr->pOptions->XAxis) delete pCurr->pOptions->XAxis;
      if(pCurr->pOptions->YAxis) delete pCurr->pOptions->YAxis;
      delete pCurr->pOptions;
    }
    if(pCurr->cAdditionalHistos){
      for(i=0;i<pCurr->cAdditionalHistos;i++)
	delete pCurr->pAdditionalHistos[i];
    }
    if(pCurr->pAdditionalHistos) free(pCurr->pAdditionalHistos);
    if(pCurr->idAdditionalHistos) free(pCurr->idAdditionalHistos);
    pTmp = pCurr->pNext;
    delete pCurr->pHisto;
    delete pCurr;
    pCurr = pTmp;
  }
  delete pVectorFactory;
  delete pHistoFactory;
}

void alFill(G4int nHisto,G4double dValue,G4double weight,G4double EvWeight,bool bCalcErrors)
{
  LPHISTONODE pFillNode;
  LPHISTONODE pDummy;
  pFillNode = FindHisto(nHisto,&pDummy);
  if(pFillNode){
    pFillNode->pHisto->fill(dValue,weight);
    int index=0;
    if(bCalcErrors){
      for(;index<pFillNode->pHisto->xAxis()->bins();index++)
	if(pFillNode->pHisto->xAxis()->binLowerEdge(index)>dValue) break;
      if((index-1 < pFillNode->pHisto->xAxis()->bins())&&(index-1>0)){
	IPoint* pPoint = pFillNode->pBinErrors->point(index-1);
	pPoint->vScalar()->setValue(pPoint->value()+EvWeight*EvWeight);
	pPoint = pFillNode->pBinEntries->point(index-1);
	pPoint->vScalar()->setValue(pPoint->value()+1);
      }
    }
  }
}

void alScaleAllHistos(G4double val)
{
  LPHISTONODE pCurr = HistoList;
  while(pCurr){
    ScaleHisto(pCurr,val);
    pCurr = pCurr->pNext;
  }
}

void alSubstract(G4int nHisto,G4int nHistoOut,G4double dMaxEnergy,G4double dMinEnergy)
{
  LPHISTONODE pCurr,pDummy;
  LPHISTONODE pOut;
  G4double value;
  pCurr = FindHisto(nHisto,&pDummy);
  pOut = FindHisto(nHistoOut,&pDummy);
  if((pCurr)&&(pOut)){
    if(pCurr->pDataSet != NULL){
      pOut->pHisto->reset();
      IVector* pVect1 = pVectorFactory->from1D(pCurr->pHisto);
      for(G4int i=0;i<pCurr->pDataSet->nPoints();i++){
	IPoint* pPoint1/*,*pPoint2*/;
	pPoint1 = pCurr->pDataSet->point(i);
	if(pPoint1->coordinate(0)<dMinEnergy) continue;
	value = alFindPointFor(pPoint1->coordinate(0),pVect1,dMaxEnergy);
	if((fabs(value)>1e-6)&&(fabs(pPoint1->value())>1e-6)){
	  pOut->pHisto->fill(pPoint1->coordinate(0),(pPoint1->value()-value)/pPoint1->value());
	}
	else{
	  pOut->pHisto->fill(pPoint1->coordinate(0),1.);
	}
      }
      delete pVect1;
    }
  }
}

double alGetMinimumEntries(int LastID)
{
  double ret,minimum=0;
  bool DoNotDeleteMinimum=false;
  LPHISTONODE pCurr = HistoList,pRes=HistoList;
  while(pCurr){
    if(pCurr->ID<=LastID){
      if(pCurr->pHisto){
	ret = pCurr->cEvents;
	if (ret == 0) return (ret);
	if(!DoNotDeleteMinimum){
	  pRes = pCurr;
	  minimum = ret;
	  DoNotDeleteMinimum = true;
	}
	else if(minimum > ret){
	  minimum = ret;
	  pRes = pCurr;
	}
      }
    }
    pCurr = pCurr->pNext;
  }
  return minimum;
}

unsigned alGetEntries(int nHisto)
{
  unsigned ret = 0;
  LPHISTONODE pHisto,pDummy;
  pHisto = FindHisto(nHisto,&pDummy);
  if(pHisto ==NULL) return 0;
  ret = pHisto->pHisto->entries();
  return ret;
}

float alCalculateIntegral(int nHisto)
{
  LPHISTONODE pHisto,pDummy;
  IVector* pVect;
  IPoint* pPoint1,*pPoint2;
  float ret=0.;

  pHisto = FindHisto(nHisto,&pDummy);
  if(pHisto==NULL) return 1.;
  pVect = pVectorFactory->from1D(pHisto->pHisto);
  if(pVect->point(0)->coordinate(0)>10){
    ret = /*pVect->point(0)->coordinate(0)*/(pVect->point(0)->coordinate(0)-10)*pVect->point(0)->value();
  }
  else{
    //delete pVect;
    ret = 0;
  }
  for(int i=0;i<pVect->nPoints()-1;i++){
    pPoint1 = pVect->point(i);
    pPoint2 = pVect->point(i+1);
    if(pPoint1->coordinate(0)>0){
      ret += (pPoint2->coordinate(0) - pPoint1->coordinate(0))*(pPoint2->coordinate(0)*pPoint2->value()+
								pPoint1->coordinate(0)*pPoint1->value())/2;
    }
    else if(pPoint2->coordinate(0)>0){
      ret += (2*pPoint2->value()-pPoint1->value())/(2*(pPoint2->coordinate(0)-pPoint1->coordinate(0)))*(pPoint2->coordinate(0)-10);
    }
  }
  //  delete pVect;
  return ret;
}

void alAddEvent(int nHisto)
{
  LPHISTONODE pHisto,pDummy;
  pHisto = FindHisto(nHisto,&pDummy);
  if(pHisto->cEvents != (unsigned long)-1)pHisto->cEvents++;
}

void alRefreshHisto(bool bAddZone)
{
  LPHISTONODE pCurr;
  unsigned zones = nZones;
  if(bAddZone){
    zones++;
    nZones++;
  }
  bOptimizeHistos = false;
  for(unsigned int i=1;i<=zones;i++){
    pCurr = HistoList;
    while(pCurr){
      if(pCurr->pViewer == i) break;
      pCurr = pCurr->pNext;
    }
    if(pCurr==NULL) break;
    if(((i==1)&&(zones!=1))||(bAddZone)){
      PlotHisto(pCurr->ID,true);
      bAddZone = false;
    }
    else PlotHisto(pCurr->ID,false);
    if(pCurr->pDataSet!=NULL) alAddData(pCurr->pDataSet);
  }
  bOptimizeHistos=true;
}

void alListHistos()
{
  LPHISTONODE pCurr = HistoList;
  while(pCurr){
    G4cout<<"  "<<pCurr->ID<<"        "<<pCurr->pHisto->title()<<"        ";
    if(pCurr->pDataSet)  G4cout<<"Data available        ";
    if(pCurr->pViewer!=0) G4cout<<"Viewed at zone: "<<pCurr->pViewer;
    G4cout<<G4endl;
    pCurr = pCurr->pNext;
  }
}

void alAddTitles(int nHisto,char* XTitle,char* YTitle)
{
  LPHISTONODE pCurr,pDummy;
  pCurr = FindHisto(nHisto,&pDummy);
  if(pCurr->pOptions->XAxis!=NULL){
    delete pCurr->pOptions->XAxis;
    pCurr->pOptions->XAxis = NULL;
  }
  if(XTitle!=NULL){
    pCurr->pOptions->XAxis = new char[strlen(XTitle)+1];
    strcpy(pCurr->pOptions->XAxis,XTitle);
  }
  if(pCurr->pOptions->YAxis){
    delete pCurr->pOptions->YAxis;
    pCurr->pOptions->YAxis = NULL;
  }
  if(YTitle){
    pCurr->pOptions->YAxis = new char[strlen(YTitle)+1];
    strcpy(pCurr->pOptions->YAxis,YTitle);
  }
}
    
void alSetBit(int nHisto, unsigned Option,bool bSet)
{
  LPHISTONODE pCurr,pDummy;
  pCurr = FindHisto(nHisto,&pDummy);
  if(Option > 8) return;
  unsigned mask = 1<<Option;
  if(bSet) pCurr->pOptions->cFlags = (pCurr->pOptions->cFlags | mask);
  else pCurr->pOptions->cFlags = (pCurr->pOptions->cFlags & (~mask));
}

void alReadOptions(unsigned int ID, char *file)
{
  char temp[512];
  char temp1[512];
  FILE* pFile = fopen(file,"rt");
  if(pFile==NULL){
    G4cout<<"Error opening: "<<file<<G4endl;
    return;
  }
  LPHISTONODE pCurr,pDummy;
  pCurr = FindHisto(ID,&pDummy);
  if(pCurr==NULL){
    fclose(pFile);
    G4cout<<"Histogram with id "<<ID<<" does not exist"<<G4endl;
    return;
  }
  while(!feof(pFile)){
    if(fscanf(pFile,"%512s%512s",temp,temp1)<2) break;
    alSetProperty(ID,G4String(temp),G4String(temp1));
  }
  fclose(pFile);
}

void alDeleteOptions(unsigned int ID, G4String szKey)
{
  LPHISTONODE pCurr,pDummy;
  LPHISTOSET pSet,pPrev;
  pCurr = FindHisto(ID,&pDummy);
  if(pCurr==NULL){
    G4cout<<"Hitogram with id "<<ID<<" does not exist"<<G4endl;
    return;
  }
  pSet = FindSetNode(pCurr->pOptions,szKey,&pPrev);
  if(pSet != NULL){
    if(pPrev == NULL)
      pCurr->pOptions->pSetOptions = pSet->pNext;
    else
      pPrev->pNext = pSet->pNext;
    delete pSet->String;
    delete pSet;
  }
}

void alListOptions(unsigned int nHisto)
{
  LPHISTONODE pCurr,pDummy;
  LPHISTOSET pSet;
  pCurr = FindHisto(nHisto,&pDummy);
  if(pCurr==NULL){
    G4cout<<"Histogram with ID "<<nHisto<<" does not exist"<<G4endl;
      return;
  }
  pSet = pCurr->pOptions->pSetOptions;
  G4cout<<"Options for histogram with ID:"<<nHisto<<G4endl;
  while(pSet){
    G4cout<<" Option name "<<pSet->String<<"     = "<<pSet->value<<G4endl;
    pSet = pSet->pNext;
  }
  G4cout<<"That is all"<<G4endl;
}

void alCreateSmooth(G4int nHisto,G4int nTargetHisto,G4double dMaxEnergy)
{
  LPHISTONODE pHisto,pTargetHisto,pDummy;
  pHisto = FindHisto(nHisto,&pDummy);
  pTargetHisto = FindHisto(nTargetHisto,&pDummy);
  if((pHisto==NULL)||(pTargetHisto==NULL)) return;
  if((pHisto->pDataSet!=NULL)&&(pTargetHisto->pDataSet==NULL)){
    IVector* pData = pHisto->pDataSet->deepClone();
    alAddData(nTargetHisto,pData);
  }
  pTargetHisto->pHisto->reset();
  IVector* pTmp = pVectorFactory->from1D(pHisto->pHisto);
  for(G4int i=0;i<pTmp->nPoints()-1;i++){
    IPoint* pPoint1 = pTmp->point(i);
    pTargetHisto->pHisto->fill(pPoint1->coordinate(0),
			       alFindPointFor(pPoint1->coordinate(0),pTmp,
					      dMaxEnergy));
  }
  delete pTmp;
  alRefreshHisto(false);
}

bool alRebin(unsigned nFirstHisto,unsigned nLastHisto,unsigned bias)
{
  LPHISTONODE pHisto = HistoList,pCurr,pDummy;
  unsigned nBins;
  unsigned nPrevBins;
  vector<float> vEdges;
  IVector* pVector;
  float* XArr,*XArr1;
  float* YArr,*YArr1;
  float* XErr,*XErr1;
  float* YErr,*YErr1;
  signed i,j;
  char* title,c=0;
  char number[20];
  float fTmp,fTmp1,fTmp2,fTmp3,fTmp4;
  for(;pHisto;pHisto = pHisto->pNext){
    if((pHisto->ID>=(signed)nFirstHisto)&&(pHisto->ID<=(signed)nLastHisto)){
      pCurr = FindHisto(nFirstHisto+bias,&pDummy);
      if((pCurr!=NULL)&&(pHisto->cAdditionalHistos!=0)){
	//rebinning pCurr with data from pHisto->pAdditionalHistos[0]
	nBins = pHisto->pAdditionalHistos[0]->xAxis()->bins();
	nPrevBins = pCurr->pHisto->xAxis()->bins();
	XArr = new float[nBins*4];
	YArr = XArr+nBins;
	XErr = YArr+nBins;
	YErr = XErr+nBins;
	HREBIN(pHisto->idAdditionalHistos[0],XArr,YArr,XErr,YErr,nBins,1,nBins);
	for(i=0;i<(signed)nBins;i++) vEdges.push_back(XArr[i]-XErr[i]);
	vEdges.push_back(XArr[nBins-1]+XErr[nBins-1]);
	title = (char*)(pCurr->pHisto->title().c_str());
	XArr1 = new float[(nPrevBins+nBins)*4];
	YArr1 = XArr1+nPrevBins;
	XErr1 = YArr1+nPrevBins;
	YErr1 = XErr1+nPrevBins;
	HREBIN(pCurr->ID,XArr1,YArr1,XErr1,YErr1,nPrevBins,1,nPrevBins);
	delete pCurr->pHisto;
	number[0]=0;
	for(i=pCurr->ID,j=0;i!=0;i /= 10,j++){
	  number[j] = (i%10)+'0';
	  number[j+1] = 0;
	}
	for(i=0;i+1<j;i++,j--){
	  number[i] ^= number[j];
	  number[j] ^= number[i];
	  number[i] ^= number[j];
	}
	pCurr->pHisto = pHistoFactory->create1DVar(number,title,vEdges);
	delete XArr;
	XArr = YErr1+nPrevBins;
	YArr = XArr+nBins;
	XErr = YArr+nBins;
	YErr = XErr+nBins;
	HREBIN(pCurr->ID,XArr,YArr,XErr,YErr,nBins,1,nBins);
	bOptimizeHistos=false;
	if(pCurr->pBinErrors){
	  IPoint*pPoint;
	  IScalar* pScalar;
	  pVector = pVectorFactory->create();
	  for(i=0;i<(signed)nBins;i++){
	    pPoint = pVectorFactory->createPoint();
	    pPoint->setDimension(1);
	    pScalar = pPoint->coordScalar(0);
	    pScalar->setValue(pCurr->pHisto->xAxis()->binCentre(i));
	    pPoint->vScalar()->setValue(0);
	    pVector->addPoint(pPoint);
	  }
	}
	IPoint* pPoint1;
	//bez 1 - binovete na novata histograma
	//s 1 - binovete na starata histograma
	// iskame s1 -> bez 1
	for(i=0;i<(signed)nPrevBins;i++){
	  for(j=0;j<(signed)nBins;j++){
	    if(((XArr1[i]-XErr1[i] >= XArr[j]-XErr[j])&&
	       (XArr1[i]-XErr1[j] <= XArr[j]+XErr[j]))||
	       ((XArr1[i]+XErr1[i] >= XArr[j]-XErr[j])&&
		(XArr1[i]+XErr1[i] <= XArr[j]+XErr[j]))){
	      //samo ako levija ili desnija kraj na starata e v bina na novata
	      //clipvame starija bin vyrhu novija
	      fTmp1 = XArr1[i]-XErr1[i];
	      fTmp2 = XArr1[i]+XErr1[i];
	      fTmp3 = XArr[j]-XErr[j];
	      fTmp4 = XArr[j]+XErr[j];
	      //fTmp1->fTmp2 - segmenta na starata
	      //fTmp3->fTmp4 - segmenta na novata
	      if(fTmp1 < fTmp3) fTmp1 = fTmp3;
	      if(fTmp2 > fTmp4) fTmp2 = fTmp4;
	      //Tegloto e chastta na starata, kojato e vytre;
	      fTmp = (fTmp2-fTmp1)/(2*XErr1[i])*(fTmp2-fTmp1)/(2*XErr[j]);
	      alFill(pCurr->ID,(fTmp1+fTmp2)/2.f,fTmp*YArr1[i]*log(XArr1[i]/eV)/log(XArr[j]/eV),1,false);
	      if(pCurr->pBinErrors){
		pPoint1 = pVector->point(j);
		pPoint1->vScalar()->setValue(pPoint1->vScalar()->value()+
					     pCurr->pBinErrors->point(i)->vScalar()->value()*fTmp);
	      }
	      break;
	    }
	  }
	}
	bOptimizeHistos = true;
	if(pCurr->pBinErrors){
	  delete pCurr->pBinErrors;
	  pCurr->pBinErrors = pVector;
	}
	delete XArr1;
	c++;
      }
    }
  }
  return (c!=0);
}

double alGetEnergy(int nHisto,double Energy)
{
  LPHISTONODE pDummy;
  LPHISTONODE pCurr = FindHisto(nHisto,&pDummy);
  float *XArr,*YArr,*XErr,*YErr,fTmp;
  int nBins;
  if(pCurr==NULL) return 0;
  nBins = pCurr->pHisto->xAxis()->bins();
  XArr = new float[nBins*4];
  YArr = XArr+nBins;
  XErr = YArr+nBins;
  YErr = XErr+nBins;
  HREBIN(pCurr->ID,XArr,YArr,XErr,YErr,nBins,1,nBins);
  if(Energy<XArr[0]-XErr[0]){
    fTmp = XArr[0];
    delete XArr;
    return (double)fTmp;
  }
  fTmp=XArr[nBins-1];
  for(int i=0;i<nBins;i++){
    if((Energy >= XArr[i]-XErr[i])&&
       (Energy <= XArr[i]+XErr[i])){
      fTmp = XArr[i];
      break;
    }
  }
  delete XArr;
  return (double)fTmp;
}

double alGetEnergyRange(int nHisto,double Energy)
{
  LPHISTONODE pDummy;
  LPHISTONODE pCurr = FindHisto(nHisto,&pDummy);
  float *XArr,*YArr,*XErr,*YErr,fTmp;
  int nBins;
  if(pCurr==NULL) return 0;
  nBins = pCurr->pHisto->xAxis()->bins();
  XArr = new float[nBins*4];
  YArr = XArr+nBins;
  XErr = YArr+nBins;
  YErr = XErr+nBins;
  HREBIN(pCurr->ID,XArr,YArr,XErr,YErr,nBins,1,nBins);
  if(Energy<XArr[0]-XErr[0]){
    fTmp = 2*XErr[0];
    delete XArr;
    return (double)fTmp;
  }
  fTmp=2*XErr[nBins-1];
  for(int i=0;i<nBins;i++){
    if((Energy >= XArr[i]-XErr[i])&&
       (Energy <= XArr[i]+XErr[i])){
      fTmp = 2*XErr[i];
      break;
    }
  }
  delete XArr;
  return (double)fTmp;
}
