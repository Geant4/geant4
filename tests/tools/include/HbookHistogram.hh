#ifndef HbookHistogram_h
#define HbookHistogram_h

#include <string.h>
#include "HbookType.hh"

extern "C"
{
  void hbook1_(const int &,const char *,const int &,const float &,const float &,const float &,int);
  void hbprof_(const int &,const char *,const int &,const float &,const float &,const float &,const float&,const char *,int,int);
  void hbook2_(const int &,const char *,const int &,const float &,const float &,const int &,const float &,const float &,const float &,int);
  void hfill_(const int &,const float &,const float &,const float &);
  void hdelet_(int &);
  int  hexistpp_(int &);
  void hcopy_(int &,int &,char *,int);
  void hreset_(int &,char *,int);
  void hbarx_(int &);
  void hbary_(int &);
  void hbookb_(int &,const char *,int &,float *,float &,int);
  float hsum_(const int &);
  void hopera_(int&,const char*,int&,int&,float &,float &,int);
  void hidopt_(int&,const char*,int);
};

class HbookHistogram : public HbookType
{
public:
  HbookHistogram(){}

  HbookHistogram(const char *n,int b,float l,float u,float s=0.)
  {
    hbook1_(id,n,b,l,u,s,strlen(n));
    hidopt_(id,"STAT",4);
  }

  HbookHistogram(const char *n,int bx,float lx,float ux,
		 int by,float ly,float uy,float s=0.)
  {
    hbook2_(id,n,bx,lx,ux,by,ly,uy,s,strlen(n));
    hidopt_(id,"STAT",4);
  }


  HbookHistogram(const char *n,int channel,float lx,float ux,
		 float ly,float uy,char *opt=" ")
  {
    hbprof_(id,n,channel,lx,ux,ly,uy,opt,strlen(n),strlen(opt));
  }

  HbookHistogram(const char *n,int b,double *bounds,float s=0.)
  {
    float *fbounds;
    fbounds=new float[b+1];
    for(int i=0;i<b+1;i++) fbounds[i]=bounds[i];
    hbookb_(id,n,b,fbounds,s,strlen(n));
    delete[] fbounds;
    hidopt_(id,"STAT",4);
    //hbarx_(GetHbookID);
  }

  HbookHistogram(const HbookHistogram &cp)
  {
    int i=cp.GetHbookID();
    if(hexistpp_(id)) hdelet_(id);
    hcopy_(i,id," ",1);
  }


  ~HbookHistogram()
  {
  } 


  void accumulate(float x,float y=0.,float weight=1.)
  {
    hfill_(id,x,y,weight);
  }

  void reset(){hreset_(id," ",1);}

  double sum()const{return hsum_(id);}

  //Operation

  HbookHistogram & operator +=(const HbookHistogram &b)
  {
    int idcp=b.GetHbookID();
    float one=1.0;
    hopera_(id,"+",idcp,id,one,one,1);
    return *this;
  }

  HbookHistogram & operator -=(const HbookHistogram &b)
  {
    int idcp=b.GetHbookID();
    float one=1.0;
    hopera_(id,"-",idcp,id,one,one,1);
    return *this;
  }

  HbookHistogram & operator *=(const HbookHistogram &b)
  {
    int idcp=b.GetHbookID();
    float one=1.0;
    hopera_(id,"*",idcp,id,one,one,1);
    return *this;
  }

  HbookHistogram & operator /=(const HbookHistogram &b)
  {
    int idcp=b.GetHbookID();
    float one=1.0;
    hopera_(id,"/",idcp,id,one,one,1);
    return *this;
  }

};

#endif
