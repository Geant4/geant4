#include "HbookManager.hh"
#include "HbookType.hh"
#include <string.h>

extern "C" 
{
  struct {int h[1000000];} pawc_;
  void hlimit_(int&);
  void hcdir_(const char*,const char*,int,int);
  void hrput_(int&,const char*,const char*,int,int);
  void hdelet_(int&);
}

HbookManager::HbookManager()
{
  int temp=1000000;
  hlimit_(temp);
  id=0;
  hbookfile=0;
}
  

HbookManager::~HbookManager()
{
  int temp=0;
  char pbug[10];


  strcpy(pbug,"N");
  hcdir_("//PAWC"," ",strlen("//PAWC"),strlen(" "));
  hrput_(temp,hbookfile,pbug,
	 strlen(hbookfile),strlen(pbug));

  delete[] hbookfile;

}


void HbookManager::Register(HbookType *hb)
{
  hb->SetID(++id);
}
  
void HbookManager::UnRegister(HbookType *hb)
{
  int temp=hb->GetHbookID();
  hdelet_(temp);
}

void HbookManager::SetFilename(const char*fn)
{
  if(hbookfile) delete[] hbookfile;
  hbookfile=new char[strlen(fn)+1];
  strcpy(hbookfile,fn);
}


HbookManager theHbookManager;

