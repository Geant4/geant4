//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

