//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

