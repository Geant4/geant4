#include "mkplenergyselector.h"

#include "mkplmessages.h"
#include "TGPicture.h"


mkplenergyselector::mkplenergyselector(const TGWindow *p, const TGWindow *main, 
				       UInt_t w, UInt_t h, 
				       mkplexpdatamanager * expdata)
  : TGTransientFrame(p, main, w, h), mkpl_father(main)
{
  mkpl_shutter = new TGShutter(this, kSunkenFrame);
  mkpl_trash = new TList;
  
  if (expdata->GetFileManager()->GetExpDEtree()->GetEntries() > 0)
    AddShutterDEItem(expdata);
  if (expdata->GetFileManager()->GetExpDAtree()->GetEntries() > 0)
    AddShutterDAItem(expdata);
  if (expdata->GetFileManager()->GetExpDDAtree()->GetEntries() > 0)
    AddShutterDDAItem(expdata);
  if (expdata->GetFileManager()->GetExpDDtree()->GetEntries() > 0)
    AddShutterDDItem(expdata);
  
  mkpl_shutter_layout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  AddFrame(mkpl_shutter, mkpl_shutter_layout);
  
  MapSubwindows();
  Resize(150, 400);
  
  // Position relative to the parent's window
  Window_t wdum;
  int ax, ay;
  //  gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(), 
  //				  (Int_t)(((TGFrame *) main)->GetWidth() - fWidth) >> 1,
  //				  (Int_t)(((TGFrame *) main)->GetHeight() - fHeight) >> 1,
  //				  ax, ay, wdum);
  gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(), 
				  (Int_t)(((TGFrame *) main)->GetWidth()),
				  (Int_t)(((TGFrame *) main)->GetHeight() - fHeight),
				  ax, ay, wdum);

  Move(ax,ay);
  
  
  SetWindowName(expdata->GetFileManager()->GetExpFile()->GetBaseName().c_str());
  
  MapWindow();
}


mkplenergyselector::~mkplenergyselector()
{
  delete mkpl_shutter;
  delete mkpl_shutter_layout;
  mkpl_trash->Delete();
  delete mkpl_trash;
}

void mkplenergyselector::CloseWindow()
{
  SendMessage(mkpl_father, mkpl_mkmsg(ENERGY_SELECTOR_MSG,ENERGY_SELECTOR_DELETED), 0, 0);
  DeleteWindow();
  return;
}

void mkplenergyselector::AddShutterDEItem(mkplexpdatamanager * expdata)
{
  const char * name = "DE"; 
  const int id = 5000;
  
  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);
  
  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  mkpl_trash->Add(item);
  TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();

  
  vector<double> energies = expdata->ListDEenergies();
  
  for (unsigned int i = 0; i < energies.size(); i++)
    {
      double E = energies[i];
      TString ts = "T = ";
      ts += E;
      ts += " MeV";
      TGButton * button = new TGTextButton(container,ts.Data(),id+1+i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }
  mkpl_shutter->AddItem(item);
  energies.clear();
  return;
}

void mkplenergyselector::AddShutterDAItem(mkplexpdatamanager * expdata)
{
  const char * name = "DA"; 
  const int id = 6000;
  
  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);
  
  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();
  mkpl_trash->Add(item);
  
  vector<double> energies = expdata->ListDAenergies();
  
  for (unsigned int i = 0; i < energies.size(); i++)
    {
      double E = energies[i];
      TString ts = "T = ";
      ts += E;
      ts += " MeV";
      TGButton * button = new TGTextButton(container,ts.Data(),id+1+i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }
  mkpl_shutter->AddItem(item);
  energies.clear();
  return;
}

void mkplenergyselector::AddShutterDDItem(mkplexpdatamanager * expdata)
{
  const char * name = "DD"; 
  const int id = 7000;

  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);
  
  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();
  mkpl_trash->Add(item);
  
  vector<double> energies = expdata->ListDDenergies();

  for (unsigned int i = 0; i < energies.size(); i++)
    {
      double E = energies[i];
      TString ts = "T = ";
      ts += E;
      ts += " MeV";
      TGButton * button = new TGTextButton(container,ts.Data(),id+1+i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }
  mkpl_shutter->AddItem(item);
  energies.clear();
  return;
}

void mkplenergyselector::AddShutterDDAItem(mkplexpdatamanager * expdata)
{
  const char * name = "DDA"; 
  const int id = 8000;

  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);
  
  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();
  mkpl_trash->Add(item);
  
  vector<double> energies = expdata->ListDDAenergies();

  for (unsigned int i = 0; i < energies.size(); i++)
    {
      double E = energies[i];
      TString ts = "T = ";
      ts += E;
      ts += " MeV";
      TGButton * button = new TGTextButton(container,ts.Data(),id+1+i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }
  mkpl_shutter->AddItem(item);
  energies.clear();
  return;
}


Bool_t mkplenergyselector::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  MKPL_MSG sm;
  Long_t typ = param1/1000;
  if (typ == 5) sm = ENERGY_SELECTOR_DE;
  else if (typ == 6) sm = ENERGY_SELECTOR_DA;
  else if (typ == 7) sm = ENERGY_SELECTOR_DD;
  else if (typ == 8) sm = ENERGY_SELECTOR_DDA;
  else return kFALSE;
  typ = param1 - typ*1000 - 1;
  SendMessage(mkpl_father, mkpl_mkmsg(ENERGY_SELECTOR_MSG,sm), typ, 0);  
  return kTRUE;
}
