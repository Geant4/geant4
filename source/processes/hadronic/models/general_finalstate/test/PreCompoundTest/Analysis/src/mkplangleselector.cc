#include "mkplangleselector.h"

#include "mkplmessages.h"

mkplangleselector::mkplangleselector(const TGWindow *p, const TGWindow * main,
				     UInt_t w, UInt_t h, 
				     mkplexpdatamanager * expdata, const int eid)
  : TGTransientFrame(p,main,w,h), mkpl_father(main), mkpl_energy_id(eid)
{
  mkpl_shutter = new TGShutter(this, kSunkenFrame);
  mkpl_trash = new TList;

  const char * name = "Angles";
  const int id = 8000;

  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);

  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();
  mkpl_trash->Add(item);

  vector<double> angles = expdata->ListDDangles(mkpl_energy_id);

  unsigned int i;
  for (i = 0; i < angles.size(); i++)
    {
      double angle = angles[i];
      TString ang("Theta = ");
      ang += angle;
      ang += " degrees";
      TGButton * button = new TGTextButton(container, ang.Data(), id + 1 + i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }
  // Add button for all angles
  TGButton * button = new TGTextButton(container, "All", id + 999);
  mkpl_trash->Add(button);
  container->AddFrame(button, layout);
  button->Associate(this);

  // Add close button
  button = new TGTextButton(container, "Close",  id + 888);
  mkpl_trash->Add(button);
  container->AddFrame(button, layout);
  button->Associate(this);

  mkpl_shutter->AddItem(item);
  angles.clear();

  mkpl_shutter_layout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  AddFrame(mkpl_shutter, mkpl_shutter_layout);

  MapSubwindows();
  Resize(180, 500);

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
  
  vector<double> energies = expdata->ListDDenergies();
  char energy[10];
  sprintf(energy,"%.1f",energies[mkpl_energy_id]);
  TString tit(energy);
  tit += " MeV";
  SetWindowName(tit);
 
  MapWindow();
}


mkplangleselector::~mkplangleselector()
{
  delete mkpl_shutter;
  delete mkpl_shutter_layout;
  mkpl_trash->Delete();
  delete mkpl_trash;
}

void mkplangleselector::CloseWindow()
{
  SendMessage(mkpl_father, mkpl_mkmsg(ANGLE_SELECTOR_MSG,ANGLE_SELECTOR_DELETED), 0, 0);
  DeleteWindow();
  return;
}

Bool_t mkplangleselector::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  MKPL_MSG sm;
  Long_t typ(0);
  if (param1==8888) sm = ANGLE_SELECTOR_CLOSE;
  else if (param1==8999) sm = ANGLE_SELECTOR_ALL;
  else 
    {
      sm = ANGLE_SELECTOR_ONE;
      typ = param1 - 1000*(param1/1000) - 1;
    }
  SendMessage(mkpl_father, mkpl_mkmsg(ANGLE_SELECTOR_MSG,sm), mkpl_energy_id, typ);
  return kTRUE;
}
