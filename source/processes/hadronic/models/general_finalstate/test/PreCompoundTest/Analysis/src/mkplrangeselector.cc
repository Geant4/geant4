#include "mkplrangeselector.h"

#include "mkplmessages.h"

mkplrangeselector::mkplrangeselector(const TGWindow * p, const TGWindow * main,
				     UInt_t w, UInt_t h,
				     mkplexpdatamanager * expdata, const int rid)
  : TGTransientFrame(p,main,w,h), mkpl_father(main), mkpl_range_id(rid)
{
  mkpl_shutter = new TGShutter(this, kSunkenFrame);
  mkpl_trash = new TList;

  const char * name = "Ranges";
  const int id = 9000;

  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);

  TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
  TGCompositeFrame * container = (TGCompositeFrame *)item->GetContainer();
  mkpl_trash->Add(item);

  vector<pair<double,double> > ranges = expdata->ListDDAranges(mkpl_range_id);

  unsigned int i;
  for (i=0; i < ranges.size(); i++)
    {
      double range_l = ranges[i].first;
      double range_h = ranges[i].second;
      TString rang("E = [");
      rang += range_l;
      rang += ',';
      rang += range_h;
      rang += "] MeV";
      TGButton * button = new TGTextButton(container, rang.Data(), id + 1 + i);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);
    }

  // Add close button
  TGButton * button = new TGTextButton(container, "Close",  id + 888);
  mkpl_trash->Add(button);
  container->AddFrame(button, layout);
  button->Associate(this);

  mkpl_shutter->AddItem(item);
  ranges.clear();

  mkpl_shutter_layout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  AddFrame(mkpl_shutter, mkpl_shutter_layout);
    
  MapSubwindows();
  Resize(180, 500);

  // Position relative to the parent's window
  Window_t wdum;
  int ax, ay;

  gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
				  (Int_t)(((TGFrame *) main)->GetWidth()),
				  (Int_t)(((TGFrame *) main)->GetHeight() - fHeight),
				  ax, ay, wdum);
  Move(ax,ay);

  vector<double> energies = expdata->ListDDAenergies();
  char energy[10];
  sprintf(energy, "%.1f",energies[mkpl_range_id]);
  TString tit(energy);
  tit += " MeV";
  SetWindowName(tit);

  MapWindow();
}

mkplrangeselector::~mkplrangeselector()
{
  delete mkpl_shutter;
  delete mkpl_shutter_layout;
  mkpl_trash->Delete();
  delete mkpl_trash;
}

void mkplrangeselector::CloseWindow()
{
  SendMessage(mkpl_father, mkpl_mkmsg(RANGE_SELECTOR_MSG,RANGE_SELECTOR_DELETED), 0, 0);
  DeleteWindow();
  return;
}

Bool_t mkplrangeselector::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  MKPL_MSG sm;
  Long_t typ(0);
  if (param1==9888) sm = RANGE_SELECTOR_CLOSE;
  else 
    {
      sm = ANGLE_SELECTOR_ONE;
      typ = param1 - 1000*(param1/1000) - 1;
    }
  SendMessage(mkpl_father, mkpl_mkmsg(RANGE_SELECTOR_MSG,sm), mkpl_range_id, typ);
  return kTRUE;
}
