#include "mkplcomparisonselector.h"

#include "mkplmessages.h"
#include "mkplcomparisonhistograms.h"


mkplcomparisonselector::mkplcomparisonselector(const TGWindow * p, const TGWindow * main,
					       UInt_t w, UInt_t h,
					       mkplsimmanager * simdata)
  : TGTransientFrame(p, main, w, h), mkpl_father(main)
{
  mkpl_shutter = new TGShutter(this, kSunkenFrame);
  mkpl_trash = new TList;

  bool test = this->AddShutterTestItem(simdata);
  bool comp = this->AddShutterCompItem(simdata);
  
  if (!test && !comp)
    {
      TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
						 5, 5, 5, 0);
      mkpl_trash->Add(layout);

      TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString("NO PLOTS"), 0);
      mkpl_trash->Add(item);

      mkpl_shutter->AddItem(item);
    }

  mkpl_shutter_layout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  AddFrame(mkpl_shutter, mkpl_shutter_layout);

  MapSubwindows();
  Resize(150, 400);

  // Position relative to the parent's window
  Window_t wdum;
  int ax, ay;
  gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
				  (Int_t)(((TGFrame *) main)->GetWidth()),
				  (Int_t)(((TGFrame *) main)->GetHeight() - fHeight),
				  ax, ay, wdum);

  Move(ax,ay);

  SetWindowName("Plot selector");

  MapWindow();
}

mkplcomparisonselector::~mkplcomparisonselector()
{
  delete mkpl_shutter;
  delete mkpl_shutter_layout;
  mkpl_trash->Delete();
  delete mkpl_trash;
}

void mkplcomparisonselector::CloseWindow()
{
  SendMessage(mkpl_father, mkpl_mkmsg(COMPARISON_SELECTOR_MSG, COMPARISON_SELECTOR_DELETED), 0, 0);
  DeleteWindow();
  return;
}

bool mkplcomparisonselector::AddShutterTestItem(mkplsimmanager * simdata)
{
  // Get the test histograms
  mkpltesthistograms * tests = simdata->GetTestHistograms();
  bool result(false);
  if (tests) 
    {
      const int id = 9000;

      TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
						 5, 5, 5, 0);
      mkpl_trash->Add(layout);

      TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString("Tests"), id);
      mkpl_trash->Add(item);
      TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();
     
      TGButton * button = new TGTextButton(container,"Baryon num.", id+1);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Charge cons.", id+2);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Energy cons.", id+3);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Momentum cons.", id+4);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Px cons.", id+5);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Py cons.", id+6);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Pz cons.", id+7);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Theta dist.", id+8);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Theta Preeq dist.", id+9);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Theta Evap dist.", id+10);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Phi dist.", id+11);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Phi Preeq dist.", id+12);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Phi Evap dist.", id+13);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Phi Nucleon dist.", id+14);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      button = new TGTextButton(container, "Type of frag.", id+15);
      mkpl_trash->Add(button);
      container->AddFrame(button, layout);
      button->Associate(this);

      mkpl_shutter->AddItem(item);
      result = true;
    }
  return result;
}

bool mkplcomparisonselector::AddShutterCompItem(mkplsimmanager * simdata)
{
  // Get the vector of comparison histograms
  vector<mkplcomparisonhistograms*> * comparisons = simdata->GetComparisonHistograms();
  int id = 100000;
  bool result(false);

  TGLayoutHints * layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
					     5, 5, 5, 0);
  mkpl_trash->Add(layout);

  for (vector<mkplcomparisonhistograms*>::iterator it = comparisons->begin();
       it != comparisons->end(); it++)
    {
      TString name = (*it)->GetName();
      
      TGShutterItem * item = new TGShutterItem(mkpl_shutter, new TGHotString(name), id);
      mkpl_trash->Add(item);
      TGCompositeFrame * container = (TGCompositeFrame *) item->GetContainer();

      TGButton * button;

      if ((*it)->ThereIsDE())
	{
	  button = new TGTextButton(container,"DE",id+1);
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}
      
      if ((*it)->ThereIsDEcm())
	{
	  button = new TGTextButton(container, "DE (cm)", id+2);
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}
      if ((*it)->ThereIsDA())
	{
	  button = new TGTextButton(container,"DA",id+3);
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}
      
      if ((*it)->ThereIsDAcm())
	{
	  button = new TGTextButton(container, "DA (cm)", id+4);
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}

      if ((*it)->ThereIsDD())
	{
	  TString label;
	  for (int a = 0; a < (*it)->GetNumOfAngles(); a++)
	    {
	      label = "DD ";
	      char angle[5];
	      sprintf(angle,"%.1f",(*it)->GetAngle(a)*180.0/TMath::Pi());
	      label += angle;
	      label +=  "deg";
	      button = new TGTextButton(container, label, id+101+a); 
	      mkpl_trash->Add(button);
	      container->AddFrame(button, layout);
	      button->Associate(this);
	    }
	  // Add ALL button
	  button = new TGTextButton(container, "ALL angles", id+101+(*it)->GetNumOfAngles());
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}

      if ((*it)->ThereIsDDcm())
	{
	  TString label;
	  for (int a = 0; a < (*it)->GetNumOfAnglesCM(); a++)
	    {
	      label = "DD ";
	      char angle[5];
	      sprintf(angle,"%.1f",(*it)->GetAngleCM(a)*180.0/TMath::Pi());
	      label += angle;
	      label +=  "deg (cm)";
	      button = new TGTextButton(container, label, id+201+a); 
	      mkpl_trash->Add(button);
	      container->AddFrame(button, layout);
	      button->Associate(this);
	    }
	  // Add ALL button
	  button = new TGTextButton(container, "ALL angles (CM)", id+201+(*it)->GetNumOfAnglesCM());
	  mkpl_trash->Add(button);
	  container->AddFrame(button, layout);
	  button->Associate(this);
	}

      if ((*it)->ThereIsDDA())
	{
	  TString label;
	  for (int r = 0; r < (*it)->GetNumOfRanges(); r++)
	    {
	      label = "DDA [";
	      char energy[10];
	      pair<double,double> range = (*it)->GetRange(r);
	      sprintf(energy,"%.1f",range.first);
	      label += energy;
	      label += ',';
	      sprintf(energy,"%.1f",range.second);
	      label += energy;
	      label +=  "] MeV";
	      button = new TGTextButton(container, label, id+301+r); 
	      mkpl_trash->Add(button);
	      container->AddFrame(button, layout);
	      button->Associate(this);
	    }
	}

      if ((*it)->ThereIsDDAcm())
	{
	  TString label;
	  for (int r = 0; r < (*it)->GetNumOfRangesCM(); r++)
	    {
	      label = "DDA [";
	      char energy[10];
	      pair<double,double> range = (*it)->GetRangeCM(r);
	      sprintf(energy,"%.1f",range.first);
	      label += energy;
	      label += ',';
	      sprintf(energy,"%.1f",range.second);
	      label += energy;
	      label +=  "] MeV (cm)";
	      button = new TGTextButton(container, label, id+401+r); 
	      mkpl_trash->Add(button);
	      container->AddFrame(button, layout);
	      button->Associate(this);
	    }
	}
  
      
      mkpl_shutter->AddItem(item);
      id += 1000;
      result = true;
    }
  
  return result;
}


Bool_t mkplcomparisonselector::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  MKPL_MSG sm;
  Long_t a(0);
  Long_t b(0);
  if (param1 < 100000)
    {
      sm = COMPARISON_SELECTOR_TEST;
      a = param1 - 9000;
    }
  else
    {
      sm = COMPARISON_SELECTOR_COMP;
      // Get the ejectile number
      a = param1/1000 - 100;
      // Get the code
      b = param1 - 1000*(a+100);
    }
  
  SendMessage(mkpl_father, mkpl_mkmsg(COMPARISON_SELECTOR_MSG,sm), a, b);
  return kTRUE;
}
