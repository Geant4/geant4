#ifndef mkplcomparisonselector_h
#define mkplcomparisonselector_h

using namespace std;

#include "TGFrame.h"
#include "TGLayout.h"
#include "TGShutter.h"
#include "TGButton.h"

#include "mkplsimmanager.h"

class mkplcomparisonselector : public TGTransientFrame
{
public:
  mkplcomparisonselector(const TGWindow *p, const TGWindow * main, UInt_t w, UInt_t h, 
			 mkplsimmanager * simdata);
  ~mkplcomparisonselector();

  

  virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);

private:
  bool AddShutterTestItem(mkplsimmanager * simdata);
  bool AddShutterCompItem(mkplsimmanager * simdata);

private:
  const TGWindow   * mkpl_father;
  TGShutter        * mkpl_shutter;
  TGLayoutHints    * mkpl_shutter_layout;
  TList            * mkpl_trash;
};

#endif
