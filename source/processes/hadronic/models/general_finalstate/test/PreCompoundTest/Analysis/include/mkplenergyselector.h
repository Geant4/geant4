#ifndef mkplenergyselector_h
#define mkplenergyselector_h

using namespace std;

#include "TGFrame.h"
#include "TGLayout.h"
#include "TGShutter.h"
#include "TGButton.h"

#include "mkplexpdatamanager.h"


class mkplenergyselector : public TGTransientFrame
{
public:
  mkplenergyselector(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, 
		     mkplexpdatamanager * expdata);
  ~mkplenergyselector();
  
  void AddShutterDEItem(mkplexpdatamanager * expdata);
  void AddShutterDAItem(mkplexpdatamanager * expdata);
  void AddShutterDDItem(mkplexpdatamanager * expdata);
  void AddShutterDDAItem(mkplexpdatamanager * expdata);

  virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);

private:

  const TGWindow * mkpl_father;
  TGShutter      * mkpl_shutter;
  TGLayoutHints  * mkpl_shutter_layout;
  TList          * mkpl_trash;
};

#endif

