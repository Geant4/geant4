#ifndef mkplangleselector_h
#define mkplangleselector_h

using namespace std;

#include "TGFrame.h"
#include "TGLayout.h"
#include "TGShutter.h"
#include "TGButton.h"

#include "mkplexpdatamanager.h"

class mkplangleselector : public TGTransientFrame
{
public:
  mkplangleselector(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
		    mkplexpdatamanager * expdata, const int eid);
  ~mkplangleselector();

  virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);

private:
  const TGWindow * mkpl_father;
  const int        mkpl_energy_id;
  TGShutter      * mkpl_shutter;
  TGLayoutHints  * mkpl_shutter_layout;
  TList          * mkpl_trash;
};

#endif
