#ifndef HbookType_hh
#define HbookType_hh

#include "HbookManager.hh"

class HbookType
{
public:

  HbookType()
  {
    theHbookManager.Register(this);
  }

  virtual ~HbookType()
  {
    theHbookManager.UnRegister(this);
  }    

  // To define this const isn´t fair but avoids 
  // anoying compiler warnings
  int GetHbookID() const
  {
    return id;
  }

protected:

  int id;

private:

  void SetID(int v)
  {
    id=v;
  }

  friend class HbookManager;
};

#endif
