#ifndef HbookManager_hh
#define HbookManager_hh

class HbookType;

class HbookManager
{
public:

  HbookManager();

  ~HbookManager();

  void Register(HbookType *);
  void UnRegister(HbookType *);

  void SetFilename(const char*);

private:
  int id;
  char *hbookfile;
};

extern HbookManager theHbookManager;

#endif
