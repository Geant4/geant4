#ifndef HEPVis_SoCounterAction_h
#define HEPVis_SoCounterAction_h 

// Inheritance :
#include <Inventor/actions/SoAction.h>

#include <Inventor/actions/SoSubAction.h>

class SoCounterAction : public SoAction {
  SO_ACTION_HEADER(SoCounterAction);
public:
  SoCounterAction();
  virtual ~SoCounterAction();  
  static void initClass(void);
  enum LookFor { NODE = 1,TYPE = 2, NAME = 3 };
  void setLookFor(LookFor);
  void setType(const SoType,SbBool = TRUE);
  void setName(const SbName);
public:
  int getCount() const;
protected:
  virtual void beginTraversal(SoNode*);
private:
  static void actionMethod(SoAction*,SoNode*);
private:
  int fCount;
  int fLookFor;
  SbName fName;
  SoType fType;
  SbBool fCheckDerived;
};

#endif