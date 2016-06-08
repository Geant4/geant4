// JHistogramFactory.h: interface for the JHistogramFactory class.
//////////////////////////////////////////////////////////////////////
#ifndef JASHISTOGRAMFACTORY
#define JASHISTOGRAMFACTORY

#if defined(G4ANALYSIS_BUILD_JAS)

#include <jni.h>

#include <IHistogramFactory.h>

class JasHistogramFactory  
: public IHistogramFactory
{
 public: // IMethods :
  virtual IHistogram1D* create1D(const std_string&,
				 int=0,double=0,double=0,
				 int=0);
  virtual IHistogram2D* create2D(const std_string&,
				 int=0,double=0,double=0,
				 int=0,double=0,double=0,
				 int=0);
  virtual void destroy(IHistogram*);
 public:
  JasHistogramFactory();
  virtual ~JasHistogramFactory();
  JNIEnv* getEnv();
  void error(const char* message);
 private:
  jobject fJob;
  void createJVM();
  void destroyJVM();
};

#endif

#endif
