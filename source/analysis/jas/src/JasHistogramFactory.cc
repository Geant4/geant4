// JasHistogramFactory.cpp: implementation of the JasHistogramFactory class.
//
//////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_BUILD_JAS

#include <jni.h>
#include <string.h>
#include <stdlib.h>

#include "G4ios.hh"

#include "JasHistogramFactory.h"
#include "JasHistogram1D.h"

static JNIEnv* env = 0;
static JavaVM* jvm = 0;       


JasHistogramFactory::JasHistogramFactory()
:fJob(NULL)
{
  createJVM();

  if(env==NULL) return;

  const char* title = "G4Job";
  // Create a Job to encapsulate the histograms
  jclass cls = env->FindClass("jas/server/HistogramServer");
  if (cls == NULL) {
    error("Could not find class jas.server.HistogramServer");
    return;
  }
  jmethodID constructor = 
    env->GetMethodID(cls,"<init>","(Ljava/lang/String;)V");
  if (constructor == NULL) {
    error("Could not find constructor");
    return;
  }
  jstring jtitle = env->NewStringUTF(title);
  fJob = env->NewObject(cls,constructor,jtitle);
  if (fJob == NULL) {
    error("Could not create job");
    return;
  }
}

JasHistogramFactory::~JasHistogramFactory()
{
  destroyJVM();
}
IHistogram1D* JasHistogramFactory::create1D(
 const std_string& name
,int,double,double
,int
)
{
  return new JasHistogram1D(this,name.c_str());
}
IHistogram2D* JasHistogramFactory::create2D(
 const std_string& name
,int,double,double
,int,double,double
,int
)
{
  return 0;
}
void JasHistogramFactory::error(const char* msg)
{
  G4cerr << "G4 JasHistogramFactory : Error: " << msg << G4std::endl;
  //exit(1);
}
JNIEnv* JasHistogramFactory::getEnv()
{
  return env;
}
void JasHistogramFactory::createJVM()
{
  JavaVMInitArgs vm_args;
  JavaVMOption options[2];
  
  char* cp = getenv("CLASSPATH");
  if(cp==NULL) {
    error("env variable CLASSPATH not setted");
    return;    
  }

  //char* cp ="c:\\program files\\java analysis studio\\lib\\hep.jar";
  char* opt = new char[strlen(cp)+100];
  strcpy(opt,"-Djava.class.path=");
  strcat(opt,cp);
  G4cout << "cp=" << opt <<G4std::endl;
  options[0].optionString = opt;
  options[1].optionString = "-verbose";
  
  vm_args.version = JNI_VERSION_1_2;
  vm_args.options = options;
  vm_args.nOptions = 1;
  vm_args.ignoreUnrecognized = 1;
  
  int rc = JNI_CreateJavaVM(&jvm, (void **)&env, &vm_args);
  delete [] opt;

  if (rc < 0)  {
    error("Failed to create Java VM");
    return;
  }
}
void JasHistogramFactory::destroy(IHistogram*){}
void JasHistogramFactory::destroyJVM()
{
  if(jvm) jvm->DestroyJavaVM();
  jvm = 0;
  env = 0;
  fJob = 0;
}

#endif
