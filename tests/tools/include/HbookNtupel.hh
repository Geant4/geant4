#ifndef HbookNtupel_h
#define HbookNtupel_h

#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "HbookType.hh"
#include "g4std/vector"

extern "C"
{
  void hropen_(int &,char *,char *,char *,int &,int &,int,int,int);
  void hbnt_ (int &,const char *,const char *,int,int);
  void hbname_(int&,char*,void*,char*,int,int);
  void hbnamc_(int&,char*,void*,char*,int,int);
  void hfnt_(int &);
  void hrend_(char *,int);
  void hrout_(int &,int &,const char*,int);
  
};

template<class T>
class HbookNtupel : public HbookType
{
public:

  HbookNtupel(const char*NN,vector<const char*>names,const char*fn,
	      const char *mode="w")
  {
    char *filen;
    
    NtupelName=strdup(NN);
    for(vector<const char*>::iterator it=names.begin();it!=names.end();it++)
      MyNames.push_back(strdup(*it));
    if(fn)
      {
	filename=strdup(fn);
	filen=new char[strlen(filename)+6];
	strcpy(filen,filename);
	strcat(filen,".ntpl");
	if(strcmp(mode,"w")==0)
	  handle=open(filen,O_CREAT | O_WRONLY, 0644);
	else
	  handle=open(filen,O_RDONLY);
	if(handle<0)
	  G4cerr << "Error opening ntupel buffer:" << strerror(errno) << G4endl;
	delete[] filen;
      }
    else
      filename=0;
  }

  ~HbookNtupel()
  {
    for(vector<char*>::iterator it=MyNames.begin();it!=MyNames.end();it++)
      free(*it);
    free(filename);
    if(handle>=0) close(handle);
  }
    
  void accumulate(const T&);

  int ReadFromTupel(T&);

  void goHbook();

  void FlushBuffer();

private:

  char *NtupelName;
  vector<char *> MyNames;
  char *filename;
  int handle;
  vector<T> buffer;
};


template<class T>
void HbookNtupel<T>::accumulate(const T& x)
{
  if(filename && buffer.size()>(65536/sizeof(T))-1) // We have filled ~64k
    {
      write(handle,&(buffer[0]),sizeof(T)*buffer.size());
      buffer.erase(buffer.begin(),buffer.end());
      buffer.push_back(x);
    }
  else
    buffer.push_back(x);
}
    

template<class T>
void HbookNtupel<T>::goHbook()
{
  int lun=99,lrec=1024,ierr=0,check,FormSize=0,i;
  char mode[]="NX";
  vector<T>::iterator it;
  char *ForHbname,*filen,*hbookfile,*lunname,*ptr;
  T buf;

  if(filename)
    {
      write(handle,&(buffer[0]),sizeof(T)*buffer.size());
      buffer.erase(buffer.begin(),buffer.end());
      close(handle);
      filen=new char[strlen(filename)+6];
      hbookfile=new char[strlen(filename)+7];
      strcpy(filen,filename);
      strcpy(hbookfile,filename);
      strcat(filen,".ntpl");
      strcat(hbookfile,".hbook");
      ptr=strtok(filename,"/");
      lunname=ptr;
      while(ptr)
	{
	  lunname=ptr;
	  ptr=strtok(0,"/");
	}
      handle=open(filen,O_RDONLY);
      if(handle<0)
	G4cerr << "Error reading buffer file:" << strerror(errno) << G4endl;
      hropen_(lun,lunname,hbookfile,mode,lrec,ierr,strlen(lunname),
	      strlen(hbookfile),strlen(mode));
      hbnt_(id,NtupelName,"D",strlen(NtupelName),1);
      delete[] filen;
      delete[] hbookfile;
    }
  else
    {
      hbnt_(id,NtupelName,"M",strlen(NtupelName),1);
      lunname=filename;
    }

  // Construct the hbook ntupel form
  for(i=0;i<MyNames.size();i++)
    FormSize+=strlen(MyNames[i])+1;
  ForHbname=new char[FormSize+1];
  for(i=1,strcpy(ForHbname,MyNames[0]);i<MyNames.size();i++)
    {
      strcat(ForHbname,",");
      strcat(ForHbname,MyNames[i]);
    }

#ifndef __linux__
  hbname_(id,lunname,&buf,ForHbname,strlen(lunname),strlen(ForHbname));
#else
  if(filename)
    hbname_(id,lunname,&buf,ForHbname,strlen(lunname),strlen(ForHbname));
  else
    hbname_(id,lunname,&buf,ForHbname,0,strlen(ForHbname));
#endif

  delete[] ForHbname;

  if(filename)
    {
      for(check=read(handle,&buf,sizeof(T));check==sizeof(T);
	  check=read(handle,&buf,sizeof(T)))
	hfnt_(id);
      close(handle);
      handle=-1;
    }
	  

  for(it=buffer.begin();it!=buffer.end();it++)
    {
      buf=*it;
      hfnt_(id);
    }


  hrout_(id,ierr," ",1);
  hrend_(lunname,strlen(lunname));
}


template<class T>
int HbookNtupel<T>::ReadFromTupel(T& buf)
{
  int check;

  check=read(handle,&buf,sizeof(T));

  if(check==sizeof(T))
    return 1;
  else if(check==0)
    return 0;
  else 
    {
      G4cerr << "Error reading Ntupel:" << strerror(errno) << G4endl;
      return -1;
    }
}

template<class T>
void HbookNtupel<T>::FlushBuffer()
{
  write(handle,&(buffer[0]),sizeof(T)*buffer.size());
  buffer.erase(buffer.begin(),buffer.end());
}
  
#endif





