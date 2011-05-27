//-----------------------------------------------------------------------------
//  $Header: /home/ppaprock/projects/migration/FINAL_geant4/CVSRepo/geant4/source/processes/hadronic/models/chiral_inv_phase_space/test/CHBook.cc,v 1.1 2000-08-17 14:17:14 mkossov Exp $
//
//  COOL Program Library
//  Copyright (C) CERES collaboration, 1996
//
//  Implementation of HBOOK wrapper classes.
//
//-----------------------------------------------------------------------------

#include "CHBook.h"
#include <stdlib.h>
#include <stdio.h>

long pawc_[SizeOfPawCommon];
long quest_[100];

const int MaxTagLength = 8;                   // HBOOK limitation

int CHObject::HBookInitialized = 0;

CHObject::CHObject()
{
   if (!HBookInitialized) {
      int size = SizeOfPawCommon;
      hlimit_(&size);
      HBookInitialized = 1;
   }
}

CHBook::CHBook()
{
   idd = 0;
   booked = 0;
   title = 0;
}

CHBook::~CHBook()
{
   delete [] title;
}

int CHBook::nextId(int wish)
{
   while(hexist_(&wish)) wish++;
   return wish;
}

int CHBook::entries()
{
   int n = 0;
   if (booked) hnoent_(&idd, &n);
   return n;
}

CHBookHistogram::CHBookHistogram() {}

CHBookHistogram::CHBookHistogram(const CHBookHistogram& old)
{
   idd = nextId(old.idd);
   title = new char[strlen(old.title)+1];
   strcpy(title, old.title);
   hcopy_((int*)&old.idd, &idd, title, strlen(title));
   booked = 1;
}

CHBookHistogram& CHBookHistogram::operator= (const CHBookHistogram& from)
{
   if (this == &from)                      // same class object
      return *this;
   else {
      delete [] title;                     // free old store
      hdelet_(&idd);
      strcpy(title, from.title);
      hcopy_((int*)&from.idd, &idd, title, strlen(title));
      booked = 1;
      return *this;
   }
}

CHBookHistogram::~CHBookHistogram() { hdelet_(&idd); }

float CHBookHistogram::max() { return hmax_(&idd); }

float CHBookHistogram::min() { return hmin_(&idd); }

float CHBookHistogram::sum() { return hsum_(&idd); }

void CHBookHistogram::setOpt(const char *chopt)
{
   hidopt_(&idd, (char*)chopt, strlen(chopt));
}

void CHBookHistogram::reset()
{
   hreset_(&idd, title, strlen(title));
}

void CHBookHistogram::print()
{
   hprint_(&idd);
}

CHBookHisto::CHBookHisto(int id, const char* text, int nbins, float x1, float x2)
{
   idd = nextId(id);
   _init(text, nbins, x1, x2);
}

CHBookHisto::CHBookHisto() {}

CHBookHisto::CHBookHisto(const char* text, int nbins, float x1, float x2)
{
   idd = nextId(1);
   _init(text, nbins, x1, x2);
}

void CHBookHisto::_init(const char* text, int nbins, float x1, float x2)
{
   title = new char[strlen(text)+1];
   strcpy(title, text);
   float words = 0;
   hbook1_(&idd, title, &nbins, &x1, &x2, &words, strlen(title));
   setOpt("STAT");
   booked = 1;
}

float CHBookHisto::mean()
{
   int icase = 1;
   int num = 0;
   char choice = 0;
   return hstati_(&idd, &icase, &choice, &num, 0);
}

float CHBookHisto::sigma()
{
   int icase = 2;
   int num = 0;
   char choice = 0;
   return hstati_(&idd, &icase, &choice, &num, 0);
}

CHBookHisto2::CHBookHisto2(int id, const char* text, int nxbins, float x1, float x2,
                                                     int nybins, float y1, float y2)
{
   idd = nextId(id);
   _init(text, nxbins, x1, x2, nybins, y1, y2);
}

CHBookHisto2::CHBookHisto2() {}

CHBookHisto2::CHBookHisto2(const char* text, int nxbins, float x1, float x2,
                                             int nybins, float y1, float y2)
{
   idd = nextId(1);
   _init(text, nxbins, x1, x2, nybins, y1, y2);
}

void CHBookHisto2::_init(const char* text, int nxbins, float x1, float x2,
                                           int nybins, float y1, float y2)
{
   title = new char[strlen(text)+1];
   strcpy(title, text);
   float words = 0;
   hbook2_(&idd, title, &nxbins, &x1, &x2, &nybins, &y1, &y2, &words, strlen(title));
   setOpt("STAT");
   booked = 1;
}

CHBookFile::CHBookFile() {}

CHBookFile::~CHBookFile()
{
}

CHBookFile::CHBookFile(const CHBookFile&) {}

CHBookFile::CHBookFile(const char* file, int lun, int forWriting)
{
  setFile(file, lun, forWriting);
}

void CHBookFile::setFile(const char* file, int lun, int forWriting)
{
   quest_[9] = ValueForIquest10;
   filename = new char[strlen(file)+1];
   strcpy(filename, file);
   context = new char[strlen("HISTOS")+1];
   strcpy(context, "HISTOS");
   mode = new char[2];
   if (forWriting)
     strcpy(mode, "N");
   else
     strcpy(mode, " ");

   int recordSize = 1024;
   hropen_(&lun, context, filename, mode, &recordSize, &rc,
           strlen(context), strlen(filename), strlen(mode));

}

CHBookFile& CHBookFile::operator=(const CHBookFile&)
{
   return *this;
}

void CHBookFile::close() {
  int icycle = 0;
  context = new char[strlen("//HISTOS")+1];
  strcpy(context, "//HISTOS");
  char blank[4];
  strcpy(blank, " ");
  hrendc_(context, strlen(context));
  delete [] filename;
  delete [] context;
  delete [] mode;
}

void CHBookFile::saveAndClose() {
  int icycle = 0;
  context = new char[strlen("//HISTOS")+1];
  strcpy(context, "//HISTOS");
  char blank[4];
  strcpy(blank, " ");
  hcdir_(context, blank, strlen(context), strlen(blank));
  int selectID = 0;
  hrout_(&selectID, &icycle, blank, strlen(blank));
  hrend_(context, strlen(context));
  hcdir_(context, blank, strlen(context), strlen(blank));    // needed
  delete [] filename;
  delete [] context;
  delete [] mode;
}

CHBookTuple::CHBookTuple() {}

CHBookTuple::CHBookTuple(const CHBookTuple&) {}

CHBookTuple::CHBookTuple(const char* name, int n) : ntags(n)
{
   idd = nextId(1);
   _init(name);
}

CHBookTuple::CHBookTuple(int id) : ntags(1)
{
   idd = id;
   _init("unnamed");
}

CHBookTuple::CHBookTuple(int id, const char* name, int n) : ntags(n)
{
   idd = nextId(id);
   _init(name);
}

void CHBookTuple::_init(const char* name)
{
   title = new char[strlen(name)+1];
   strcpy(title, name);
   itag = 0;
   tags = new char[MaxTagLength*ntags+ntags];
   for (int i=0; i<ntags; i++)
      sprintf(&tags[i*MaxTagLength],"%-8.8s", "unknown");
}

CHBookTuple::~CHBookTuple() { delete [] tags; }

CHBookTuple& CHBookTuple::operator=(const CHBookTuple&)
{
  return *this;
}

CHBookTuple& CHBookTuple::setTag(const char *tag)
{
   if (booked) {
      cerr << "CHBookTuple::setTag():" << endl;
      cerr << "\tWARNING" << endl;
      cerr << "\tTuple already booked. No further tags can be added." << endl;
   }
   else {
      if (itag < ntags) {
         sprintf(&tags[itag*MaxTagLength],"%-8.8s", tag);
         itag++;
      }
   }
   return *this;
}

CHBookTuple& CHBookTuple::operator<< (const char *tag)
{
   return this->setTag(tag);
}

void CHBookTuple::book()
{
   if (!booked) {
      char *context = new char[strlen("HISTOS")+1];
      strcpy(context, "HISTOS");
      int nPrime = 4096;
      hbookn_(&idd, title, &ntags, context, &nPrime, tags,
	      strlen(title), strlen(context), MaxTagLength);
      booked = 1;
      delete [] context;
   }
}

void book(CHBookTuple& tuple) { tuple.book(); }

void CHBookTuple::operator<< (void (*pf)(CHBookTuple& tuple))
{
   (*pf)(*this);
}
