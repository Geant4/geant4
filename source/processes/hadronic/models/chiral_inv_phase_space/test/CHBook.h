//-----------------------------------------------------------------------------
//  $Header: /home/ppaprock/projects/migration/FINAL_geant4/CVSRepo/geant4/source/processes/hadronic/models/chiral_inv_phase_space/test/CHBook.h,v 1.1 2000-08-17 14:17:14 mkossov Exp $
//
//  COOL Program Library
//  Copyright (C) CERES collaboration, 1996
//
//  HBOOK wrapper classes.
//
//  Note that all initialization of HBOOK and ZEBRA is already
//  performed within the classes. No additional compile flags
//  are needed. Link packlib only. If the /pawc/ common is too
//  small to cope your requirements increase the value of the 
//  constant 'SizeOfPawCommon' to whatever you need and recompile.
//
//  Member function of the various objects:
//
//  CHBookHisto (1-dim histogram):
//        CHBookHisto(const char* name, int nbins, float x1, float x2)
//        CHBookHisto(int id, const char* name, int nbins, float x1, float x2)
//        CHBookHisto(const CHBookHisto&)
//        CHBookHisto& operator= (const CHBookHisto&)
//        void  CHBookHisto::fill(float x, float weight = 1)
//        void  CHBookHisto::fastFill(float x, float weight = 1)
//        int   CHBookHisto::id()
//        int   CHBookHisto::entries()
//        float CHBookHisto::max()
//        float CHBookHisto::min()
//        float CHBookHisto::sum()
//        void  CHBookHisto::setOpt(const char*)
//        void  CHBookHisto::print()
//        void  CHBookHisto::reset()
//        float CHBookHisto::mean()
//        float CHBookHisto::sigma()
//
//  CHBookHisto2 (2-dim histogram):
//        CHBookHisto2(const char* name, int nxbins, float x1, float x2,
//                                       int nybins, float y1, float y2)
//        CHBookHisto2(int id, const char* name, int nxbins, float x1, float x2
//                                               int nybins, float y1, float y2)
//        CHBookHisto2(const CHBookHisto2&)
//        CHBookHisto2& operator= (const CHBookHisto2&)
//        void  CHBookHisto2::fill(float x, float weight = 1)
//        void  CHBookHisto2::fastFill(float x, float weight = 1)
//        int   CHBookHisto2::id()
//        int   CHBookHisto2::entries()
//        float CHBookHisto2::max()
//        float CHBookHisto2::min()
//        float CHBookHisto2::sum()
//        void  CHBookHisto2::setOpt(const char*)
//        void  CHBookHisto2::print()
//        void  CHBookHisto2::reset()
//
//  CHBookTuple (RWN tuples):
//        CHBookTuple(int id, const char* name, int ntags)
//        CHBookTuple(const char* name, int ntags)
//        CHBookTuple& CHBookTuple::setTag(const char* tag)
//        void  CHBookTuple::book()
//        int   CHBookTuple::id()
//        int   CHBookTuple::length()
//        int   CHBookTuple::entries()
//        void  CHBookTuple::fill(float *vec)
//        CBoolean CHBookTuple::getEvent(int eventNumber, float *vec)
//
//  CHBookFile (ZEBRA RZ file):
//        CHBookFile(const char* filename, int lun = 10)
//        int isGood() const
//        void CHBookFile::saveAndClose()
//
//
//  Author: Thomas Ullrich
//  Last update: 13/03/96 initial version
//               25/03/96 added 2-dim histos, no cfortran needed any more
//                        many new member functions added, copy constructors
//                        and assignment operators defined.
//-----------------------------------------------------------------------------
#ifndef CHBOOK_H
#define CHBOOK_H

#include <iostream.h>
#include <string.h>
#include <stdlib.h>

const int SizeOfPawCommon  = 1000000;
const int ValueForIquest10 = 30000000; // max number of records in an hbook file
                                       // (default recl 1024)

extern "C" {
void  hlimit_(int*);
void  hbook1_(int*, char*, int*, float*, float*, float*, int);
void  hbook2_(int*, char*, int*, float*, float*, int*, float*, float*, float*, int);
void  hfill_(int*, float*, float*, float*);
void  hf1_(int*, float*, float*);
void  hf2_(int*, float*, float*, float*);
void  hdelet_(int*);
void  hropen_(int*, char*, char*, char*, int*, int*, int, int, int);
void  hcdir_(char*, char*, int, int);
void  hrout_(int*, int*, char*, int);
void  hrin_(int*, int*, int*);
void  hrend_(char*, int);
void  hrendc_(char*, int);
void  hbookn_(int*, char*, int*, char*, int*, char*, int, int, int);
void  hfn_(int*, float*);
void  hgnpar_(int*, char*, int);
void  hgnf_(int*, int*, float*, int*);
void  hgn_(int*, int*, int*, float*, int*);
void  hnoent_(int*, int*);
void  hprint_(int*);
void  hreset_(int*, char*, int);
void  hcopy_(int*, int*, char*, int);
void  hindex_();
int   hexist_(int*);
float hmax_(int*);
float hmin_(int*);
float hsum_(int*);
float hstati_(int*, int*, char*, int*, int);
float hidopt_(int*, char*, int);
}

//
//   'Virtual' base class.
//   Although instances can be created they are of no use.
//   The only purpose is to ensure that HBOOK is initialized
//   whenever a derived class is used.
//

class CHObject {
public:
   CHObject();
private:
   static int HBookInitialized;
};

//
//   Base class common to all histograms and tuple.
//
class CHBook : public CHObject {
public:
   CHBook();
   virtual ~CHBook();
   int id() const { return idd; }
   int entries();

protected:
   int isBooked() const { return booked; }
   int nextId(int);

protected:
   int idd;
   int booked;
   char *title;
};

//
//  Base class common to all histograms (1-2 dim)
//
class CHBookHistogram : public CHBook {
public:
   CHBookHistogram();
   CHBookHistogram(const CHBookHistogram&);
   CHBookHistogram& operator= (const CHBookHistogram&);
   virtual ~CHBookHistogram();

public:
   float max();
   float min();
   float sum();
   void setOpt(const char*);
   void print();               // ascii print out
   void reset();               // reset contents
};


//
//  One-dimensional histogram class
//
//  Create/book and fill HBOOK histograms.
//  The usual id may be defined by the user by passing the id
//  as first argument. However, if the specified id already exists
//  the next higher free id is used. The object will find its own
//  id if non is specified.
//  The destructor ~CHBookHisto() deletes the histogram from the
//  (ZEBRA) memory by invoking (HDELET).
//  See class CHBookFile for how to store histograms on file.
//
//  Example:
//
//  CHBookHisto  ptHisto("pt-dist", 100, 0., 2.);
//  ....
//  ptHisto.fill(pt);
//  cout << ptHisto.mean() << endl;
//  cout << ptHisto.entries() << endl;
//
//  Note that fill() uses the HBOOK option 'STAT', in order to
//  calculate the mean and sigma from the passed values directly.
//  This feature is ignored by fastFill(). In this case mean and
//  sigma are derived from the binned data.
//  (See also HF1, HFILL in the HBOOK manual).
//
class CHBookHisto : public CHBookHistogram {
public:
   CHBookHisto(const char* name, int nbins, float x1, float x2);
   CHBookHisto(int id, const char* name, int nbins, float x1, float x2);
   inline void fill(float x, float weight = 1);
   inline void fastFill(float x, float weight = 1);
   float mean();
   float sigma();

private:
   CHBookHisto();
   void _init(const char*, int, float, float);
};

//
//  Two-dimensional histogram class
//
//  Example:
//
//  CHBookHisto2  xyHisto("xy plane", 100, 0., 2., 200, 0., 3.);
//  ....
//  xyHisto.fill(x, y);
//  cout << xyHisto.sum() << endl;
//
//  See CHBookHisto for more info.
//
class CHBookHisto2 : public CHBookHistogram {
public:
   CHBookHisto2(const char* name, int nxbins, float x1, float x2,
                                         int nybins, float y1, float y2);
   CHBookHisto2(int id, const char* name, int nxbins, float x1, float x2,
                                                 int nybins, float y1, float y2);
   inline void fill(float x, float y, float weight = 1);
   inline void fastFill(float x, float y, float weight = 1);

private:
   CHBookHisto2();
   void _init(const char*, int, float, float, int, float, float);
};

//
//  Create and fill RWN tuples.
//
//  A CHBookFile must be instantiated in advance in order
//  to store the tuple on file. (See class CHBookFile below).
//  The usual id may be defined by the user by passing the id
//  as first argument. However, if the specified id already exists
//  the next higher free id is used. The object will find its own
//  id if non is specified.
//  The tuple is defined in 3 steps.
//  - the constructor only defines the # of tags, the title and
//    (optional) the id
//  - the individual tags are defined via the setTag() method
//  - method book() finally books the referring tuple.
//
//  Setting tags:
//  The order is important. Tags not defined by setTag() are stored
//  as "unknown". Max tag lenght is 8 character.
//
//  Example:
//
//  CHBookTuple mytuple("foo",4);
//  mytuple.setTag("x").setTag("y").setTag("z").book()
//  defines at tuple with the following tags
//  1    x
//  2    y
//  3    z
//  4    unknown
//  float vec[4];
//  ...
//  mytuple.fill(vec);                // fill the referring tuple
//
//  Instead of using the setTag() and book() member one may also use
//  the 'put to' operator << in the following way:
//  mytuple << "x" << "y" << "z" << book;
//
class CHBookTuple : public CHBook {
public:
   CHBookTuple(int id); //dd used only for reading
   CHBookTuple(int id, const char* name, int ntags);
   CHBookTuple(const char* name, int ntags);
   ~CHBookTuple();

public:
   CHBookTuple& setTag(const char* tag);
   void book();
   void remove() { hdelet_(&idd); }
   inline void fill(float *vec);
   inline int getEvent(int, float*);
   CHBookTuple& operator<< (const char *tag);
   void operator<< (void (*)(CHBookTuple&));
   int length() const { return ntags; }
   
private:
   void _init(const char*);
   CHBookTuple();
   CHBookTuple(const CHBookTuple&);
   CHBookTuple& operator= (const CHBookTuple&);

private:
   int ntags;
   char *tags;
   int itag;
   int nevents;
};
void book(CHBookTuple&);

//
//  Defines an HBOOK RZ file
//  After creating one instance all ntuples are automatically 
//  written to file. Histograms (and the residual tuples still 
//  in memory are dumped to file by invoking the saveAndClose()
//  method. Note that a CHBookFile must be defined prior to 
//  objects of type CHBookTuple.
//
class CHBookFile : public CHObject {
public:
   CHBookFile();
   CHBookFile(const char* filename, int lun = 10, int forWriting=1);
   ~CHBookFile();
   void setFile(const char* filename, int lun = 10, int forWriting=1);
   void close();
   void saveAndClose();
   int isGood() const { return rc == 0; }

private:
   CHBookFile(const CHBookFile&);
   CHBookFile& operator= (const CHBookFile&);

private:
   int rc;
   char *mode;
   char *context;
   char *filename;
   int lun;
};


//
//  Definition of inline member functions
//

void CHBookHisto::fill(float x, float weight)
{
   float dummy = 0;
   hfill_(&idd, &x, &dummy, &weight);
}

void CHBookHisto::fastFill(float x, float weight)
{
   hf1_(&idd, &x, &weight);
}

void CHBookHisto2::fill(float x, float y, float weight)
{
   hfill_(&idd, &x, &y, &weight);
}

void CHBookHisto2::fastFill(float x, float y, float weight)
{
   hf2_(&idd, &x, &y, &weight);
}

void CHBookTuple::fill(float *tuple) { hfn_(&idd, tuple); }



int CHBookTuple::getEvent(int eventNumber, float *vector)
{
  int ierror = 0;
//   static int id=0;
   if (eventNumber == 0)
   {
     int highCycle=999999;
     int zero=0;
     hrin_(&idd,&highCycle,&zero);

//     hgn_(&idd,&id,&eventNumber, vector, &ierror);

     hnoent_(&idd,&nevents);
     char *errText = "CHBookTuple::getEvent";
     hgnpar_(&idd, errText, strlen(errText));
   }

//     hgn_(&idd,&id,&eventNumber,vector, &ierror);
   if (eventNumber >= 0 && eventNumber < nevents)
   {
     int enr = eventNumber + 1; // FORTRAN metric
     hgnf_(&idd, &enr, vector, &ierror);
   }
   
   return (ierror == 0);
}


#endif /* CHBOOK_H */
