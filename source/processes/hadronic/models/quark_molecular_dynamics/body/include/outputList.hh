#include <iostream.h>
#include "newBinning.hh"

class outputList
{
  int N;
  ostream** files;
  bool* erase;
public:
  outputList(int n);
  outputList(int n,char* file);
  ~outputList();
  int Number() const { return N; }
  void add(int i,ostream*);
  ostream& operator[](int i) { return *files[i]; }
};

class PlotList : public outputList
{
  Binning** plots;
public:
  PlotList(int n) : outputList(n),plots(new Binning*[n]) {}
  PlotList(int n,char* file) : outputList(n,file),plots(new Binning*[n]) {}
  ~PlotList() { delete [] plots; }
  void add(int i,Binning& x,ostream* f = 0) { plots[i] = &x; if (f) outputList::add(i,f); }
  Binning& content(int i) { return *plots[i]; }
};

