#ifndef Histo_h
#define Histo_h

#include "globals.hh"
#include "g4std/iostream"

// Due to RW-STL transition problems this simpl histo class works only
// with C++-arrays

class odHisto
{
public:

  odHisto()
    {
      content=0;
      overflow=0;
      underflow=0;
      sumw=0;
      sumx=0;
      sumx2=0;
      lowerbound=0;
      binwidth=0;
      entries=0;
      nbins=0;
    }
      
  
  odHisto(int,double,double);

  odHisto(const odHisto&);

  odHisto& operator=(const odHisto&);

  ~odHisto()
    {
      delete[] content;
    }
  
  void accumulate(double x,double weight=1.0);
  
  double Overflow() const
    {
      return overflow;
    }

  double Underflow() const
    {
      return underflow;
    }

  double Sum() const
    {
      double s=0;
      for(int i=0;i<nbins;i++)
	s+=content[i];
      return s;
    }

  double Average() const
    {
      return sumx/sumw;
    }

  double RMS() const
    {
      double av2=Average();
      return sqrt(sumx2/sumw-av2*av2);
    }

  int Entries() const
    {
      return entries;
    }

  void output(G4std::ostream& o) const
    {
      //The statistical information
      o << "# Bins " << nbins << G4endl;
      o << "# Entries " << entries << G4endl;
      o << "# Overflow " << overflow << G4endl;
      o << "# Underflow " << underflow << G4endl;
      o << "# Sum " << sumw << G4endl;
      o << "# Average " << Average() << G4endl;
      o << "# RMS " << RMS() << G4endl;
      int i;
      for(i=0;i<nbins;i++)
	o << i*binwidth+lowerbound << ' ' << content[i] << G4endl;
      o << i*binwidth+lowerbound << ' ' << content[i-1] << G4endl;
    }
      
  void input(G4std::istream& in)
    {
      char c,dummy[80];
      double average,rms,bla;
      //The statistical information
      in >> c >> dummy >> nbins;
      in >> c >> dummy >> entries;
      in >> c >> dummy >> overflow;
      in >> c >> dummy >> underflow;
      in >> c >> dummy >> sumw;
      in >> c >> dummy >> average;
      in >> c >> dummy >> rms;
      if(content) delete[] content;
      content=new double[nbins];
      if(nbins>0)
	in >> lowerbound >> content[0];
      if(nbins>1)
	{
	  in >> bla >> content[1];
	  binwidth=bla-lowerbound;
	}
      for(int i=2;i<nbins;i++)
	in >> bla >> content[i];
      sumx=average*sumw;
      sumx2=(rms*rms+average*average)*sumw;
    }

  double compare(const odHisto&) const;

private:

  double *content;
  double overflow;
  double underflow;
  double sumw;
  double sumx;
  double sumx2;
  double lowerbound;
  double binwidth;
  int entries;
  int nbins;
};

#endif





