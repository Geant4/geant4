#include <string.h>
#include <math.h>

#include "Histo.hh"

odHisto::odHisto(int bins,double l,double h): overflow(0.0),underflow(0.0),
sumw(0.0),sumx(0.0),sumx2(0.0),entries(0),nbins(bins)
{
  lowerbound=l;
  binwidth=(h-l)/bins;
  content=new double[nbins];
}

odHisto::odHisto(const odHisto&h)
{
  overflow=h.overflow;
  underflow=h.underflow;
  sumw=h.sumw;
  sumx=h.sumx;
  sumx2=h.sumx2;
  entries=h.entries;
  nbins=h.nbins;
  lowerbound=h.lowerbound;
  binwidth=h.binwidth;
  content=new double[nbins];
  memcpy(content,h.content,nbins*sizeof(double));
}

odHisto& odHisto::operator=(const odHisto&h)
{
  overflow=h.overflow;
  underflow=h.underflow;
  sumw=h.sumw;
  sumx=h.sumx;
  sumx2=h.sumx2;
  entries=h.entries;
  if(content)
    {
      if(nbins!=h.nbins)
	{
	  delete[] content;
	  content=new double[h.nbins];
	}
    }
  else
    {
      content=new double[h.nbins];
    }
  nbins=h.nbins;
  lowerbound=h.lowerbound;
  binwidth=h.binwidth;
  memcpy(content,h.content,nbins*sizeof(double));
  return *this;
}

void odHisto::accumulate(double x,double weight)
{
  int bin; 

  bin=(int) ((x-lowerbound)/binwidth);
  
  entries++;
  sumx+=x*weight;
  sumx2+=x*x*weight;
  sumw+=weight;

  if(bin<0) 
    underflow+=weight;
  else if(bin>=nbins)
    overflow+=weight;
  else
    content[bin]+=weight;

}
  
double odHisto::compare(const odHisto& ah) const
{
  if(nbins!=ah.nbins) return -1.0;

  if(fabs(lowerbound-ah.lowerbound)>fabs(lowerbound)*1.e-5) return -1.0;
  
  if(fabs(binwidth-ah.binwidth)>binwidth*1.e-5) return -1.0;

  // Ok same booking

  int i;
  double dbin=0.;
  double d,e;

  for(i=0;i<nbins;i++)
    {
      if(content[i]>0.)
	e=content[i];
      else
	e=1.;
      e=sqrt(e);
      d=ah.content[i]-content[i];
      dbin+=d*d/e;
    }

  return dbin/nbins;
}
  










