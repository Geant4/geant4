#include "g4std/fstream"
#include "Histo.hh"


int main(int argc,char* argv[])
{
  
  G4std::ifstream in(argv[1]),comp(argv[2]);

  odHisto hin,hcomp;
  char c,bla[80];
  int nhisto[2];
  
  in >> c >> bla >> nhisto[0];
  comp >> c >> bla >> nhisto[1];

  if(nhisto[0]!=nhisto[1]) exit(-1);

  for(int i=0;i<nhisto[0];i++)
    {
      hin.input(in);
      hcomp.input(comp);
      if(hin.compare(hcomp)>1.) exit(-1);
    }

  exit(0);
}
  
