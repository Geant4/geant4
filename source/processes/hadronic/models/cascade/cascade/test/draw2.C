{
gROOT->Reset();

#include <iostream.h>
#include <ctype.h>

ifstream in;
int nlines = 0, i=0;
double x, y, *u, *v;
char array1[128], array2[128], array3[128], array4[128], c;

in.open("draw2.out", ios::in); // data must be in current directory

//let's read the number of rows
while(1){
  in >> x >> y;
  if(!in.good()) break;
  nlines++;
}

in.close();
in.open("draw2.out", ios::in);

cout << "The name of the canvas:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    array1[i]=c;
    i++;
  }
  if(i<128)
    array1[i]='\0';

cout << "The title of the figure:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    array2[i]=c;
    i++;
  }
  if(i<128)
    array2[i]='\0';

cout << "The title of the x-axis:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    array3[i]=c;
    i++;
  }
  if(i<128)
    array3[i]='\0';

cout << "The title of y-axis:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    array4[i]=c;
    i++;
  }
  if(i<128)
    array4[i]='\0';

c1=new TCanvas("c1", array1);
c1->SetFillColor(0);

c1->SetGridx();
c1->SetGridy();

u = new double[nlines];
v = new double[nlines];

i=0;
while (1) {
  in >> u[i] >> v[i];
  if (!in.good()) break;
  i++;
}

gr1=new TGraph(nlines, u, v);
gr1->SetFillColor(0);
//gr1->SetMarkerStyle(20);
gr1->SetTitle(array2);
//gr1->Draw();
gr1->Draw("ACP"); 
gr1->GetXaxis()->SetTitle(array3);
gr1->GetYaxis()->SetTitle(array4);
gr1->GetXaxis()->CenterTitle();
gr1->GetYaxis()->CenterTitle();
c1->SetBorderSize(0);

in.close();
}

/*
   c1->SetFillColor(42);
   c1->SetGridx();
   c1->SetGridy();

  gr = new TGraph(n,x,y);
   gr->SetFillColor(19);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(20);
   gr->SetTitle("a simple graph");
   gr->Draw("ACP");
*/
