#include <iostream.h>
#include <string.h>

int main(){

gROOT->Reset();

 ifstream in;
 int numberOfVectors, size1, size2, size3, size4;
 double *index1, *index2, *index3, *index4, *xsec1, *xsec2, *xsec3, *xsec4;
 char name1[128], name2[128], name3[128], name4[128], c;

 in.open("xsecs.out", ios::in);

  cout << "How many data vectors?" << endl;
  cin >> numberOfVectors;

  if(numberOfVectors==3){
    cout << "Give the sizes of the data vectors:" << endl;
    cin >> size1 >> size2 >> size3;
   // cout << size1 << size2 << size3;
    index1 = new double[size1];
    index2 = new double[size2];
    index3 = new double[size3];
    xsec1 = new double[size1];
    xsec2 = new double[size2];
    xsec3 = new double[size3];

    int i;
    for(i=0; i < size1; i++) index1[i]=i;
    for(i=0; i < size2; i++) index2[i]=i;
    for(i=0; i < size3; i++) index3[i]=i;

    for(i=0; i < size1; i++) in >> xsec1[i];
    for(i=0; i < size2; i++) in >> xsec2[i];
    for(i=0; i < size3; i++) in >> xsec3[i];

  cout << "Give the first title:" << endl;      
  
  c=getchar();
  i=0;
  while((c=getchar())!='\n' && i<127){
    name1[i]=c;
    i++;
  }
  name1[i]='\0';

  cout << "Give the second title:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){   
    name2[i]=c;
    i++;
  }
  name2[i]='\0';

  cout << "Give the third title:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){   
    name3[i]=c;
    i++;
  }
  name3[i]='\0';

  c1=new TCanvas("c1", "cross sections in HETC");
  c1->Divide(1,3);

  c1->cd(1);
  gr1=new TGraph(size1,index1,xsec1);
  gr1->SetTitle(name1);
  gr1->Draw("ACP");
  gr1->GetXaxis()->SetTitle("Index");
  gr1->GetYaxis()->SetTitle("Cross section / barn");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();

  c1->cd(2);
  gr2=new TGraph(size2,index2,xsec2);
  gr2->SetTitle(name2);
  gr2->Draw("ACP");
  gr2->GetXaxis()->SetTitle("Index");
  gr2->GetYaxis()->SetTitle("Cross section / barn");
  gr2->GetXaxis()->CenterTitle();
  gr2->GetYaxis()->CenterTitle();

  c1->cd(3);
  gr3=new TGraph(size3,index3,xsec3);
  gr3->SetTitle(name3);
  gr3->Draw("ACP");
  gr3->GetXaxis()->SetTitle("Index");
  gr3->GetYaxis()->SetTitle("Cross section / barn");
  gr3->GetXaxis()->CenterTitle();
  gr3->GetYaxis()->CenterTitle();
 }
  if(numberOfVectors==4){
    cout << "Give the sizes of the data vectors:" << endl;
    cin >> size1 >> size2 >> size3 >> size4;
    cout << size1 <<" "<< size2 <<" "<< size3 <<" "<< size4 << endl;   

    index1 = new double[size1];
    index2 = new double[size2];
    index3 = new double[size3];
    index4 = new double[size4];
    xsec1 = new double[size1];
    xsec2 = new double[size2];
    xsec3 = new double[size3];
    xsec4 = new double[size4];
 
    int i;

    for(i=0; i < size1; i++) index1[i]=i;  
    for(i=0; i < size2; i++) index2[i]=i;
    for(i=0; i < size3; i++) index3[i]=i;
    for(i=0; i < size4; i++) index4[i]=i;

    
    for(i=0; i < size1; i++){ in >> xsec1[i]; cout << "xsec1: " << xsec1[i] << endl;}
    for(i=0; i < size2; i++) in >> xsec2[i];
    for(i=0; i < size3; i++) in >> xsec3[i];
    for(i=0; i < size4; i++){ in >> xsec4[i]; cout << "xsec4: " << xsec4[i] << endl; }
  
   
    cout << "Give the first title:" << endl;

  c=getchar();  
  i=0;
  while((c=getchar())!='\n' && i<127){
    name1[i]=c;
    i++;
  }
  name1[i]='\0';
  
  cout << name1 << endl;
  
  cout << "Give the second title:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    name2[i]=c;
    i++;
  }
  name2[i]='\0';
   
  cout << "Give the third title:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    name3[i]=c;
    i++;
  }
  name3[i]='\0';

  cout << "Give the fourth title:" << endl;
  i=0;
  while((c=getchar())!='\n' && i<127){
    name4[i]=c;
    i++;
  }
  name4[i]='\0';

  
  c1=new TCanvas("c1", "cross sections in HETC");
  c1->Divide(1,4);
  
  c1->cd(1); 
  gr1=new TGraph(size1,index1,xsec1); 
  gr1->SetTitle(name1);
  gr1->Draw("ACP");
  gr1->GetXaxis()->SetTitle("Index");
  gr1->GetYaxis()->SetTitle("Cross section / barn");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();

  c1->cd(2);
  gr2=new TGraph(size2,index2,xsec2);
  gr2->SetTitle(name2);
  gr2->Draw("ACP");
  gr2->GetXaxis()->SetTitle("Index");
  gr2->GetYaxis()->SetTitle("Cross section / barn");
  gr2->GetXaxis()->CenterTitle();
  gr2->GetYaxis()->CenterTitle();

  c1->cd(3);
  gr3=new TGraph(size3,index3,xsec3);
  gr3->SetTitle(name3);
  gr3->Draw("ACP");
  gr3->GetXaxis()->SetTitle("Index");
  gr3->GetYaxis()->SetTitle("Cross section / barn");
  gr3->GetXaxis()->CenterTitle();
  gr3->GetYaxis()->CenterTitle();

  c1->cd(4);
  gr4=new TGraph(size4,index4,xsec4);
  gr4->SetTitle(name4);
  gr4->Draw("ACP");
  gr4->GetXaxis()->SetTitle("Index");
  gr4->GetYaxis()->SetTitle("Cross section / barn");
  gr4->GetXaxis()->CenterTitle();
  gr4->GetYaxis()->CenterTitle();
    
  }
 
  return 0;
}


