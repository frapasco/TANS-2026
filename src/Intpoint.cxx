#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include "Intpoint.h"

using namespace TMath; 
using namespace std;

ClassImp(Intpoint)


#include "TMath.h"
#include "Intpoint.h"

using namespace TMath; 
using namespace std;

ClassImp(Intpoint)

//_____________________________________________________________________________________________
//default constructor
Intpoint :: Intpoint() : TObject(){
  fX=0;
  fY=0;
  fZ=0;
  fLabeltrack=0;
} 

//______________________________________________________________________________________________
//costruttore standard, usato per memorizzare le coordinate dei punti di intersezione con beam pipe e rivelatori
Intpoint :: Intpoint (double X, double Y, double Z, int label) : TObject(){
  fX=X;
  fY=Y;
  fZ=Z;
  fLabeltrack=label;
}

//______________________________________________________________________________________________
//costruttore per generare noise  uniforme sui due rivelatori
Intpoint :: Intpoint (double R, double L, int label) : TObject(){
  fLabeltrack=label;
  double phi = gRandom->Rndm()*2*Pi();  
  fX=R*Cos(phi);
  fY=R*Sin(phi);
  fZ=-L/2+L*gRandom->Rndm();
}

//______________________________________________________________________________________________
double Intpoint::Getx() const{ 
  return fX;
}

double Intpoint::Gety() const{ 
  return fY;
}

double Intpoint::Getz() const{ 
  return fZ;
}

//Get phi: angolo phi del punto di intersezione in coordinate cilindriche
double Intpoint::Getphi() const{ 
  if (fX==0&&fY==0) {cout<<"Errore: x=0 e y=0 in intersezione"<<endl; return -1000;}
  if (fX==0&&fY>0) return Pi()/2;
  if (fX==0&&fY<0) return 3*Pi()/2;
  if (fX>0&&fY>=0) return ATan(fY/fX);
  if (fX>0&&fY<0) return ATan(fY/fX)+Pi()*2;
  if (fX<0&&fY>0) return ATan(fY/fX)+Pi();
  if (fX<0&&fY<=0) return ATan(fY/fX)+Pi();
  cout<<"Errore in Intpoint::Getphi()"<<endl;
  return 0;  
}
 
int Intpoint::Getlabel() const{ 
  return fLabeltrack;
} 

//______________________________________________________________________________________________
//funzione per smearing 
void Intpoint::Smearing(double R, double sigz, double sigt){
  if (sigz<=0||sigt<=0) {cout<<"Errore inizializzazione parametri"<<endl; return;}
  double sz = sigz; //cm. Std dev ricostruzione in z
  double st = sigt; //cm. Std dev ricostruzione trasversale (R*phi)
  
  double zrec=gRandom->Gaus(0,sz); //smearing in z
  double arec=gRandom->Gaus(0,st); // smearing trasversale
  //double zrec=Sqrt(-2*Log(gRandom->Rndm()))*Cos(2*Pi()*gRandom->Rndm())*sz; //smearing in z con Box-Muller
  //double arec=Sqrt(-2*Log(gRandom->Rndm()))*Cos(2*Pi()*gRandom->Rndm())*st; //smearing trasversale con Box-Muller
  
  double phirec = this->Getphi()+arec/R; //angolo phi dopo smearing
  
  //coordinate dopo smearing
  fX=R*Cos(phirec);
  fY=R*Sin(phirec);
  fZ=fZ+zrec;
 
}

//______________________________________________________________________________________________
//distruttore
Intpoint:: ~Intpoint(){ 
}
