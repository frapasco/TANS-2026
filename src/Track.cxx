#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "Track.h"

using namespace TMath; 
using namespace std;

ClassImp(Track)

//_____________________________________________________________________________________________
//costruttore di default 
Track :: Track() : TObject(){
  fC1=0;
  fC2=0;
  fC3=0;
}

//_____________________________________________________________________________________________
//costruttore standard 
Track :: Track (double theta1, double phi1) : TObject(){
  fTheta=theta1;
  fPhi=phi1;
  fC1=Sin(fTheta)*Cos(fPhi);
  fC2=Sin(fTheta)*Sin(fPhi);
  fC3=Cos(fTheta);	
}

//_____________________________________________________________________________________________
//costruttore per creare retta per la ricostruzione di tracklet
Track :: Track (double x, double y, double z) : TObject(){
  double norm =1/Sqrt(x*x+y*y+z*z); 
  fC1= x*norm;
  fC2= y*norm;
  fC3= z*norm;	
}  

//_____________________________________________________________________________________
//Get per theta e phi
double Track :: Gettheta () const{ 
  return fTheta;
}	

double Track :: Getphi () const{ 
  return fPhi;
}  

//_____________________________________________________________________________________________
//Determinazione parametro t della retta
double Track :: Parameter (double x0, double y0, double R){ 
  double prod = x0*fC1+y0*fC2;
  double C = fC1*fC1+fC2*fC2;
  double Delta = (prod*prod)-C*((x0*x0)+(y0*y0)-(R*R));
  
  if (Delta<0){ //delta<0 se (x0,y0,z0) esterno a cilindro raggio R
   	cout<<"Errore, discriminante negativo."<<endl;	
	  return 0;
  }
  double t;
  double t1 = (-prod+Sqrt(Delta))/C;
  double t2 = (-prod-Sqrt(Delta))/C;

  if(t1>0) t=t1; //selezione semiretta voluta (t>0) 
  else t=t2;

  return t;
}

//_____________________________________________________________________________________________
//Valutazione punto di intersezione fra retta e cilindro raggio R 
bool Track :: Intersection(double x0, double y0, double z0, double &X, double &Y, double &Z, double R, double L){  
  double tt;
  tt =  Parameter(x0, y0, R); //valutazione tt (parametro retta)
  
  //coordinate punto intersezione
  X=x0+fC1*tt; 
  Y=y0+fC2*tt;
  Z=z0+fC3*tt;
  
  //per la simulazione, si valuta se Z del punto di hit rientra nelle dimensioni dei rivelatori
  if(Abs(Z)>L/2){ 
    return false;
  }
  else return true;
}

//___________________________________________________________________________________________
//intersezione tracklet con piano x=0
void Track::Intersection0(double x0, double y0, double z0, double &Z){
  double tt=-x0/fC1; //mi serparametro retta che corrisponde al suo punto di intersezione col piano x=0 (tt=(0-x0)/fC1)	
  Z=z0+fC3*tt;	//coordinata Z del punto di intersezione della tracklet col piano x=0
}

//_____________________________________________________________________________________________
//multiple scattering
void Track:: MultipleScattering(bool scattering, double sigma){
  if(scattering==true){	 //accensione/spegnimento multiple scatter con bool
	  if (sigma<=0) {cout<<"Errore nel settare sigma di gaussiana"<<endl; return ;} 
	  double thetaprimo = gRandom->Gaus(0,sigma); //theta e phi di scatter rispetto alla direzione della traccia prima del multiple scatter: theta gaussiano, phi uniforme
	  double phiprimo = 2*Pi()*gRandom->Rndm(); 

	  //matrice di trasformazione fra sistemi riferimento (SR del laboratorio: fasci lungo z; SR di traccia di un prodotto: z'=direzione traccia)
	  double M[3][3];
	  M[0][0]=-Sin(fPhi);
	  M[0][1]=-Cos(fTheta)*Cos(fPhi);
	  M[0][2]=Sin(fTheta)*Cos(fPhi);
	  M[1][0]=Cos(fPhi);
	  M[1][1]=-Cos(fTheta)*Sin(fPhi);
	  M[1][2]=Sin(fTheta)*Sin(fPhi);
	  M[2][0]=0;
	  M[2][1]=Sin(fTheta);
	  M[2][2]=Cos(fTheta);

    //coordinate del versore del nuovo vettore (direzione traccia post scattering) nel sistema di riferimento della traccia (coseni direttori)
  	double usecondo[3];
	  usecondo[0]=Sin(thetaprimo)*Cos(phiprimo);
	  usecondo[1]=Sin(thetaprimo)*Sin(phiprimo);
	  usecondo[2]=Cos(thetaprimo);

	  //coordinate u[i] del versore del nuovo vettore (traccia post scattering) nel sistema di riferimento del laboratorio (coseni direttori)
	  double u[3];
	  for(int i=0;i<3;i++){
	    u[i]=0; 
	    for(int j=0;j<3;j++){
		  u[i]+=M[i][j]*usecondo[j];
      }
    }
	  //theta e phi post scatter
    if (u[0]==0) {
	   if (u[1]==0) {cout<<"Errore nei coseni direttori: ATan non definita"<<endl; return;}
	   if (u[1]>0) fPhi=Pi()/2;
     if (u[1]<0) fPhi=3*Pi()/2;
	  }
	  if (u[0]>0){
	   if (u[1]>=0)	fPhi=ATan(u[1]/u[0]);
	   else fPhi=ATan(u[1]/u[0])+2*Pi();
	  }
	  if (u[0]<0){
	   if(u[1]>0) fPhi=ATan(u[1]/u[0])+Pi(); 
	   else fPhi=ATan(u[1]/u[0])+Pi(); 
	  } 
	  if(u[2]>0) fTheta=ATan(Sqrt(u[0]*u[0]+u[1]*u[1])/u[2]); 
	  else {
	    if(u[2]<0) fTheta=ATan(Sqrt(u[0]*u[0]+u[1]*u[1])/u[2])+Pi();
	    else fTheta=Pi()/2;
	  }
  }
}

//_____________________________________________________________________________________________
//distruttore
Track :: ~Track(){	 
}
