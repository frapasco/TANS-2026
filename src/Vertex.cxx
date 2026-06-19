#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "Vertex.h"

using namespace std;
using namespace TMath ;

ClassImp(Vertex)

//Nota: il file contenente le distribuzioni per eta e molteplicità (e i relativi istogrammi) viene manipolato dalle funzioni interne alla classe, 
//cercando così di ottenere la maggiore protezione possibile del suo contenuto dall'utente (minore rischio che una modifica/un errore da parte di un utente rovini 
//il file e porti a impossibilità (o assurdità) di estrazioni random dagli istogrammi e conseguenti possibili problemi con l'esecuzione del programma
//_____________________________________________________________________________________________
//costruttore di default (vertice=origine, 0 prodotti)
Vertex :: Vertex() : TObject(){  
  fX=0;     
  fY=0;   
  fZ=0;
  fMulti=0;   
}

//_____________________________________________________________________________________________
//costruttore standard 
Vertex :: Vertex(TString multv, TH1D *multiHisto, double sigmax, double sigmay, double sigmaz, int u1, int u2) : TObject(){
  fMulti=0;
  fX=Var(sigmax);    //x,y,z gaussiani attorno a origine
  fY=Var(sigmay);   
  fZ=Var(sigmaz);
  if (multv=="No"||multv=="NO"||multv=="no") this->Multunif(u1, u2); //con stringa si sceglie come estrarre molteplicità (uniforme o funzione)
  else fMulti=(int) multiHisto->GetRandom(); //estrazione random secondo la distribuzione rappresentata dall'istogramma
}

//_____________________________________________________________________________________________
//copy constructor 
Vertex :: Vertex(const Vertex &v) : 
  TObject(v),
  fX(v.fX),
  fY(v.fY),
  fZ(v.fZ),
  fMulti(v.fMulti)
  {
  }
  
//_____________________________________________________________________________________________
//Set delle coordiante (Box-Muller): variazione gaussiana di x, y e z rispetto all'origine
double Vertex::Var(double s){ 
  if(s>0){
  	double u1=gRandom->Rndm(); 
  	double u2=gRandom->Rndm();   
  	return Sqrt(-2*Log(u1))*Cos(2*Pi()*u2)*s; //estrazione random da gaussiana con media=0 (origine), sigma=s
  }
  else{  
    cout<<"Errore: sigma <=0"<<endl; //dovrebbe sempre essere maggiore di 0 per controlli in simulazione, verifica errori di memoria e di battitura e cambiamenti in macro simulazione
    return 0;
  }
}

//_____________________________________________________________________________________________
//Getters
double Vertex :: Getx() const{
	return fX;
}

double Vertex :: Gety() const{
	return fY;
}

double Vertex :: Getz() const{
	return fZ;
}

double Vertex :: Getm() const{
	return fMulti;
}

double Vertex :: Getphi() const{
	return fPhi;
}

double Vertex :: Gettheta() const{
	return fTheta;
}

//_____________________________________________________________________________________________
//Setters
void Vertex :: SetPhi(double phi){
  fPhi=phi;
}

void Vertex :: SetTheta(double theta){
  fTheta=theta;
}

void Vertex :: SetMulti(double multi){
  fMulti=multi;
}


//_______________________________________________________________
//metodo per molteplicità da distribuzione uniforme fra u1 e u2
void Vertex :: Multunif(int u1, int u2){
 if (u2>=u1)  fMulti=u1+(u2-u1)*gRandom->Rndm();
 else {
   cout<<"Errore in scelta parametri"<<endl;
   return;
 }
}


//_______________________________________________________________
//inizializza la direzione iniziale di un prodotto dal vertice primario, min e max per ridurre il range della distribuzione al range di proprio interesse (riducendo il numero di tracce fuori rivelatore a causa del loro theta)
//attenzione in scelta max e min per non perdere una frazione troppo grande della distribuzione (molte estrazioni "a vuoto"); con -2\2 e la distribuzione assegnata, si preservano i significativi picchi centrali e non si perdono troppe estrazioni
void Vertex :: Initialdir(double min, double max, TH1D *etaHisto){ 
  if(max<min){ //scambio estremi in modo che siano in ordine giusto
    double o=max; 
    max=min; 
    min=o;
  }
  if(max==min){
    cout<<"Errore: etamax=etamin"<<endl; 
    return;
  }
  fPhi = gRandom->Rndm()*2*Pi(); //phi uniforme
  double eta=0.;
  do{
    eta =(double) etaHisto->GetRandom();
  }
  while (eta>max||eta<min); //estrae finchè non trova un valore nel range max--min voluto
  fTheta=2*ATan(Exp(-eta)); //theta come funzione di eta
} 

//_______________________________________________________
//distruttore      
Vertex :: ~Vertex(){  
}
