#ifndef VERTEX_H
#define VERTEX_H

#include "TObject.h"
#include "TString.h"

//////////////////////////////////////////////////////////////////////////////////
//Classe per generare vertice primario e direzione prodotti
//////////////////////////////////////////////////////////////////////////////////

class Vertex : public TObject {  

public:  
  //costruttori e distruttore
  Vertex();                
  Vertex(TString multv, TH1D *multiHisto, double sigmax, double sigmay, double sigmaz, int u1=1, int u2=90); //multv=scelta distribuzione molteplicità
  Vertex(const Vertex &v); //copy constructor
  virtual ~Vertex();
            
  void Initialdir(double min, double max, TH1D *etaHisto);//estrazione direzione prodotti di collisione
  
  //get
  double Getx() const;
  double Gety() const;
  double Getz() const;
  double Getm() const;
  double Getphi() const; 
  double Gettheta() const; 
  
  //setters
  void SetPhi(double phi);
  void SetTheta(double theta);
  void SetMulti(double multi);

private:
  double fX, fY, fZ;//coordinate vertice
  int fMulti;//numero di particelle cariche da collisione (molteplicità)
  double fPhi,fTheta;//angoli polare e azimutale di prodotti
  
  double Var(double s); //estrazione gaussiana con Box-Muller
  void Multunif(int u1, int u2); //estrazione molteplicità uniforme fra u1 e u2
 
  ClassDef(Vertex,1);
};  
        
#endif
