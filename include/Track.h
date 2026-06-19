#ifndef TRACK_H
#define TRACK_H

#include "TObject.h"
#include "TString.h"

////////////////////////////////////////////////////////////////////
//Classe per creare traccia e punto di intersezione
////////////////////////////////////////////////////////////////////

class Track: public TObject {

public:
  //costruttori e distruttore
  Track();
  Track(double theta1, double phi1);
  Track(double x, double y, double z);
  virtual  ~Track();
  
  //get
  double Gettheta() const;
  double Getphi() const;
   
  bool Intersection(double x0, double y0, double z0, double &X, double &Y, double &Z, double R, double L); //intersezione con cilindro raggio R 
  void Intersection0(double x0, double y0, double z0, double &Z); //intersezione tracklet con piano x=0
  void MultipleScattering(bool scattering, double sigma=0.001); //multiple scatter (bool=true per accenderlo)
 
private:
  double Parameter(double x0, double y0, double R); //calcolo di parametro t per valutare intersezione
  double fC1, fC2, fC3;//coseni direttori
  double fTheta,fPhi;//theta e phi di retta
 
  ClassDef(Track,1)
};

#endif
