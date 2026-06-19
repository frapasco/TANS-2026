#ifndef INTPOINT_H
#define INTPOINT_H

#include "TObject.h"

////////////////////////////////////////////////////////////////////
//Classe per salvare le coordiante del punto di intersezione
////////////////////////////////////////////////////////////////////

class Intpoint: public TObject {

public:
  //costruttori e distruttori
  Intpoint();
  Intpoint(double X, double Y, double Z, int label=-1);
  Intpoint(double R, double L, int label=-1);
  virtual ~Intpoint();
   
  //get
  double Getx() const;
  double Gety() const;
  double Getz() const;
  double Getphi() const;
  int Getlabel() const;
  
  void Smearing(double R, double sigz=0.012,double sigt=0.003); //smearing di punti di hit

private:
  double fX, fY, fZ;//coordinate del punto di intersezione
  int fLabeltrack;//numero della traccia a cui appartiene l'intersezione
  ClassDef(Intpoint,1)
};

#endif
