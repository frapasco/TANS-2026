// root includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>

void SetmyIsto(TH1D *isto, int mstyle, int mcolor, int lcolor, TString tit,
               TString xa, TString ya, double xm, double xM, double ym,
               double yM);

void analysis() {
  TFile *hfile = TFile::Open("istogrammi.root",
                             "READ"); // apertura e verifica esistenza file
  if (hfile == NULL) {
    cout << "Errore, non esiste il file" << endl;
    return;
  }
  hfile->ls();             // lettura contenuto file
  gStyle->SetOptFit(1111); // opzioni inserimento informazioni fit in grafici

  TH1D *restot = (TH1D *)hfile->Get(
      "isto_restot"); // estrazione istogramma e verifica sua esistenza
  if (restot == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(
      restot, 4, 4, 4, "residui con tutte le molteplicita'",
      "z_{rec}-z_{true} [cm]", "numero di eventi", -0.08, 0.08, 0,
      1010); // setting di caratteristiche istogramma (stile e colore marker,
             // colore linea, titolo e titoli assi, range in x e y)
  restot->Fit("gaus");
  restot->Draw("pe");

  TH1D *res1 = (TH1D *)hfile->Get(
      "histo_res1"); // estrazione istogramma e verifica sua esistenza
  if (res1 == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res1, 4, 4, 4, "residui con 3<=molteplicita'<=5",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.08, 0.08, 0,
            150); // setting di caratteristiche istogramma (stile e colore
                  // marker, colore linea, titolo e titoli assi)
  res1->Fit("gaus");
  res1->Draw("pe");

  TH1D *res2 = (TH1D *)hfile->Get(
      "histo_res2"); // estrazione istogramma e verifica sua esistenza
  if (res2 == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res2, 4, 4, 4, "residui con 45<=molteplicita'<=55",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.03, 0.03, 0,
            120); // setting di caratteristiche istogramma (stile e colore
                  // marker, colore linea, titolo e titoli assi)
  //  SetmyIsto(res2, 4,4,4,"residui con 45<=molteplicita'<=55",
  //  "z_{rec}-z_{true} [cm]", "numero di eventi", -0.03, 0.03, 0, 90);
  //  //setting di caratteristiche istogramma (stile e colore marker, colore
  //  linea, titolo e titoli assi)
  res2->Fit("gaus");
  res2->Draw("pe");

  TH1D *eff1m = (TH1D *)hfile->Get(
      "isto_eff1sigma"); // estrazione istogramma e verifica sua esistenza
  if (eff1m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(eff1m, 4, 4, 4,
            "efficienza (molteplicita'): z_{true} entro 1 #sigma da 0",
            "molteplicità", "efficienza", 0, 55, 0.7,
            1.05); // setting di caratteristiche istogramma (stile e colore
                   // marker, colore linea, titolo e titoli assi)
  eff1m->Draw("pe");

  TH1D *eff3m = (TH1D *)hfile->Get(
      "isto_eff3sigma"); // estrazione istogramma e verifica sua esistenza
  if (eff3m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; //creazione nuova canvas
  SetmyIsto(eff3m, 4, 4, 4,
            "efficienza (molteplicita'): z_{true} entro 3 #sigma da 0",
            "molteplicità", "efficienza", 0, 55, 0.7,
            1.05); // setting di caratteristiche istogramma (stile e colore
                   // marker, colore linea, titolo e titoli assi)
  eff3m->Draw("pe");

  // un paio di istogrammi di residui bin per bin per vedere se il programma ha
  // funzionato correttamente

  TH1D *res1m = (TH1D *)hfile->Get(
      "isto_residui1[5]"); // estrazione istogramma e verifica sua esistenza
  if (res1m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res1m, 4, 4, 4,
            "residui con 4.5<=molteplicita'<=5.5: z_{true} entro 1 #sigma da 0",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.06, 0.06, 0,
            50); // setting di caratteristiche istogramma (stile e colore
                 // marker, colore linea, titolo e titoli assi)
  res1m->Draw("pe");

  TH1D *res2m = (TH1D *)hfile->Get(
      "isto_residui1[13]"); // estrazione istogramma e verifica sua esistenza
  if (res2m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(
      res2m, 4, 4, 4,
      "residui con 20.5<=molteplicita'<=25.5: z_{true} entro 1 #sigma da 0",
      "z_{rec}-z_{true} [cm]", "numero di eventi", -0.04, 0.04, 0,
      90); // setting di caratteristiche istogramma (stile e colore marker,
           // colore linea, titolo e titoli assi)
  res2m->Draw("pe");

  TH1D *res3m = (TH1D *)hfile->Get(
      "isto_residui3[5]"); // estrazione istogramma e verifica sua esistenza
  if (res3m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res3m, 4, 4, 4,
            "residui con 4.5<=molteplicita'<=5.5: z_{true} entro 3 #sigma da 0",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.06, 0.06, 0,
            60); // setting di caratteristiche istogramma (stile e colore
                 // marker, colore linea, titolo e titoli assi)
  res3m->Draw("pe");

  TH1D *res4m = (TH1D *)hfile->Get(
      "isto_residui1[13]"); // estrazione istogramma e verifica sua esistenza
  if (res4m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(
      res4m, 4, 4, 4,
      "residui con 20.5<=molteplicita'<=25.5: z_{true} entro 3 #sigma da 0",
      "z_{rec}-z_{true} [cm]", "numero di eventi", -0.04, 0.04, 0,
      100); // setting di caratteristiche istogramma (stile e colore marker,
            // colore linea, titolo e titoli assi)
  res4m->Draw("pe");

  TH1D *risol1m = (TH1D *)hfile->Get(
      "isto_ris1sigma"); // estrazione istogramma e verifica sua esistenza
  if (risol1m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(risol1m, 4, 4, 4,
            "risoluzione (molteplicita'): z_{true} entro 1 #sigma da 0",
            "molteplicità", "risoluzione [cm]", 0, 55, 0.005,
            0.04); // setting di caratteristiche istogramma (stile e colore
                   // marker, colore linea, titolo e titoli assi)
  risol1m->Draw("pe");

  TH1D *risol3m = (TH1D *)hfile->Get(
      "isto_ris3sigma"); // estrazione istogramma e verifica sua esistenza
  if (risol3m == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(risol3m, 4, 4, 4,
            "risoluzione (molteplicita'): z_{true} entro 3 #sigma da 0",
            "molteplicità", "risoluzione [cm]", 0, 55, 0.005,
            0.04); // setting di caratteristiche istogramma (stile e colore
                   // marker, colore linea, titolo e titoli assi)
  risol3m->Draw("pe");

  //////////////////////////////////////////////////////////////////////////////////////////

  TH1D *eff1z = (TH1D *)hfile->Get(
      "isto_effztrue"); // estrazione istogramma e verifica sua esistenza
  if (eff1z == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(eff1z, 4, 4, 4, "efficienza (z_{true})", "z_{true} [cm]",
            "efficienza", -18, 18, 0.82,
            1); // setting di caratteristiche istogramma (stile e colore marker,
                // colore linea, titolo e titoli assi)
  eff1z->Draw("pe");

  // un paio di istogrammi residui bin per bin per vadere se sono "sensati"

  TH1D *res1z = (TH1D *)hfile->Get(
      "isto_residuiztrue[2]"); // estrazione istogramma e verifica sua esistenza
  if (res1z == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res1z, 4, 4, 4, "residui con -9 cm<z_{true}<=-7 cm",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.06, 0.06, 0,
            70); // setting di caratteristiche istogramma (stile e colore
                 // marker, colore linea, titolo e titoli assi)
  res1z->Draw("pe");

  TH1D *res2z = (TH1D *)hfile->Get(
      "isto_residuiztrue[6]"); // estrazione istogramma e verifica sua esistenza
  if (res2z == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(res2z, 4, 4, 4, "residui con -1 cm<z_{true}<=1 cm",
            "z_{rec}-z_{true} [cm]", "numero di eventi", -0.06, 0.06, 0,
            180); // setting di caratteristiche istogramma (stile e colore
                  // marker, colore linea, titolo e titoli assi)
  res2z->Draw("pe");

  TH1D *risol1z = (TH1D *)hfile->Get(
      "isto_risztrue"); // estrazione istogramma e verifica sua esistenza
  if (risol1z == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  new TCanvas; // creazione nuova canvas
  SetmyIsto(risol1z, 4, 4, 4, "risoluzione (z_{true})", "z_{true} [cm]",
            "risoluzione [cm]", -18, 18, 0.018,
            0.027); // setting di caratteristiche istogramma (stile e colore
                    // marker, colore linea, titolo e titoli assi)
  risol1z->Draw("pe");
}

void SetmyIsto(TH1D *isto, int mstyle, int mcolor, int lcolor, TString tit,
               TString xa, TString ya, double xm, double xM, double ym,
               double yM) {
  isto->SetMarkerStyle(mstyle);
  isto->SetMarkerColor(mcolor);
  isto->SetLineColor(lcolor);
  isto->SetTitle(tit);
  isto->GetXaxis()->SetTitle(xa);
  isto->GetYaxis()->SetTitle(ya);
  isto->GetXaxis()->SetRangeUser(xm, xM); // range di asse x
  isto->GetYaxis()->SetRangeUser(ym, yM); // range di asse y
}
