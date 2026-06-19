// root includes
#include "TClonesArray.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include <TH1D.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// headers includes
#include "Intpoint.h"
#include "Track.h"
#include "Vertex.h"

using namespace std;
// tutte le dimensioni si intendono in cm
// variabili globali
const double kRbp = 3.04; // raggio beam pipe (a metà spessore)
const double kR1 = 4.01;  // raggio primo rivelatore (T1) (a metà spessore)
const double kR2 = 7.01;  // raggio secondo rivelatore (T2) (a metà spessore)
const double kL = 27;     // lunghezza rivelatori

const double kSigmax = 0.01; // spread gaussiano (sigma) delle coordinate del
                             // vertice attorno a (0,0,0)
const double kSigmay = 0.01;
const double kSigmaz = 5.3;

const int kEvents = 100000; // numero vertici primari (PV) generati
const double kEtamin =
    -2.; // limiti in pseudorapidità entro cui si generano le tracce dei
         // prodotti (limiti imposti su istogramma di distribuzione eta)
const double kEtamax = +2.;

const TString kMul =
    "sì"; // scelta distribuzione di molteplicità("sì"=distribuzione assegnata;
          // "No","NO","no"=uniforme)
const bool kMs = true;         // multiple scattering (true=attivo)
const int kStampasim = 10000;  // ogni quanti eventi fa stampe
const unsigned int kSeed = 18; // seed per estrazione random (TRandom3)

const TString kFile = "simulazione.root"; // file per salvare tree

////////////////////////////////////////////////////////////////////
// Funzione per generare i vertici primari e memorizzare i punti
// di impatto delle particelle sui diversi piani di rivelazione.
////////////////////////////////////////////////////////////////////

void simulazione() {
  if (kR1 <= kRbp || kR2 <= kR1 || kRbp <= 0 || kEvents <= 0 ||
      kEtamax <= kEtamin || kL <= 0 || kSigmaz <= 0 || kSigmax <= 0 ||
      kSigmay <= 0 || kStampasim <= 0) {
    cout << "Errore nel setting dei parametri" << endl;
    return;
  }
  TRandom3 *rand = new TRandom3(kSeed);
  delete gRandom;
  gRandom = rand;

  TString filename = "kinem.root"; // nome file (e istogrammi al suo interno)
                                   // per distribuzioni eta e molteplicità
  TString muldis = "hmul";
  TString etadis = "heta";
  char a;
  if (kMul == "No" || kMul == "NO" || kMul == "no")
    a = 'N';
  else
    a = 'S';

  // reading the file to extract the kinematic histos
  TFile *f = new TFile(filename);
  if (f == NULL) {
    cout << "Errore, non esiste il file" << endl;
    return;
  }
  TH1D *multiHisto = (TH1D *)f->Get(
      muldis); // estrazione e verifica esistenza istogramma in file
  if (multiHisto == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }
  TH1D *etaHisto =
      (TH1D *)f->Get(etadis); // estrazione istogramma per distribuzione in eta
  if (etaHisto == NULL) {
    cout << "Errore, non esiste l'istogramma" << endl;
    return;
  }

  // creazione array per punti di intersezione
  TClonesArray *bp = new TClonesArray("Intpoint", 100);     // beam pipe
  TClonesArray *clone1 = new TClonesArray("Intpoint", 100); // primo rivelatore
  TClonesArray *clone2 =
      new TClonesArray("Intpoint", 100); // secondo rivelatore
  TClonesArray &intbp = *bp;
  TClonesArray &int1 = *clone1;
  TClonesArray &int2 = *clone2;

  int nparticle; // numero particelle prodotte

  // file per passare dati al codice per la ricostruzione
  FILE *hdata = fopen("data.txt", "w");
  fprintf(hdata, "%f %f %f %f %d %c", (float)kR1, (float)kR2, (float)kL,
          (float)kSigmaz, kEvents, a);
  fclose(hdata);

  // creazione puntatore a oggetto Vertex
  Vertex *vertice = NULL;

  // creazione TFile, TTree e Branches
  TFile *hfile = TFile::Open(kFile, "RECREATE");
  TTree *tree = new TTree("tree", "eventi");
  tree->Branch("Vertext", &vertice); // PV
  tree->Branch("tclonebp", &bp);     // intersezione BP
  tree->Branch("tclone1", &clone1);  // intersezione T1
  tree->Branch("tclone2", &clone2);  // intersezione T2

  // variabili per salvare coordinate hit
  double X1, Y1, Z1;
  double X2, Y2, Z2;
  double Xbp, Ybp, Zbp;
  int lost1, lost2, lostall; // tracce non viste da T1, T2, entrambi
  int vertexnohit = 0;       // numero vertici in totale senza nessun hit

  tree->SetDirectory(
      hfile); // assicurazione che il tree appartenga al file voluto

  // loop sul numero di vertici
  for (int i = 0; i < kEvents; i++) {
    lost1 = 0;
    lost2 = 0;
    lostall = 0;

    if (i % kStampasim == 0)
      cout << "\n \n  EVENTO  " << i + 1 << endl;
    vertice =
        new Vertex(kMul, multiHisto, kSigmax, kSigmay, kSigmaz); // creo vertice
    nparticle = vertice->Getm();
    if (i % kStampasim == 0)
      cout << "Creato vertice" << endl;
    // loop sul numero di particelle cariche generate
    for (int j = 0; j < nparticle; j++) { // loop su particelle prodotte
      vertice->Initialdir(
          kEtamin, kEtamax,
          etaHisto); // creazione direzione iniziale 1 particella
      Track *tbp =
          new Track(vertice->Gettheta(), vertice->Getphi()); // track verso bp

      // valutazione intersezione con BP
      (void)(tbp->Intersection(
          vertice->Getx(), vertice->Gety(), vertice->Getz(), Xbp, Ybp, Zbp, kRbp, kL));
                    // cast a void (scarta il bool), lunghezza infinita, non
                    // conta se interseca entro +\-L/2
          new (intbp[j]) Intpoint(Xbp, Ybp, Zbp, j);
      tbp->MultipleScattering(kMs);                          // multiple scatter
      Track *t1 = new Track(tbp->Gettheta(), tbp->Getphi()); // traccia verso T1

      // valutazione intersezione con T1 e se esiste proseguo per T2 da T1 con
      // MS
      if (t1->Intersection(Xbp, Ybp, Zbp, X1, Y1, Z1, kR1, kL)) {
        new (int1[j - lost1]) Intpoint(X1, Y1, Z1, j);
        t1->MultipleScattering(kMs);
        Track *t2 = new Track(t1->Gettheta(), t1->Getphi()); // traccia verso T2

        // valutazione intersezione con T2
        if (t2->Intersection(X1, Y1, Z1, X2, Y2, Z2, kR2, kL)) {
          new (int2[j - lost2]) Intpoint(X2, Y2, Z2, j);
        } else {
          lost2++;
        }

        delete t2;
      } else { // valutazione intersezione con T2 per tracce fuori da T1
        if (t1->Intersection(Xbp, Ybp, Zbp, X2, Y2, Z2, kR2, kL)) {
          new (int2[j - lost2]) Intpoint(X2, Y2, Z2, j);
        } else {
          lost2++;
          lostall++;
        }
        lost1++;
      }

      delete t1;
      delete tbp;
    }

    if (i % kStampasim == 0)
      cout << "Valutati hit" << endl;

    // riempimento Tree
    tree->Fill();
    if (lostall == nparticle)
      vertexnohit++;
    clone1->Clear(); // pulizia TClone per evento successivo
    clone2->Clear();
    bp->Clear();
    delete vertice; // eliminazione Vertex per evento dopo
  }
  cout << endl;
  cout << "Su " << kEvents << " vertici generati, " << vertexnohit
       << " non hanno lasciato hit in nessuno dei due tracker" << endl;
  hfile->cd();   // selezione di hfile come directory
  tree->Write(); // scrittura TTree
  delete bp;
  delete clone1;
  delete clone2;
  f->Close();
  hfile->Close(); // chiusura file

  cout << "Simulazione completata" << endl;
}
