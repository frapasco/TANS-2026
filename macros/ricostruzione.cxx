// root includes
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TStopwatch.h>

// headers include
#include "Intpoint.h"
#include "Track.h"
#include "Vertex.h"

using namespace std;
const int kStamparic = 10000; // ogni quanti eventi fare stampe
const int kLimit1 = 17; // elemento di array di limiti bin istogramma (per
                        // grafici in funzione di molteplicità) a cui fermarsi
                        // per molteplicità da distribuzione assegnata
const int kLimit2 = 21; // analogo per molteplicità uniforme
const double kNoisefrac = 1.2; // numero punti noise/molteplicità
const double kPhimax =
    0.01; // sfasamento massimo fra due punti fra cui creare tracklet (valutato
          // per eccesso fra punti di una stessa traccia da
          // simulazione+ricostruzione indipendenti)
const double kRange = 0.1; // range, attorno al picco dell'istogramma delle
                           // intersezioni delle tracklet, entro cui mediare le
                           // intersezioni per ricostruire z del vertice

const TString kSim = "simulazione.root"; // file da cui leggere eventi
const TString kRic =
    "ricostruzione.root";                // file per salvare tree ricostruzione
const TString kIsto = "istogrammi.root"; // file per salvare istogrammi
const string kData = "data.txt"; // file da cui leggere caratteristiche run

///////////////////////////////////////////////////////////////////
// Funzione per rappresentare i risultati (residui, efficienza vs
// molteplicità, risoluzione vs molteplicità, efficienza vs ztrue,
// risoluzione vs ztrue)
////////////////////////////////////////////////////////////////////

void plot(const vector<double> *const vztruep,
          const vector<double> *const vzrecp, const vector<int> *const vzmultip,
          const vector<double> *const vzrecrmsp, const double *const multibin,
          const int arrayLenghtMulti, TH1D *tot1multi, TH1D *tot3multi,
          const int limit, const double sigmaz, const double *const zbin,
          const int arraylenghtz, TH1D *totz, const int sizetrue);

////////////////////////////////////////////////////////////////////
// Funzione per aggiungere smearing ai punti di intersezione e generare
// noise e ricostruire la posizione del vertice primario a partire
// dalle tracklet.
////////////////////////////////////////////////////////////////////

void ricostruzione() {
  double step = 0.05; // step istogrammi intersezioni tracklet con z
  double R1, R2, L, sigmaz;
  int events;
  char a;
  // letturas di data.txt per recuperare valori di L, R1, R2, sigmaz, events e
  // scelta su come estrarre molteplicità
  ifstream In = ifstream(kData);
  if (!In.is_open()) {
    cout << "Errore, non esiste il file" << endl;
    return;
  }

  In >> R1 >> R2 >> L >> sigmaz >> events >> a;
  In.close();
  TString mul;
  if (a == 'N')
    mul = "No";
  else
    mul = "sì";

  // apertura file di simulazione e creazione tree di lettura
  TFile *hfile = TFile::Open(kSim, "READ");
  if (hfile == NULL) {
    cout << "Errore, non esiste il file" << endl;
    return;
  }
  hfile->ls(); // lettura di contenuto file
  TTree *treeInput = (TTree *)hfile->Get("tree");
  if (treeInput == NULL) {
    cout << "Errore, non esiste il tree" << endl;
    return;
  }

  // creazione file e tree di ricostruzione
  TFile *hfileOutput = TFile::Open(kRic, "RECREATE");
  TTree *treeRec = new TTree("treeRec", "ricostruzione");
  treeRec->SetDirectory(hfileOutput);

  // istogrammi per ricostruzione PV in z
  if (step <= 0 || kStamparic <= 0 || kNoisefrac <= 0 || kPhimax <= 0 ||
      kRange <= 0) {
    cout << "Errore nei parametri forniti " << endl;
    return;
  }

  double zest = R2 * (L / (R2 - R1)) - (L / 2); // z estremo ottenibile da tracklet e centratura nel S.d.R. 
  double zmin = -zest - step / 2;
  double zmaxbin = zest + step / 2;
  int nbin = (int)(zmaxbin - zmin) / step;
  step = (zmaxbin - zmin) /
         nbin; // riadattamento step per coprire tutto il range con nbin

  // creazione istogramma e creazione TClone per punti di intersezione
  TH1D *zintHisto =
      new TH1D("zrec", "intersezioni in z; z [cm]; numero di eventi", nbin,
               zmin, zmaxbin);
  TClonesArray *intpRec1 = new TClonesArray("Intpoint", 200); // T1
  TClonesArray *intpRec2 = new TClonesArray("Intpoint", 200); // T2
  TClonesArray &int1 = *intpRec1;
  TClonesArray &int2 = *intpRec2;

  Vertex *verRec = NULL;        // creazione puntatore a Vertex
  vector<double> Zintersection; // vector per intersezioni Tracklet ricostruite
                                // per un singolo vertice

  // creazione Branches per treeRec
  treeRec->Branch("VertexMCtruth", &verRec); // PV
  treeRec->Branch("zintersectionvector",
                  &Zintersection);          // vector di intersezioni Tracklet
  treeRec->Branch("int1enoise", &intpRec1); // intersezione + noise T1
  treeRec->Branch("int2enoise", &intpRec2); // intersezione + noise T2

  // preparazione lettura treeInput
  TClonesArray *clone1 = new TClonesArray("Intpoint", 100);
  TClonesArray *clone2 = new TClonesArray("Intpoint", 100);
  for (int i = 0; i < 10; i++) {
    new ((*clone1)[i]) Intpoint();
    new ((*clone2)[i]) Intpoint();
  } // per non connettere il branch a oggetto di dimensione 0

  Vertex *verInput = new Vertex(); // creazione oggetto Vertex

  treeInput->GetBranch("tclone1")->SetAutoDelete(kFALSE);
  treeInput->SetBranchAddress("tclone1",
                            &clone1); // indirizzamento per lettura TClone

  treeInput->GetBranch("tclone2")->SetAutoDelete(kFALSE);
  treeInput->SetBranchAddress("tclone2", &clone2);

  TBranch *b3 = treeInput->GetBranch("Vertext");
  b3->SetAddress(&verInput); // lettura Vertex

  clone1->Clear(); // pulizia vettori clone1 e clone2
  clone2->Clear();

  // preparazione istogramma per efficienza(molteplicità) e array limiti bin.
  double multibin[] = {-0.5, 0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,
                       7.5,  9.5,  11.5, 14.5, 17.5, 20.5, 25.5, 30.5,
                       40.5, 50.5, 60.5, 70.5, 80.5, 90.5, 91.5};
  int arrayLenghtMulti = sizeof(multibin) / sizeof(multibin[0]) - 1;
                            // numero bin isto, deve essere 1 in meno di
                            // elementi di vettore per costruttore
  int limit; // dice molteplicità massima isto con multibin[limit]
  if (kLimit1 > kLimit2 || kLimit2 > arrayLenghtMulti || kLimit1 <= 0 ||
      kLimit2 <= 0) {
    cout << "Errore su kLimit1 e kLimit2" << endl;
    return;
  }
  if (mul == "No")
    limit = kLimit2;
  else
    limit = kLimit1; // ci fermiamo prima in molteplicità se non molteplicità
                     // uniforme

  TH1D *tot1multi = new TH1D("tot1multi", "tot1multi", arrayLenghtMulti,
               multibin); // totale eventi con z simulato entro 1 sigma da 0 per
                          // efficienza in funzione di molteplicità
  TH1D *tot3multi = new TH1D("tot3multi", "tot3multi", arrayLenghtMulti,
                             multibin); // entro 3 sigma

  // analogo per efficienza(ztrue=z del vertice primario)
  double zbin[] = {-16, -12, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 12, 16};
  int arraylenghtz = sizeof(zbin) / sizeof(zbin[0]) - 1; 
                        // numero bin isto, deve essere 1 in meno di elementi di
                        // vettore per costruttore
  TH1D *totz = new TH1D("totz", "totz", arraylenghtz,
                        zbin); // totale eventi per efficienza con ztrue

  // vector per fare istogrammi finali
  vector<double> *vztrue = new vector<double>; // z di vertice simulato
  vector<double> *vzrec = new vector<double>;  // z ricostruito
  vector<int> *vzmulti = new vector<int>;      // molteplicità (da simulazione)
  vector<double> *vzrecrms =
      new vector<double>;  // errore (RMS) su z ricostruito
  vztrue->reserve(events); // allocazione di memoria adeguata a evitare
                           // frammentazioni di memoria
  vzrec->reserve(events);
  vzmulti->reserve(events);
  vzrecrms->reserve(events);

  // variabili utili
  int nlines1, nlines2;
  int out1; // contatori per Intpoint fuori da T1 o T2  dopo smearing
  int out2;
  int noTrack = 0; // numero di tracklet per cui non si trova intersezione sensata
  int tottraccia = 0; // totale tracklet generate
  double x1, y1, z1, x2, y2, z2;
  double Z;               // z di intersezione tracklet - piano x=0
  int element1, element2; // numero elementi in int1 e int 2 (senza noise)
  double zrec=0. , zrecrms=0.;   // z ricostruita e suo errore
  double zmax;     // z a cui si trova il primo massimo dell'istogramma delle
                   // intersezioni delle tracklet
  int sizevec = 0; // dimensioni vector allocati su heap
TStopwatch time;

  // loop sul numero di eventi (entrate di treeInput)
  for (int ev = 0; ev < treeInput->GetEntries(); ev++) {
    if (ev % kStamparic == 0)
      cout << "\n \n EVENTO  " << ev + 1 << endl;
    Zintersection.clear(); // pulizia vettore intersezioni
    Zintersection.reserve(200);
    out1 = 0;
    out2 = 0;
    element1 = 0;
    element2 = 0;
    delete verRec;           // eliminazione verRec
    zintHisto->Reset("ICES"); // svuotamento istogramma da evento precedente
    treeInput->GetEntry(ev);   // lettura evento ev
    nlines1 = clone1->GetEntriesFast(); // numero eventi in clone1
    nlines2 = clone2->GetEntriesFast(); // numero eventi in clone2
    verRec = new Vertex(*verInput);       // copia di verInput in verRec

    // lettura punti di intersezione da treeInput e aggiunta dello smearing
    // loop su tracce simulate
    for (int i = 0; i < nlines1; i++) {
      int t = i - out1;

      // estrazione di un' intersezione per volta da T1
      Intpoint *i1 = (Intpoint *)clone1->At(i);
      i1->Smearing(R1);

      // riempimento treeRec int 1 con le coordinate con smearing
      if (TMath::Abs(i1->Getz()) < L / 2) { // il punto rimane dentro T1
        x1 = i1->Getx();
        y1 = i1->Gety();
        z1 = i1->Getz();
        new (int1[t]) Intpoint(x1, y1, z1, i1->Getlabel());
        element1++;
      } else
        out1++;
    }
    for (int i = 0; i < nlines2; i++) {
      int k = i - out2;

      // estrazione di un'intersezione per volta da T2

      Intpoint *i2 = (Intpoint *)clone2->At(i);
      i2->Smearing(R2);

      // riempimento treeRec int 2 con le coordinate con smearing
      if (TMath::Abs(i2->Getz()) < L / 2) { // il punto rimane dentro T2
        x2 = i2->Getx();
        y2 = i2->Gety();
        z2 = i2->Getz();
        new (int2[k]) Intpoint(x2, y2, z2, i2->Getlabel());
        element2++;
      } else
        out2++;
    }
    if (ev % kStamparic == 0)
      cout << "Applicato smearing" << endl;
    if (ev % kStamparic == 0)
      cout << "Hit usciti da T1: " << out1 << "  Hit usciti da T2: " << out2 << endl;

    // aggiunta del noise
    int Nnoise = (int)(verInput->Getm()) * kNoisefrac; // numero di punti di noise, assunto identico per T1 e T2
    for (int j = element1; j < element1 + Nnoise; j++) {
      new (int1[j]) Intpoint(R1, L);
    }
    for (int k = element2; k < element2 + Nnoise; k++) {
      new (int2[k]) Intpoint(R2, L);
    }
    if (ev % kStamparic == 0)
      cout << "Applicato noise" << endl;

    // inzio ricostruzione PV
    int Nrec1 = intpRec1->GetEntriesFast(); // numero di elementi (hit+noise) in int1 e int2
    int Nrec2 = intpRec2->GetEntriesFast();

    int ntraccia = 0; // contatore tracklet fatte per un vertice
    for (int i = 0; i < Nrec1; i++) {
      for (int j = 0; j < Nrec2; j++) {
        Intpoint *i1 = (Intpoint *)intpRec1->At(i);
        Intpoint *i2 = (Intpoint *)intpRec2->At(j);
        if (TMath::Abs(i2->Getphi() - i1->Getphi()) < kPhimax) { // creazione tracklet per punti compatibili entro phimax
          x1 = i1->Getx();
          x2 = i2->Getx();
          y1 = i1->Gety();
          y2 = i2->Gety();
          z1 = i1->Getz();
          z2 = i2->Getz();
          Track *t = new Track(x2 - x1, y2 - y1, z2 - z1);
          ntraccia++;
          tottraccia++;
          Z = 200000000; // per controllo di funzionamento funzione intersezione
          t->Intersection0(x2, y2, z2, Z); // aggiorna Z con intersezione tracklet - piano x=0
          if (Z > -10000000 && Z < 10000000) { // si riempie vettore e istogramma con
                              // intersezioni riuscite
            zintHisto->Fill(Z);
            Zintersection.push_back(Z);
          } else
            noTrack++;
          delete t;
        } // fine if su phi1-phi2<kPhimax
      }
    } // fine for sui due TClone di intersezioni per fare tracklet

    if (ev % kStamparic == 0)
      cout << "Create tracklet e riempito istogramma" << endl;
    if (zintHisto->GetEntries() != 0) { // if per chiedere istogramma non vuoto
      int binmax = zintHisto->GetMaximumBin(); // numero del bin col primo
                                              // massimo dell'istogramma
      double zmax =
          (binmax - 1) * step + zmin + step / 2; // zmin è limite inferiore di range dell'ascissa di
                    // istogramma; centro di bin (in z) dove si ha primo massimo

                 
      // calcolo zrec come media dei valori di intersezione (aggiunti
      // precedentemente in vector) in un range attorno a zmax
      double sum = 0.;
      double scarti = 0; // somma di scarti quadratici sui valori mediati per calcolo rms su zrec
      int denom = 0;
      int size = Zintersection.size();
      for (int i = 0; i < size; i++) {
        if (TMath::Abs(Zintersection.at(i) - zmax) < kRange) { // scelta range da osservazione degli istogrammi
          sum += Zintersection.at(i);
          scarti = scarti + ((Zintersection.at(i) - zrec) *
                               (Zintersection.at(i) - zrec));
          denom++;
        }
      }
      if (denom > 0) {
        zrec = sum / denom;
        if (denom > 1)
          zrecrms =
              TMath::Sqrt(scarti / (denom - 1)); // calcolo deviazione standard
        else {
          zrecrms = 0.05;
        } // 0.05 scelto come errore plausibile nel caso di una media calcolata
          // su un singolo valore
        
        if (ev % kStamparic == 0)
          cout << "zrec= " << zrec << " +- " << zrecrms << endl;
        vzrec->push_back(zrec); // riempimento vector
        vztrue->push_back(verInput->Getz());
        vzmulti->push_back(verInput->Getm());
        vzrecrms->push_back(zrecrms);
        sizevec++; // contatore per tenere conto delle dimensioni dei vector
      } // fine if denom>0
    } // fine richiesta istogramma non vuoto
      
    if (ev % kStamparic == 0) {
      new TCanvas;
      zintHisto->DrawCopy();
    }
    if (verInput->Getz() < 1 * sigmaz && verInput->Getm() < multibin[limit])
      tot1multi->Fill(verInput->Getm());
    if (verInput->Getz() < 3 * sigmaz && verInput->Getm() < multibin[limit])
      tot3multi->Fill(verInput->Getm());
    totz->Fill(verInput->Getz()); // riempimento istogrammi con eventi totali
                                // simulati (efficienza=buoni/totali)

    // riempimento tree e pulizia TClone
    treeRec->Fill();
    clone1->Clear();
    clone2->Clear();
    intpRec1->Clear();
    intpRec2->Clear();

  } // fine loop su eventi

  cout << "\n \n Zintersection non trovata correttamente per " << noTrack
       << " Tracklet su " << tottraccia << " Trackelet costruite" << endl;

  plot(vztrue, vzrec, vzmulti, vzrecrms, multibin, arrayLenghtMulti, tot1multi,
       tot3multi, limit, sigmaz, zbin, arraylenghtz, totz, sizevec);

  delete vztrue; // eliminazione vector
  delete vzrec;
  delete vzmulti;
  delete vzrecrms;

  hfileOutput->cd();
  treeRec->Write(); // scrittura tree su file e chiusura file

  delete intpRec1;
  delete intpRec2;
  delete clone1;
  delete clone2;
  delete verInput;
  delete verRec;

  hfileOutput->Close();
  hfile->Close();

  cout << "Ricostruzione completata" << endl;
  time.Stop();
  time.Print();
}

void plot(const vector<double> *const vztruep,
          const vector<double> *const vzrecp, const vector<int> *const vzmultip,
          const vector<double> *const vzrecrmsp, const double *const multibin,
          const int arrayLenghtMulti, TH1D *tot1multi, TH1D *tot3multi,
          const int limit, const double sigmaz, const double *const zbin,
          const int arraylenghtz, TH1D *totz, const int sizetrue) {

  TFile *histo = TFile::Open(kIsto, "RECREATE"); // file per salvare istogrammi
  const vector<double> &vztrue = *vztruep;
  const vector<double> &vzrec = *vzrecp;
  const vector<int> &vzmulti = *vzmultip;
  const vector<double> &vzrecrms = *vzrecrmsp;

  char nome[200]; // nome e titolo istogrammi
  char titolo[200];

  // istogrammi dei residui
  TH1D *isto_res[3];
  int min_mult1 = 3; // range molteplicità per istogrammi
  int max_mult1 = 5;
  int min_mult2 = 45;
  int max_mult2 = 55;
  // residui con tutte le molteplicità
  snprintf(nome, sizeof(nome), "isto_restot");
  snprintf(titolo, sizeof(titolo),  "Residui con tutte le molteplicità");
  isto_res[0] = new TH1D(nome, titolo, 500, -0.1, 0.1);
  // residui con molteplicità tra 3 e 5
  snprintf(nome, sizeof(nome), "histo_res1");
  snprintf(titolo, sizeof(titolo),  "%i <= molteplicita' <= %i", min_mult1, max_mult1);
  isto_res[1] = new TH1D(nome, titolo, 500, -0.1, 0.1);
  // residui con molteplicità tra 45 e 55
  snprintf(nome, sizeof(nome), "histo_res2");
  snprintf(titolo, sizeof(titolo),  "%i <= molteplicita' <= %i", min_mult2, max_mult2);
  isto_res[2] = new TH1D(nome, titolo, 500, -0.1, 0.1);

  for (int i = 0; i < sizetrue; i++) {
    int multi = vzmulti[i];
    if (TMath::Abs(vztrue[i] - vzrec[i]) < 3 * vzrecrms[i]) { // così riempiamo solo con gli eventi ricostruiti
                           // discretamente e non con quelli con residui enormi
                           // dovuti al prendere il picco sbagliato
                           // dell'istogramma
      isto_res[0]->Fill(vzrec[i] - vztrue[i]); //istogramma con residui di tutte le molteplicità
      if (multi >= min_mult1 && multi <= max_mult1)
        isto_res[1]->Fill(vzrec[i] - vztrue[i]);
      if (multi >= min_mult2 && multi <= max_mult2)
        isto_res[2]->Fill(vzrec[i] - vztrue[i]);
    }
  }

  // istogrammi eventi ben ricostruiti per efficienza vs molteplicità
  TH1D *isto_ok[2];
  snprintf(nome, sizeof(nome), "isto_ok1sigma");
  snprintf(titolo, sizeof(titolo),  "ok 1 sigma");
  isto_ok[0] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);
  snprintf(nome, sizeof(nome), "isto_ok3sigma");
  snprintf(titolo, sizeof(titolo),  "ok 3 sigma");
  isto_ok[1] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);

  for (int i = 0; i < sizetrue; i++) {
    if (TMath::Abs(vztrue[i] - vzrec[i]) < 3 * vzrecrms[i] &&
        vztrue[i] < 1 * sigmaz && vzmulti[i] < multibin[limit])
      isto_ok[0]->Fill(vzmulti[i]);
    if (TMath::Abs(vztrue[i] - vzrec[i]) < 3 * vzrecrms[i] &&
        vztrue[i] < 3 * sigmaz && vzmulti[i] < multibin[limit])
      isto_ok[1]->Fill(vzmulti[i]);
  }

  // istogrammi efficienza vs molteplicità
  TH1D *isto_eff[2];
  // z fra [-sigma;sigma]
  snprintf(nome, sizeof(nome), "isto_eff1sigma");
  snprintf(titolo, sizeof(titolo),  "efficienza 1 sigma");
  isto_eff[0] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);
  // z fra [-3 sigma;3 sigma]
  snprintf(nome, sizeof(nome), "isto_eff3sigma");
  snprintf(titolo, sizeof(titolo),  "efficienza 3 sigma");
  isto_eff[1] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);

  for (int i = 1; i < arrayLenghtMulti + 1; i++) {
    double e;
    double error;
    if (tot1multi->GetBinContent(i) != 0) {
      e = isto_ok[0]->GetBinContent(i) /
          tot1multi->GetBinContent(i); // efficienza =buoni/totali
      error = TMath::Sqrt(e * (1 - e) / tot1multi->GetBinContent(i));
      if (1 / tot1multi->GetBinContent(i) > error)
        error = 1 / tot1multi->GetBinContent(
                        i); // per evitare errore =0 per efficienza e=1 o e=0
      isto_eff[0]->SetBinContent(i, e);
      isto_eff[0]->SetBinError(i, error);
    }
    if (tot3multi->GetBinContent(i) != 0) {
      e = isto_ok[1]->GetBinContent(i) / tot3multi->GetBinContent(i);
      error = TMath::Sqrt(e * (1 - e) / tot3multi->GetBinContent(i));
      if (1 / tot3multi->GetBinContent(i) > error)
        error = 1 / tot3multi->GetBinContent(i);
      isto_eff[1]->SetBinContent(i, e);
      isto_eff[1]->SetBinError(i, error);
    }
  }

  // istogrammi residui bin per bin per risoluzione vs. molteplicità per ztrue
  // entro 1 sigma da 0
  vector<TH1D*> isto_residui1(arrayLenghtMulti);
  // analogo per 3 sigma
  vector<TH1D*> isto_residui3(arrayLenghtMulti);
  for (int i = 0; i < arrayLenghtMulti; i++) {
    snprintf(nome, sizeof(nome), "isto_residui1[%i]", i);
    snprintf(titolo, sizeof(titolo),  "residui 1 sigma per %f<molteplicita'<%f", multibin[i],
            multibin[i + 1]);
    isto_residui1[i] = new TH1D(nome, titolo, 500, -0.1, 0.1);
    snprintf(nome, sizeof(nome), "isto_residui3[%i]", i);
    snprintf(titolo, sizeof(titolo),  "residui 3 sigma per %f<moteplicita'<%f", multibin[i],
            multibin[i + 1]);
    isto_residui3[i] = new TH1D(nome, titolo, 500, -0.1, 0.1);

    for (int j = 0; j < sizetrue; j++) {
      if (TMath::Abs(vztrue[j] - vzrec[j]) < 3 * vzrecrms[j] &&
          vztrue[j] < 1 * sigmaz && vzmulti[j] < multibin[limit] &&
          vzmulti[j] > multibin[i] && vzmulti[j] <= multibin[i + 1])
        isto_residui1[i]->Fill(vzrec[j] - vztrue[j]);
      if (TMath::Abs(vztrue[j] - vzrec[j]) < 3 * vzrecrms[j] &&
          vztrue[j] < 3 * sigmaz && vzmulti[j] < multibin[limit] &&
          vzmulti[j] > multibin[i] && vzmulti[j] <= multibin[i + 1])
        isto_residui3[i]->Fill(vzrec[j] - vztrue[j]);
    }
  }

  // istogrammi risoluzione vs molteplicità 1 sigma e 3 sigma
  TH1D *isto_ris[2];
  // z fra [-sigma;sigma]
  snprintf(nome, sizeof(nome), "isto_ris1sigma");
  snprintf(titolo, sizeof(titolo),  "risoluzione 1 sigma");
  isto_ris[0] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);
  // z fra [-3 sigma;3 sigma]
  snprintf(nome, sizeof(nome), "isto_ris3sigma");
  snprintf(titolo, sizeof(titolo),  "risoluzione 3 sigma");
  isto_ris[1] = new TH1D(nome, titolo, arrayLenghtMulti, multibin);

  for (int i = 1; i < arrayLenghtMulti + 1; i++) {
    double RMS; // risoluzione è RMS istogrammi residui bin per bin
    double error;
    if (isto_residui1[i - 1]->GetRMS() != 0) {
      RMS = isto_residui1[i - 1]->GetRMS();
      error = isto_residui1[i - 1]->GetRMSError();

      isto_ris[0]->SetBinContent(i, RMS);
      isto_ris[0]->SetBinError(i, error);
    }
    if (isto_residui3[i - 1]->GetRMS() != 0) {
      RMS = isto_residui3[i - 1]->GetRMS();
      error = isto_residui3[i - 1]->GetRMSError();

      isto_ris[1]->SetBinContent(i, RMS);
      isto_ris[1]->SetBinError(i, error);
    }
  }

  // analogo per istogrammi in funzione ztrue

  // istogrammi eventi ben ricostruiti per efficienza vs ztrue
  TH1D *isto_okztrue[1];
  snprintf(nome, sizeof(nome), "isto_okztrue");
  snprintf(titolo, sizeof(titolo),  "ok ztrue");
  isto_okztrue[0] = new TH1D(nome, titolo, arraylenghtz, zbin);

  for (int i = 0; i < sizetrue; i++) {
    if (TMath::Abs(vztrue[i] - vzrec[i]) < 3 * vzrecrms[i])
      isto_okztrue[0]->Fill(vztrue[i]);
  }

  // istogrammi efficienza vs ztrue
  TH1D *isto_effztrue[1];
  snprintf(nome, sizeof(nome), "isto_effztrue");
  snprintf(titolo, sizeof(titolo),  "efficienza vs ztrue ");
  isto_effztrue[0] = new TH1D(nome, titolo, arraylenghtz, zbin);

  for (int i = 1; i < arraylenghtz + 1; i++) {
    double e;
    double error;
    if (totz->GetBinContent(i) != 0) {
      e = isto_okztrue[0]->GetBinContent(i) / totz->GetBinContent(i);
      error = TMath::Sqrt(e * (1 - e) / totz->GetBinContent(i));
      if (1 / totz->GetBinContent(i) > error)
        error = 1 / totz->GetBinContent(i);
      isto_effztrue[0]->SetBinContent(i, e);
      isto_effztrue[0]->SetBinError(i, error);
    }
  }

  // istogramma per verificare distribuzione residui complessiva su tutti gli
  // ztrue
  TH1D *isto_ressomma = new TH1D("isto_ressomma",
                                 "somma di istogrammi residui bin per bin di "
                                 "ztrue; z_rec-z_true [cm]; numero di eventi",
                                 500, -0.1, 0.1);

  // istogrammi residui bin per bin per risoluzione vs ztrue
  vector<TH1D*> isto_residuiztrue(arraylenghtz);

  for (int i = 0; i < arraylenghtz; i++) {
    snprintf(nome, sizeof(nome), "isto_residuiztrue[%i]", i);
    snprintf(titolo, sizeof(titolo),  "residui  per %f<ztrue <=%f", zbin[i], zbin[i + 1]);
    isto_residuiztrue[i] = new TH1D(nome, titolo, 500, -0.1, 0.1);
    for (int j = 0; j < sizetrue; j++) {
      if (TMath::Abs(vztrue[j] - vzrec[j]) < 3 * vzrecrms[j] &&
          vztrue[j] > zbin[i] && vztrue[j] <= zbin[i + 1]) {
        isto_residuiztrue[i]->Fill(vzrec[j] - vztrue[j]);
      }
    }
    isto_ressomma->Add(isto_residuiztrue[i]);
  }

  // istogrammi risoluzione vs ztrue
  TH1D *isto_risztrue[1];
  snprintf(nome, sizeof(nome), "isto_risztrue");
  snprintf(titolo, sizeof(titolo),  "risoluzione vs ztrue ");
  isto_risztrue[0] = new TH1D(nome, titolo, arraylenghtz, zbin);

  for (int i = 1; i < arraylenghtz + 1; i++) {
    double RMS;
    double error;
    if (isto_residuiztrue[i - 1]->GetRMS() != 0) {
      RMS = isto_residuiztrue[i - 1]->GetRMS();
      error = isto_residuiztrue[i - 1]->GetRMSError();
      isto_risztrue[0]->SetBinContent(i, RMS);
      isto_risztrue[0]->SetBinError(i, error);
    }
  }

  histo->Write(); // scrittura istogrammi su file
  histo->ls();    // controllo contenuto file
  histo->Close();
}
