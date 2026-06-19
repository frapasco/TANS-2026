# TANS-UniTO
## Abstract

Lo scopo del software è simulare lo stato finale di un evento di scattering tra due fasci di particelle, con produzione di particelle cariche in un collision diamond gaussiano rispetto all'origine e la seguente propagazione nell'apparato sperimentale costituito da una beam pipe e due tracciatori cilindrici.
Nella seconda parte viene implementato un algoritmo di ricostruzione di tipo FAST2 che permette l'analisi dei dati effettuando un vertex finding con cui si ricostruisce il punto di creazione delle particelle secondarie a partire dalle hit delle stesse sui rivelatori.
Si esegue infine un'analisi dei risultati e delle performance di ricostruzione.

## Istruzioni per la compilazione e l'esecuzione del codice

Per l'esecuzione della simulazione e della ricostruzione con conseguente analisi dei risultati bisogna lanciare i seguenti comandi da un terminale con una finestra di root:

```
	.L Compileclass.C
	Compileclass()
	simulazione()
	ricostruzione()
	analysis()
```

N.B. Il codice è stato compilato ed eseguito con:
```
ROOT Version: 6.36.10
Built for macosxarm64
```
