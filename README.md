IDENTIFICAZIONE SISTEMI INCERTI - DOPPIO PENDOLO INVERSO
STUDENTE: PAOLO  

In questa cartella sono presenti vari script Matlab, per eseguire il progetto è sufficente aprire solo quello denominato come
"Progetto_ISI_finale", tutti gli altri verranno eseguiti automaticamente all'interno di questo, anche il Simulink si aprirà da solo.

° "Progetto_ISI_finale" è diviso in 9 sezioni per lo più indipendenti:

	1. "Definizione dei parametri fisici" --> è la prima da avviare (apre anche lo schema Simulink da avviare manualmente eventualmente 
		abilitando gli outliers di misura) 

	2. "Estrazioni variabili da Simulink" --> va eseguito per estrarre tutte le ulteriori variabili necessarie per l'esecuzione dei
		 filtri successivi

(Queste due sezioni vanno SEMPRE rieseguite ogni volta che si vuole eseguire una nuova sezione da 4 a 8 per evitare sovrapposizioni
 di variabili nel workspace)

	3. "Simulazione doppio pendolo" --> non è necessario eseguirla e non necessita nel caso di rieseguire la 1. e 2., serve per vedere 
		una animazione di come si muove il pendolo

	4. "EKF" --> applica il filtro di Kalman esteso per stimare lo stato, vengono subuto graficati i risultati ed è chiesto se si voglia
		vedere una animazione del pendolo stimato sovrapposta ad una animazione analoga a quella proposta dalla sezione 3.

	5. "UKF" --> (ricorda di rieseguire 1. e 2. prima di avviarla) effettua la stima dinamica tramite un Unscented Kalman Filter

	6. "Smoother"--> (ricorda di rieseguire 1. e 2. prima di avviarla) esegue una stima in avanti con EKF e ricorsione in indietro per 
		risolvere un problema di regolarizzazione

	7. "EKF con diversi tempi di campionamento" --> (ricorda di rieseguire 1., modificare almeno un tempo di campionamento di un sensore
		 e rieseguire 2. prima di avviarla) risolve il problema afforontato nella sezione 4 ma considerando la possibilità di avere 
		sensori non sincroni tra loro nè con filtro

	8. "EKF con outliers" --> (ricorda di rieseguire 1. e 2. prima di avviarla, nel Simulink abilita gli outliers con i manual switch)
		risolve il problema affrontato nella sezione 4 ma con la possibilità di eventuali outliers di misura da reiettare

	9. "Index funzioni" --> raccoglie tutte le funzoni usate nello script

° Tutti gli altri script presenti servono per effettuare i vari grafici e simulazioni per i vari filtri EKF, UKF e Smoother. 

° "Matrici linearizzato" è servito per trovare la espressione analitica dei jacobiani ma è obsoleto.

° "EKF_outliers" può essere pure ignorato in quanto è stato il banco di prova per la sezione 8
