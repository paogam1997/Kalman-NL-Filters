%% Subplot_confronto
%Questo script è usato per l'EKF per graficare i risultati importanti
%ottenuti ed avviare la eventuale simulazione di confronto

%===== Grafici andamento reale vs stimato =====%

figure; %apre nuova figura vuota
%Nella stessa figure inserisco 4 grafici quindi sono specificati i subplot
%con relativa posizione dove inserire ogni volta il nuovo grafico
subplot (2, 2, 1)                                                               %crea 4 sezioni e metti il grafico nella prima
plot(time,teta1, 'b-'); grid on;xlabel('tempo [s]');ylabel('\theta_1 [rad]')    %grafica andamento di teta1 [rad] nel tempo [s]
hold on; plot(time1,log_results.x_hat(:,1), 'r--')                              %nello stesso grafico, inserisci anche l'andamento stimato della teta1
title('\theta_1'); legend('vero', 'stima');                                     %Inserimento del titolo e legenda per il grafico fatto

subplot (2, 2, 2)                                                               %inserisci il grafico sottostante nella seconda sezione
plot(time,teta2, 'b-'); grid on;xlabel('tempo [s]');ylabel('\theta_2 [rad]')    %grafica andamento di teta2 [rad] nel tempo [s]
hold on; plot(time1,log_results.x_hat(:,2), 'r--')                              %su stesso grafico traccia andamento stimato della teta2
title('\theta_2'); legend('vero', 'stima');                                     %Inserimento del titolo e legenda

subplot (2, 2, 3)                                                               %inserisci il grafico sottostante nella terza sezione
plot(time,dteta1, 'b-'); grid on;xlabel('tempo [s]');ylabel('\omega_1 [rad/s]') %grafica andamento di dteta1 [rad/s] nel tempo [s]
hold on; plot(time1,log_results.x_hat(:,3), 'r--')                              %su stesso grafico traccia andamento stimato della dteta1
title('\omega_1'); legend('vero', 'stima');                                     %Inserimento del titolo e legenda

subplot (2, 2, 4)                                                               %inserisci il grafico sottostante nella quarta sezione
plot(time,dteta2, 'b-'); grid on;xlabel('tempo [s]');ylabel('\omega_2 [rad/s]') %grafica andamento di dteta2 [rad/s] nel tempo [s]
hold on; plot(time1,log_results.x_hat(:,4), 'r--')                              %su stesso grafico traccia andamento stimato della dteta2
title('\omega_2'); legend('vero', 'stima');                                     %Inserimento del titolo e legenda

%===== Andamnento dell'errore di stima (reale-stimato) =====%

figure; %Prendi una nuova figura vuota
subplot (2, 2, 1)                                                                                              %metti il grafico nella prima sezione
plot(time1,teta1(2:end)-log_results.x_hat(:,1)); grid on;xlabel('tempo [s]');ylabel('\theta_1 tilde [rad]')    %grafico errore commesso nella stima di teta1 (qua si vede bene effetto di distubo di integrazione)
title('errore di stima su \theta_1')                                                                           %Inserimento del titolo

subplot (2, 2, 2)                                                                                              %metti il grafico nella seconda sezione
plot(time1,teta2(2:end)-log_results.x_hat(:,2) ); grid on;xlabel('tempo [s]');ylabel('\theta_2 tilde [rad]')   %grafico errore commesso nella stima di teta2 (qua si vede bene effetto di distubo di integrazione)
title('errore di stima su \theta_2')                                                                           %Inserimento del titolo
 
subplot (2, 2, 3)                                                                                              %metti il grafico nella terza sezione
plot(time1,dteta1(2:end)-log_results.x_hat(:,3)); grid on;xlabel('tempo [s]');ylabel('\omega_1 tilde [rad/s]') %grafico errore commesso nella stima di dteta1
title('errore di stima su \omega_1')                                                                           %Inserimento del titolo

subplot (2, 2, 4);                                                                                             %metti il grafico nella quarta sezione
plot(time1,dteta2(2:end)-log_results.x_hat(:,4)); grid on;xlabel('tempo [s]');ylabel('\omega_2 tilde [rad/s]') %grafico errore commesso nella stima di dteta2
title('errore di stima su \omega_2')                                                                           %Inserimento del titolo

%===== Grafico innovazione per controllare sia un rumore bianco =====%

figure;                                                                    %nuova figura
subplot(2,1,1)                                                             %inserisci grafico su prima delle due sezioni
plot(time1,log_results.e(:,1));xlabel('tempo [s]');ylabel('e(1) [m]')      %grafica andamento di prima componente della innovazione nel tempo
title('innovazione y_d');                                                  %inserimento titolo
subplot(2,1,2)                                                             %inserisci grafico nella seconda sezione (di sotto)
autocorr(squeeze(log_results.e(:,1,:)));                                   %per verificare la bianchezza di e se ne grafica l'autocorrelazione (impulso in 0 per un esatto wn)

figure;                                                                    %nuova figura
subplot(2,1,1)                                                             %inserisci grafico su prima delle due sezioni
plot(time1, log_results.e(:,2));xlabel('tempo [s]');ylabel('e(2) [rad]')   %grafica andamento di seconda componente della innovazione nel tempo
title('innovazione y_ alfa');                                              %inserimento titolo
subplot(2,1,2)                                                             %inserisci grafico nella seconda sezione (di sotto)
autocorr(squeeze(log_results.e(:,2,:)));                                   %per verificare la bianchezza di e se ne grafica l'autocorrelazione (impulso in 0 per un esatto wn)

%===== Richiesta per avviare la simulazione di confronto =====%
prompt="Avviare simulazione di confronto Reale vs Stima? Y/N [Y]:";
txt=input(prompt,"s");
if isempty(txt)
    txt='Y';                          %se premi direttamente invio è visto come un consenso
end

if txt=='Y'
    disp('Ok simulo')
    run("Simulazione_RealevsStima.m") %se desiderato esegue lo script per graficare sovrapposti il sistema reale e stimato
else
    disp('Ok niente simulazione')
end

%=== Altro ===%
 %plot(time1,squeeze(log_results.tr(1,:,:)  %valore singolare massimo della
 %matrice di guadagno di correzione