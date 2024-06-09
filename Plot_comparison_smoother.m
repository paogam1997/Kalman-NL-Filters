%% Plot_comparison_smoother
%Similmente allo script "Subplot_confronto" vado a graficare i risultati di
%stimatore e regolarizzatore per vedere la bontà del codice proposto

%Creo una nuova figura vuota
figure; hold on; grid on; box on     


%===== Grafici andamento reale vs stimato vs regolarizzato =====%

subplot (2, 2, 1);                                                              %metti il grafico nella prima sezione
plot(time,teta1, 'b-'); grid on;xlabel('tempo [s]');ylabel('\theta_1 [rad]')    %grafica andamento di teta1 [rad] nel tempo [s]
hold on; plot(time1,squeeze(log_EKF_s.x_hat_correction(1,1,:)), 'r--')          %nello stesso grafico, inserisci anche l'andamento stimato della teta1
plot(time1,squeeze(log_EKF_s.x_hat_smoothed(1,1,:)), 'g -.')                    %ed anche l'andamento derivante dalla regolarizzazione
title('\theta_1')                                                               %Inserimento del titolo

subplot (2, 2, 2);                                                              %inserisci il grafico sottostante nella seconda sezione
plot(time,teta2, 'b-');grid on;xlabel('tempo [s]');ylabel('\theta_2 [rad]')     %grafica andamento di teta2 [rad] nel tempo [s]
hold on; plot(time1,squeeze(log_EKF_s.x_hat_correction(2,1,:)), 'r--')          %nello stesso grafico, inserisci anche l'andamento stimato della teta2
plot(time1,squeeze(log_EKF_s.x_hat_smoothed(2,1,:)), 'g-.')                     %ed anche l'andamento derivante dalla regolarizzazione
title('\theta_2')                                                               %Inserimento del titolo

subplot (2, 2, 3);                                                              %inserisci il grafico sottostante nella terza sezione
plot(time,dteta1, 'b');grid on;xlabel('tempo [s]');ylabel('\omega_1 [rad/s]')   %grafica andamento di dteta1 [rad/s] nel tempo [s]
hold on; plot(time1,squeeze(log_EKF_s.x_hat_correction(3,1,:)), 'r--')          %nello stesso grafico, inserisci anche l'andamento stimato della dteta1
plot(time1,squeeze(log_EKF_s.x_hat_smoothed(3,1,:)), 'g-.')                     %ed anche l'andamento derivante dalla regolarizzazione
title('\omega_1')                                                               %Inserimento del titolo

subplot (2, 2, 4);                                                              %inserisci il grafico sottostante nella quarta sezione
plot(time,dteta2, 'b-'); grid on;xlabel('tempo [s]');ylabel('\omega_2 [rad/s]') %grafica andamento di dteta2 [rad/s] nel tempo [s]
hold on; plot(time1,squeeze(log_EKF_s.x_hat_correction(4,1,:)), 'r--')          %nello stesso grafico, inserisci anche l'andamento stimato della dteta2
plot(time1,squeeze(log_EKF_s.x_hat_smoothed(4,1,:)), 'g-.')                     %ed anche l'andamento derivante dalla regolarizzazione
title('\omega_2')                                                               %Inserimento del titolo

legend('reale','EKF','SMOOTHER')                                                %inserisco una unica legenda per semplicità

%===== Andamnento dell'errore di stima (reale-stimato) =====%

figure;                                                                                                                                %Creo una nuova figura vuota
subplot (2, 2, 1)                                                                                                                      %metti il grafico nella prima sezione
plot(time1,teta1(2:end)-squeeze(log_EKF_s.x_hat_correction(1,1,:)),'r'); grid on;xlabel('tempo [s]');ylabel('\theta_1 tilde [rad]')    %grafico errore commesso nella stima di teta1 (qua si vede bene effetto di distubo di integrazione)
title('errore di stima su \theta_1')                                                                                                   %Inserimento del titolo

subplot (2, 2, 2)                                                                                                                      %inserisci il grafico sottostante nella seconda sezione
plot(time1,teta2(2:end)-squeeze(log_EKF_s.x_hat_correction(2,1,:)),'r' ); grid on;xlabel('tempo [s]');ylabel('\theta_2 tilde [rad]')   %grafico errore commesso nella stima di teta2 (qua si vede bene effetto di distubo di integrazione)
title('errore di stima su\theta_2')                                                                                                    %Inserimento del titolo
 
subplot (2, 2, 3)                                                                                                                      %inserisci il grafico sottostante nella terza sezione
plot(time1,dteta1(2:end)-squeeze(log_EKF_s.x_hat_correction(3,1,:)),'r'); grid on;xlabel('tempo [s]');ylabel('\omega_1 tilde [rad/s]') %grafico errore commesso nella stima di dteta1
title('errore di stima su \omega_1')                                                                                                   %Inserimento del titolo

subplot (2, 2, 4);                                                                                                                     %inserisci il grafico sottostante nella quarta sezione
plot(time1,dteta2(2:end)-squeeze(log_EKF_s.x_hat_correction(4,1,:)),'r'); grid on;xlabel('tempo [s]');ylabel('\omega_2 tilde [rad/s]') %grafico errore commesso nella stima di dteta2
title('errore di stima su \omega_2')                                                                                                   %Inserimento del titolo

%===== Andamnento differenza stimato - regolarizzato =====%

figure;                                                                                                                              %Creo una nuova figura vuota
subplot (2, 2, 1)                                                                                                                    %metti il grafico nella prima sezione
plot(time1,squeeze(log_EKF_s.x_hat_correction(1,1,:))-squeeze(log_EKF_s.x_hat_smoothed(1,1,:)),'k'); grid on;xlabel('tempo [s]');ylabel('err_\theta_1  [rad]')    %grafico errore tra stima e regolarizzazione di teta1
title('differenza stima - regolarizzazione su \theta_1')                                                                             %Inserimento del titolo

subplot (2, 2, 2)                                                                                                                    %inserisci il grafico sottostante nella seconda sezione
plot(time1,squeeze(log_EKF_s.x_hat_correction(2,1,:))-squeeze(log_EKF_s.x_hat_smoothed(2,1,:)),'k' ); grid on;xlabel('tempo [s]');ylabel('err_\theta_2  [rad]')   %grafico errore tra stima e regolarizzazione di teta2
title('differenza stima - regolarizzazione su \theta_2')                                                                              %Inserimento del titolo
 
subplot (2, 2, 3)                                                                                                                    %inserisci il grafico sottostante nella terza sezione
plot(time1,squeeze(log_EKF_s.x_hat_correction(3,1,:))-squeeze(log_EKF_s.x_hat_smoothed(3,1,:)),'k'); grid on;xlabel('tempo [s]');ylabel('err_\omega_1  [rad/s]') %grafico errore tra stima e regolarizzazione di dteta1
title('differenza stima - regolarizzazione su \omega_1')                                                                             %Inserimento del titolo

subplot (2, 2, 4);                                                                                                                   %inserisci il grafico sottostante nella quarta sezione
plot(time1,squeeze(log_EKF_s.x_hat_correction(4,1,:))-squeeze(log_EKF_s.x_hat_smoothed(4,1,:)),'k'); grid on;xlabel('tempo [s]');ylabel('err_\omega_2  [rad/s]') %grafico errore tra stima e regolarizzazione di dteta2
title('differenza stima - regolarizzazione su \omega_2')                                                                             %Inserimento del titolo

%===== Andamnento dell'errore di regolarizzazione (reale-regolarizzato) =====%

figure;                                                                                                                              %Creo una nuova figura vuota
subplot (2, 2, 1)                                                                                                                    %metti il grafico nella prima sezione
plot(time1,teta1(2:end)-squeeze(log_EKF_s.x_hat_smoothed(1,1,:)),'g'); grid on;xlabel('tempo [s]');ylabel('\theta_1 tilde [rad]')    %grafico errore commesso nella regolarizzazione di teta1
title('errore di regolarizzazione su \theta_1')                                                                                      %Inserimento del titolo

subplot (2, 2, 2)                                                                                                                    %inserisci il grafico sottostante nella seconda sezione
plot(time1,teta2(2:end)-squeeze(log_EKF_s.x_hat_smoothed(2,1,:)),'g' ); grid on;xlabel('tempo [s]');ylabel('\theta_2 tilde [rad]')   %grafico errore commesso nella regolarizzazione di teta2
title('errore di regolarizzazione su \theta_2')                                                                                      %Inserimento del titolo
 
subplot (2, 2, 3)                                                                                                                    %inserisci il grafico sottostante nella terza sezione
plot(time1,dteta1(2:end)-squeeze(log_EKF_s.x_hat_smoothed(3,1,:)),'g'); grid on;xlabel('tempo [s]');ylabel('\omega_1 tilde [rad/s]') %grafico errore commesso nella regolarizzazione di dteta1
title('errore di regolarizzazione su \omega_1')                                                                                      %Inserimento del titolo

subplot (2, 2, 4);                                                                                                                   %inserisci il grafico sottostante nella quarta sezione
plot(time1,dteta2(2:end)-squeeze(log_EKF_s.x_hat_smoothed(4,1,:)),'g'); grid on;xlabel('tempo [s]');ylabel('\omega_1 tilde [rad/s]') %grafico errore commesso nella regolarizzazione di dteta2
title('errore di regolarizzazione su \omega_2')                                                                                      %Inserimento del titolo

%===== Grafico innovazione per controllare sia un rumore bianco =====%

figure;                                                                    %Creo una nuova figura vuota
subplot(2,1,1)                                                             %inserisci grafico su prima delle due sezioni
plot(time1,log_EKF_s.e(:,1));xlabel('tempo [s]');ylabel('e(1) [m]')        %grafica andamento di prima componente della innovazione nel tempo
title('innovazione y_d');                                                  %Inserimento del titolo
subplot(2,1,2)                                                             %inserisci grafico nella seconda sezione (di sotto)
autocorr(squeeze(log_EKF_s.e(:,1,:)));                                     %per verificare la bianchezza di e se ne grafica l'autocorrelazione (impulso in 0 per un esatto wn)

figure;                                                                    %Creo una nuova figura vuota
subplot(2,1,1)                                                             %inserisci grafico su prima delle due sezioni
plot(time1,log_EKF_s.e(:,2));xlabel('tempo [s]');ylabel('e(2) [rad]')      %grafica andamento di seconda componente della innovazione nel tempo
title('innovazione y_ alfa');                                              %Inserimento del titolo
subplot(2,1,2)                                                             %inserisci grafico nella seconda sezione (di sotto)
autocorr(squeeze(log_EKF_s.e(:,2,:)));                                     %per verificare la bianchezza di e se ne grafica l'autocorrelazione (impulso in 0 per un esatto wn)

%===== Richiesta per avviare la simulazione di confronto tra i 3 casi =====%
prompt="Avviare simulazione di confronto Reale vs Stima vs Ricostruito? Y/N [Y]:";
txt=input(prompt,"s");
if isempty(txt)
    txt='Y';                       %se premi direttamente invio è visto come un consenso
end

if txt=='Y'
    disp('Ok simulo')
    run("Simulazione_Smoother.m")  %se desiderato esegue lo script per graficare sovrapposti il sistema reale, stimato e regolarizzato
else
    disp('Ok niente simulazione')
end

%=== Altro ===%
 %plot(time1,squeeze(log_results.tr(1,:,:)  %valore singolare massimo della
 %matrice di guadagno di correzione