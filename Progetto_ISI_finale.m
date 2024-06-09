%% 1.Definizione dei parametri fisici

clear all
close all
clc

%===Parametri usati sono gli stessi implementati per il proetto di CSI====%

tau_m1=0.3;                %[Nm] coppia di attuazione sul primo giunto
D=2;                       %[m] distanza orizzontale dal muro
l1=0.5;                    %[m] lunghezza prima asta
l2=0.4;                    %[m] lunghezza seconda asta
h1=0.35;                   %[m] posizione del centro di massa del primo link
h2=0.181;                  %[m] posizione del centro di massa del secondo link
m1=3.25;                   %[kg] massa prima asta
m2=1.90;                   %[kg] massa seconda asta
I1=0.654;                  %[kg*m^2] momento inerzia rispetto al centro di massa del primo link
I2=0.117;                  %[kg*m^2] momento inerzia rispetto al centro di massa del primo link
c1=8*10^(-1);              %[N*m*s] coeff. attrito viscoso sul primo giunto
c2=9*10^(-1);              %[N*m*s] coeff. attrito viscoso sul secondo giunto
g=9.81;                    %[m/s^2] accelerazione di gravità
%Salvo i parametri dinamici da riusare in seguito nelle funzioni
param.M = [m1, m2];
param.l = [l1, l2];
param.h = [h1, h2];
param.I = [I1, I2];
param.C = [c1, c2];

%=== Matrici della dinamica ===%

G=[1;0];
T=[1 -1;0 1];
N = [c1+c2, -c2;
    -c2,  c2];

%=== Rumore su sensori ===%

mean_y_d1 = 0 ;          %[m] media rumore su sensore 2
mean_y_alfa = 0 ;        %[rad] media rumore su sensore 1
sigma_y_d1 = 0.05 ;      %[m] deviazione standard rumore su sensore 2, 5 cm
sigma_y_alfa = 0.05 ;    %[rad] deviazione standard rumore su sensore 1,  gradi circa

%Seguendo le ipotesi di Kalman per sviluppare il filtro modello tutte le non
%idealità nella dinamica e nel modello di osservazione come rumori bianchi
%indipendenti a media nulla e varianza Q e R (rispettivamente):
dev_std_dtheta = 1;                                           %[rad/s] dev. stnd. di w1 e w2 errore causato dall'integrazione è pari a dt*w1 e dt*w2
dev_std_tau_d=0.3;                                            %[N*m] dev.stnd. di w3 e w4 disturbi di processo
Q = blkdiag(dev_std_dtheta^2*eye(2),dev_std_tau_d^2*eye(2));  % matrice di covarianza vettore di disturbo w (assunto wn)
dev_std_v_d=sigma_y_d1;                                       %richiamo la dv.stnd. rumore di misura per i due sensori
dev_std_v_alfa=sigma_y_alfa;
R = diag([dev_std_v_d^2, dev_std_v_alfa^2]);                  % matrice di covarianza vettore dei rumori v=[v_d;v_alfa] (assunto wn)

%=== Inizializzazione della simulazione ===%

%Definizione delle condizioni iniziali del sistema reale
teta1_0=0.5;  %[rad]
teta2_0=0.5;  %[rad]
dteta1_0=0;   %[rad/s]
dteta2_0=0;   %[rad/s]
%Definizione delle condizioni iniziali per il filtro NL
x_hat=[teta1_0; teta2_0; dteta1_0; dteta2_0];   %informazione a priori [rad;rad;rad/s;rad/s]
%      x_hat=[0;0;0;0];                              %prova per filtro che parte da condizione diversa rispetto a sistema
dev_std_0=0.05;                                 %deviazione standard iniziale, legata ad affidabilità della info a priori
P=dev_std_0^2*eye(4);                           %matrice di MSE iniziale per avviare l'algoritmo

dt = 0.001;                                     %[s] tempo campionamento del filtro
dt_sens_y_d1 = 0.001;                           %[s] tempo campionamento del sensore 2
dt_sens_y_alfa = 0.001;                         %[s] tempo campionamento del sensore 1
StopTime=10;                                    %[s] tempo di simulazione
t_max = StopTime;                               %leggero cambio di notazione

%=== Parametri per caso con outliers ===%
amp_out_yd_min=2;         amp_out_yd_max=5;           %[m] ampiezza max e min outlier sensore 1
amp_out_yalfa_min= pi/2;  amp_out_yalfa_max=pi;       %[rad] ampiezza max e min outliers sensore 2
period_out_yd= 0.5+(10-0.5)*rand(1,1);                %[s] ogni quanto si manifesta outliers su sensore 1, come numero random tra 0.5 e 10 secondi
period_out_yalfa= 0.5+(10-0.5)*rand(1,1);             %[s] ogni quanto si manifesta outliers su sensore 2, come numero random tra 0.5 e 10 secondi
pulse_width_out_yd= 3+(20-3)*rand(1,1);               %[%period] per quanto dura l'impulso dato come outliers, come numero random tra 3 e 20
pulse_width_out_yalfa= 3+(20-3)*rand(1,1);            %[%period] per quanto dura l'impulso dato come outliers, come numero random tra 3 e 20
num_out_s1=0;
num_out_s2=0;                                         %Inizializza il conteggio degli outliers
num_out=0;                                            %Inizializza il conteggio degli outliers per il caso "Come in teoria"
soglia_Maha=5;                                        %soglia per identificare un outliers, se si parte da c.i. distanti col filtro rispetto al reale, deve esere alzata sennò taglia fuori misure che non sono outliers ma non troppo sennò non riconosce i veri outliers

open_system('ISI_Simulink');                    %apre il file dello schema Simulink così lo si può eseguire

%% 2.Estrazione variabili da Simulink
teta1 = out.teta1.data;   %estraggo i valori assunti da teta1 durante la simulazione (tutte le variabili salvate come TimeSeries)
teta2 = out.teta2.data;   %estraggo i valori assunti da teta2

y_d1 = out.y_d1.data;     %estraggo i valori letti dal primo sensore
y_alfa = out.y_alfa.data; %estraggo i valori letti dal secondo sensore

dteta1 = out.dteta1.data; %estraggo i valori assunti da dteta1
dteta2 = out.dteta2.data; %estraggo i valori assunti da dteta2

time=out.time;            %estraggo anche il vettore dei tempi dal clock della simulazione campionato ogni dt
time1=time(2:end);        %stesso vettore dei tempi dove escludo il tempo 0 iniziale (serve per graficare  le stime ed errori che non sono definiti per t=0 non andando a salvare le condizioni iniziali dello stimatore)


%% 3.Simulazione doppio pendolo
%Compare una richiesta nella Comand Window per scegliere se avviare una
%simulazione rappresentativa della evoluzione del doppio pendolo
prompt="Avviare simulazione doppio pendolo? Y/N [Y]:";
txt=input(prompt,"s");
% Oltre a 'Y' anche premere solo invio è visto come un consenso
if isempty(txt)
    txt='Y';
end
if txt=='Y'
    disp('Ok simulo')

    %=== Simulazione doppio pendolo ===%
    figure(); hold on; xlim([-1,1]);xlabel('orizzontale [m]') ; ylim([-1,1]); ylabel('verticale [m]'); grid on; axis equal; %Creazione della figura vuota
    p0=[0;0];                                                           %Definizione dell'origine
    plot_origine = plot(p0(1),p0(2),'or');                              %Plot dell'origine
    N_=size(teta2,1);                                                   %Quanti campioni plotta
    a1=[-1;0];a2=[1;0];a3=[0;-1];a4=[0;1]; plot_lim=plot(a1(1),a1(2),'k',a2(1),a2(2),'k',a3(1),a3(2),'k',a4(1),a4(2),'k'); %punti a estremi del grafico, servono per tenere fissa la camera durante la simulazione
    for k=1:N_
        q=[teta1(k),teta2(k)];                                          %variabili di giunto aggiornate per ogni iterazione del ciclo
        [p1,p2]=cinematica(p0,q(1),q(2),l1,l2);                         %sfrutto la funzione "cinematica" per calcolare le posizioni dei due giunti
        %Vado a plottare i due link (dandone gli estremi) e i due giunti
        plot_link1=plot([p0(1),p1(1)],[p0(2),p1(2)],'b');
        plot_giunto2=plot(p1(1),p1(2),'or');
        plot_link2=plot([p1(1),p2(1)],[p1(2),p2(2)],'b');
        plot_ee=plot(p2(1),p2(2),'*k');
        timer=text(-0.9,0.9,"Timer: "+num2str(time(k),2));              %inserisco nel plot un piccolo timer per aver conto del tempo di simulazione
        pause(0.001);                                                   %passano 0.001 secondi (frame rate) prima di cancellare i grafici attuali e sostituirli con quelli della nuova iterazione del ciclo for
        %cancello tutti i grafici attualmente presenti a meno dell'origine
        delete(plot_link1);
        delete(plot_giunto2);
        delete(plot_link2);
        delete(plot_ee);
        delete(timer);
    end

else
    disp('Ok niente simulazione') %messaggio che compare se non simuli
end

%% 4.EKF

k=0; %Inizializzo la variabile k che conta le iterazioni del ciclo for (ad ogni ciclo predice e corregge con la misura attuale)
tic
for t = dt:dt:t_max
    k = k + 1; %Aggiorno k

    %============ Fase di Predizione =============%

    %Calcolo le matrici jacobiane della dinamica tc (f) (F=df/dx ;D_predict=df/dw)
    [F,D_predict] = jacobiano_f(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    %Aggiono F e D_predict per trovare le jacobiane della dinamica
    %discretizzata x_k+1= x_k+dt*f == F*x_k+ D_predict*w_k  (sono queste le
    %jacobiane che voglio in realtà)
    F = eye(4)+dt*F;            %jacobiana di x_k+1 (din.td) ripetto a x
    D_predict = dt*D_predict;   %jacobiana di x_k+1 (din.td) ripetto a w
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td ('eq_dinamica' è la dinamica tc) (è x_k|k-1)
    x_hat = x_hat + dt * eq_dinamica(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    P = F*P*F'+D_predict*Q*D_predict'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%

    %Calcolo le matrici jacobiane del modello di osservazione H=dh/dx ;M=dh/dv
    [H,M] = jacobiano_h(x_hat,param);
    %Calcolo la innovazione e riporto la seconda componente (legata a misura angolare) tra -pi e pi
    e = [y_d1(k); y_alfa(k)] - [D-l1*sin(x_hat(1))-l2*sin(x_hat(2)); pi+x_hat(1)-x_hat(2)];
    e = [e(1); atan2(sin(e(2)),cos(e(2)))];               %innovazione [m; rad]
    S = H*P*H' + M*R*M';                                  %matrice di covarianza della innovazione
    L = P*H'*S^(-1);                                      %guadagno di correzione
    x_hat = x_hat + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
    P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    log_results.x_hat(k,1:4) = x_hat;
    log_results.P(:,:,k) = P;
    log_results.e(k,1:2) = e;
    log_results.tr(:,:,k)=svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end

toc
log_EKF = log_results;                                     %rinomino il luogo di salvataggio per evitare fraintendimenti futuri se dovessi salvare altro

%=== Faccio partire uno script apposito per graficare i risultati ottenuti ed
%eventualmente avviare una simulazione di confronto ===%
run("Subplot_confronto.m");


%% 5.UKF

% !!! Rieseguire prima sezione 1., Simulink e 2. prima di avviarla (così
% sei sicuro di non avere variabili doppie) !!! %
tic

k = 0;                %inizializzo variabile per contare cicli del filtro

for t = dt:dt:t_max   %eseguo il filtro dopo ogni dt fino al tempo t_max

    k = k + 1;        %aggiorno variabile per conteggio cicli


    %============== Fase di Predizione =========================%
    %applico UT(x'_k|k)--> momenti di x_k+1|k con x'=[x;w] avendo disturbi
    %(4) non additivi (x' di dimensione 4+4)

    %Inizializzo algoritmo con media e covarianza a priori di x'
    prior_mean=[x_hat;0;0;0;0];              %media a priori di x'
    prior_cov=[P,zeros(4,4);zeros(4,4),Q];   %matrice covarianza a priori per x'

    %Definisco i parametri per la UT
    alpha = 1;                               % trasformata UT non scalata
    beta = 2;                                % valore ottimo per il caso gaussiano (min. varianza errore stima)
    kappa = 0;                               % lambda = 0 [peso del sigma pint centrale]
    enne = size(prior_mean,1);               %dimensione di x'
    lambda = alpha^2*(enne + kappa) - enne;  %parametro in funzione di cui calcolo i pesi dei sigma point

    %Calcolo dei pesi
    w0_pred = lambda/(lambda + enne);             %peso centrale su mendia
    wc_pred(1) = w0_pred + 1 - alpha^2 + beta;    %peso centrale da usare per stima momenti secondo ordine
    wc_pred(2:2*enne+1) = 1/2/(enne+lambda);      %tutti i pesi su altri sigma point sono tutti uguali e tali da dare somma 1

    wm_pred = wc_pred;                                 %pesi su sigma-point per i=\ 0 sono gli stessi per le medie campionarie definenti i momenti del primo e secondo ordine dei sigma-point propagati
    wm_pred(1) = w0_pred;                              %sovrascrivo il peso sul primo sigma-point per la media rimettendolo pari a w0

    %Fattorizzo la matrice di covarianza tramite la SVD
    [U,Sigma,~] = svd(prior_cov);
    GAMMA = U*Sigma^(1/2);

    % genero i 2n+1 sigma points
    sigma_points = prior_mean;              %primo s.p. concide con media a priori
    for i = 1:size(GAMMA,2)                 %tutti gli altri s.p. sono simmetrici rispetto alla media spostandomi lungo le direzioni date da fattorizzazione di covarianza a priori
        sigma_points(:,i+1)      = prior_mean + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points(:,i+1+enne) = prior_mean - sqrt(enne+lambda)*GAMMA(:,i);
    end

    %propagazione dei sigma points generati tramite eq_td
    for i = 1:size(sigma_points,2)
        %uso variabili locali per brevità
        x1 = sigma_points(1,i);
        x2 = sigma_points(2,i);
        x3 = sigma_points(3,i);
        x4 = sigma_points(4,i);

        w1 = sigma_points(5,i);
        w2 = sigma_points(6,i);
        w3 = sigma_points(7,i);
        w4 = sigma_points(8,i);
        %applico la dinamica td per propagare i sigma points
        [propagated_sigma_points] = eq_td (x1, x2, x3, x4, tau_m1, w1, w2, w3, w4, param, g, N, T, G, dt);
        %salvo assieme tutti i loro valori
        save.propagated_sigma_points(:,i) = propagated_sigma_points;
        propagated_sigma_points = save.propagated_sigma_points;

    end

    % Calcolo la media pesata dei s.p. propagati per ottenere approssimazione del momento del primo ordine di x_k+1|k (NOTE: mean of angles)
    measurements_mean(1,1) = atan2( sin(propagated_sigma_points(1,:))*wm_pred' , cos(propagated_sigma_points(1,:))*wm_pred');
    measurements_mean(2,1) = atan2( sin(propagated_sigma_points(2,:))*wm_pred' , cos(propagated_sigma_points(2,:))*wm_pred');
    measurements_mean(3,1) = propagated_sigma_points(3,:)*wm_pred';
    measurements_mean(4,1) = propagated_sigma_points(4,:)*wm_pred';

    %Trovo lo scarto rispetto alla media campionaria dei s.p.p. (NOTE: difference of angles)
    measurements_tilde=[];       %inizializzo la variabile
    for i=1:size(sigma_points,2)
        measurements_tilde(1,i) = atan2(sin(propagated_sigma_points(1,i) -measurements_mean(1,1)),cos(propagated_sigma_points(1,i) -measurements_mean(1,1)));
        measurements_tilde(2,i) = atan2(sin(propagated_sigma_points(2,i) -measurements_mean(2,1)),cos(propagated_sigma_points(2,i) -measurements_mean(2,1)));
        measurements_tilde(3,i) = propagated_sigma_points(3,i)-measurements_mean(3,1);
        measurements_tilde(4,i) = propagated_sigma_points(4,i)-measurements_mean(4,1);
    end

    %Calcolo il momento del secondo ordine di x_k+1|k come covarianza
    %campionaria pesata
    measurements_cov = zeros(4,4);   %inizializzo la variabile
    for i = 1:size(sigma_points,2)
        measurements_cov = measurements_cov+ wc_pred(i)*measurements_tilde(:,i)*measurements_tilde(:,i)';
    end

    %rinomino i momenti stimati così da usarli per la fase successiva di
    %correzione
    x_hat=measurements_mean;
    P=measurements_cov;


    %=============== Fase di Correzione ==============%
    %applico UT(x_k+1|k)--> momenti di y_k+1|k+1 (a meno del contributo del rumore di misura additivo)

    %Inizializzo algoritmo con media e covarianza a priori di x trovate dal
    %passo di predizione precedente
    prior_mean =x_hat;
    prior_cov = P;

    %Definisco i parametri per la UT
    alpha = 1;                               % trasformata UT non scalata
    beta = 2;                                % valore ottimo per il caso gaussiano (min. varianza errore stima)
    kappa = 0;                               % lambda = 0 [peso del sigma pint centrale]
    enne = size(prior_mean,1);               %dimensione di x'
    lambda = alpha^2*(enne + kappa) - enne;  %parametro in funzione di cui calcolo i pesi dei sigma point

    %Calcolo dei pesi
    w0_corr = lambda/(lambda + enne);                  %peso centrale su mendia
    wc_corr(1) = w0_corr + 1 - alpha^2 + beta;         %peso centrale da usare per stima momenti secondo ordine
    wc_corr(2:2*enne+1) = 1/2/(enne+lambda);           %parametro in funzione di cui calcolo i pesi dei sigma point
    wm_corr = wc_corr;                                 %pesi su sigma-point per i=\ 0 sono gli stessi per le medie campionarie definenti i momenti del primo e secondo ordine dei sigma-point propagati
    wm_corr(1) = w0_corr;                              %sovrascrivo il peso sul primo sigma-point per la media rimettendolo pari a w0

    %Fattorizzo la matrice di covarianza tramite la SVD
    [U,Sigma,~] = svd(prior_cov);
    GAMMA = U*Sigma^(1/2);

    % genero i 2n+1 sigma points
    sigma_points_corr = prior_mean;   %primo s.p. concide con media a priori
    for i = 1:size(GAMMA,2)           %tutti gli altri s.p. sono simmetrici rispetto alla media spostandomi lungo le direzioni date da fattorizzazione di covarianza a priori
        sigma_points_corr(:,i+1)      = prior_mean + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points_corr(:,i+1+enne) = prior_mean - sqrt(enne+lambda)*GAMMA(:,i);
    end

    %propagazione dei sigma points generati tramite modello di osservazione
    for i = 1:size(sigma_points_corr,2)
        propagated_sigma_points_corr(1,i) = D - l1*sin(sigma_points_corr(1,i)) - l2*sin(sigma_points_corr(2,i));
        propagated_sigma_points_corr(2,i) = pi + sigma_points_corr(1,i) - sigma_points_corr(2,i);
        propagated_sigma_points_corr(2,i)=atan2(sin(propagated_sigma_points_corr(2,i)),cos(propagated_sigma_points_corr(2,i))); %riporto tra -pi e pi la seconda componente che è angolare
    end

    % Calcolo la media pesata dei s.p. propagati per ottenere approssimazione del momento del primo ordine di y_k+1|k+1 (NOTE: mean of angles)
    virtual_measurements_mean(1,1) = propagated_sigma_points_corr(1,:)*wm_corr';
    virtual_measurements_mean(2,1) = atan2( sin(propagated_sigma_points_corr(2,:))*wm_corr' , cos(propagated_sigma_points_corr(2,:))*wm_corr');

    %Trovo lo scarto rispetto alla media campionaria dei s.p.p. (NOTE: difference of angles)
    for i=1:size(sigma_points_corr,2)
        virtual_measurements_tilde(1,i) = propagated_sigma_points_corr(1,i)-virtual_measurements_mean(1,1);
        virtual_measurements_tilde(2,i) = atan2(sin(propagated_sigma_points_corr(2,i) -virtual_measurements_mean(2,1)),cos(propagated_sigma_points_corr(2,i) -virtual_measurements_mean(2,1)));
    end

    %Calcolo il momento del secondo ordine di y_k+1|k+1 (a meno di efetti del rumore) come covarianza
    % e cross-covarianza campionaria pesata
    virtual_measurements_cov = zeros(2,2);  %inizializzo la variabile
    cross_cov = zeros(4,2);                 %inizializzo la variabile
    for i = 1:size(sigma_points_corr,2)
        virtual_measurements_cov = virtual_measurements_cov + wc_corr(i)*virtual_measurements_tilde(:,i)*virtual_measurements_tilde(:,i)';
        cross_cov = cross_cov + wc_corr(i)*(sigma_points_corr(:,i) - prior_mean)*virtual_measurements_tilde(:,i)';
    end

    %Salvo e correggo i momenti così da inglobare il rumore di misura
    %additivo
    S=virtual_measurements_cov+R;     %correzioni per tenere conto del rumore v additivo
    y_hat=virtual_measurements_mean;  %rumore a media nulla non agisce su momento primo ordine
    Pxy=cross_cov;                    %rumore non agisce neanche sul momento incrociato del secondo ordine

    %trovati i momenti di y_k+1|k+1 applico la stima BLUE per aggiornare la
    %stima dello stato con la misura più recente
    e = [y_d1(k); y_alfa(k)] - y_hat;         % innovazione
    e = [e(1); atan2(sin(e(2)),cos(e(2)))];   %seconda componente è angolare quindi la riporto tra -pi e pi
    L = Pxy*S^(-1);                           %guadagno di correzione
    %stima BLUE
    x_hat = x_hat + L*e;                      %stato corretto in funzione del predetto e innovazione
    P = P - L*S*L';                           %varianza sulla stima corretta diminuisce avendo aggiunto nuove info

    %Salvo i risultati per avere i dati pronti da usare per il passo di predizione del
    %nuovo ciclo di predizione
    log_results.x_hat(k,1:4) = x_hat;
    log_results.P(:,:,k) = P;
    log_results.e(k,1:2) = e;
    log_results.tr(:,:,k)=svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end
toc

%=== Faccio partire uno script apposito per graficare i risultati ottenuti ed
%eventualmente avviare una simulazione di confronto ===%
run("Subplot_confrontoUKF.m");

%% 6.SMOOTHER

% !!! Rieseguire prima sezione 1., Simulink e 2. prima di avviarla (così
% sei sicuro di non avere variabili doppie) !!! %

%Prima fase da eseguire è una stima in avanti usando un EKF, segue la fase
%di regolarizzazione con la ricorsione in indietro per ricostruire lo stato
%al passo k usando tutte le osservazioni N già disponibili.


%================== Fase in avanti con EKF ==================%

%Il seguente script è analogo a quello trovato nella sezione 4.EKF, viene
%riportato in modo tale che questa sezione 6.Smoother sia il più indipendente possibile
k=0;

for t = dt:dt:t_max

    k = k + 1;

    %========== Fase di Predizione ==========%

    [F,D_predict] = jacobiano_f(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    F = eye(4)+dt*F;
    D_predict = dt*D_predict;
    x_hat = x_hat + dt * eq_dinamica(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    P = F*P*F'+D_predict*Q*D_predict';
    %Per corretto funzionamento dello smoother devo salvare anche le variabili
    %nella fase di predizione (Aggiunta rispetto sezione 4.)%
    log_EKF_s.F(:,:,k) = F;
    log_EKF_s.x_hat_prediction(:,1,k) = x_hat;
    log_EKF_s.P_prediction(:,:,k) = P;

    %========== Fase di Correzione ==========%

    [H,M] = jacobiano_h(x_hat,param);
    e = [y_d1(k); y_alfa(k)] - [D-l1*sin(x_hat(1))-l2*sin(x_hat(2)); pi+x_hat(1)-x_hat(2)];
    e = [e(1); atan2(sin(e(2)),cos(e(2)))];
    S = H*P*H' + M*R*M';
    L = P*H'*S^(-1);
    x_hat = x_hat + L*e;
    P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';

    log_EKF_s.e(k,1:2) = e;
    log_EKF_s.x_hat_correction(:,1,k) = x_hat;
    log_EKF_s.P_correction(:,:,k) = P;
    log_results.tr(:,:,k)=svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end


% ============Inizializzazione dello Smoother ============%

%Devo imporre che  x_n|n de P_n|n siano uguali agli ultimi valori di x e P corretti
n=size(log_EKF_s.x_hat_correction,3);
log_EKF_s.x_hat_smoothed(:,1,n) =log_EKF_s.x_hat_correction(:,1,n) ; %è più corretto inizializzare con n e non k
log_EKF_s.P_smoothed(:,:,n) = log_EKF_s.P_correction(:,:,n);


%============= Fase in indietro con Smoother =============%

%Devo fare una ricorsione da enne-1 ; ogni -1 ; fino a 1

for k = n-1:-1:1

    % Ck = P_k|k * F_k+1' * P_k+1|k^(-1) % ATTENZIONE non è Fk+1 ma Fk !!!
    Ck = log_EKF_s.P_correction(:,:,k) ...
        * log_EKF_s.F(:,:,k)'/ ...
        log_EKF_s.P_prediction(:,:,k+1);
    log_EKF.Ck(:,:,k)=Ck;  %Calcolo e salvo la matrice C_k necessaria per la ricorsione all'indietro

    % x_k|n = x_k|k + Ck*(x_k+1|n - x_k+1|k)  Calcolo e salvo la stima dello stato al passo k
    log_EKF_s.x_hat_smoothed(:,1,k) = log_EKF_s.x_hat_correction(:,1,k) +Ck*( log_EKF_s.x_hat_smoothed(:,1,k+1) - log_EKF_s.x_hat_prediction(:,1,k+1) );
    log_EKF.diff(:,:,k)=log_EKF_s.x_hat_smoothed(:,1,k+1)-log_EKF_s.x_hat_prediction(:,1,k+1);

    % P_k|n = P_k|k + Ck*(P_k+1|n - P_k+1|k)*Ck'  Calcolo e salvo la
    % matrice MSE trovata al passo k
    log_EKF_s.P_smoothed(:,:,k) = log_EKF_s.P_correction(:,:,k)...
        + Ck*log_EKF_s.P_smoothed(:,:,k+1)...
        - log_EKF_s.P_prediction(:,:,k+1)*Ck';
end

%=== Eseguo un apposito script per l'analisi dei risultati ottenuti e
%simulazione per confrontare le evoluzioni reale, stimata e ricostruita del
%sistema ===%
run("Plot_comparison_smoother")

%% 7. EKF con diversi tempi di campionamento (modifica un dt_sens_ e rimanda la simulazione)

k=0;                  %conta iterazioni filtro
s1=0;                 %conta misure sensore 1
s2=0;                 %conta misure sensore 2
ts1=out.y_d1.Time;    %estraggo i tempi a cui il sensore 1 ha dato delle misure
ts2=out.y_alfa.Time;  %estraggo i tempi a cui il sensore 2 ha dato delle misure

for t = dt:dt:t_max
    k = k + 1; %Aggiorno k

    %============ Fase di Predizione =============%

    %Calcolo le matrici jacobiane della dinamica tc (f) (F=df/dx ;D_predict=df/dw)
    [F,D_predict] = jacobiano_f(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    %Aggiono F e D_predict per trovare le jacobiane della dinamica
    %discretizzata x_k+1= x_k+dt*f == F*x_k+ D_predict*w_k  (sono queste le
    %jacobiane che voglio in realtà)
    F = eye(4)+dt*F;            %jacobiana di x_k+1 (din.td) ripetto a x
    D_predict = dt*D_predict;   %jacobiana di x_k+1 (din.td) ripetto a w
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td ('eq_dinamica' è la dinamica tc) (è x_k|k-1)
    x_hat = x_hat + dt * eq_dinamica(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    P = F*P*F'+D_predict*Q*D_predict'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%

    %Calcolo le matrici jacobiane del modello di osservazione H=dh/dx ;M=dh/dv
    [H,M] = jacobiano_h(x_hat,param);
    [i,~,val1]=find((k-1)*dt<ts1 & ts1<=k*dt);  %i=posizione elemento di ts1 che rispetta la specifica (c'è una misura da sensore 1 arrivata durante il dt damiterazione precedente ad ora), val=0 se non trova tale indice
    [j,~,val2]=find((k-1)*dt<ts2 & ts2<=k*dt);  %j=posizione elemento di ts2 che rispetta la specifica (c'è una misura da sensore 2 arrivata durante il dt damiterazione precedente ad ora), val=0 se non trova tale indice
    if isempty(val1)                            %nel caso NON ci sia misura 1 disponibile può usare solo la riga di correzione relativa al sensore 2
        if isempty(val2)                        %se non ci fosse nenache la misura 2 non corregge proprio in questa iterazione
            x_hat = x_hat;
            P = P;
            e=[0;0]; L=zeros(4,2);              %non facendo correzione e ed L non ci sono, le devo porre a 0 per evitare errori nel salvataggio (quando esegue questa sezione non troverebbe e ed L da salvare altrimenti)
        else                                    % ha solo misura 2
            H = H(2,:);
            M = M(2,:);
            e_2 = (y_alfa(j)) - (pi+x_hat(1)-x_hat(2)); e_2 = atan2(sin(e_2),cos(e_2)); e=[0;e_2]; %devo mettere la prima componente a 0 sennò si salvano assieme scalari e vettori (quando usa entrambe le misure)
            S = H*P*H' + M*R*M';
            L = P*H'*S^(-1);
            x_hat = x_hat + L*e_2;
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';

        end
    else                                       %la misura 1 è disponibile, sicuramente uso almeno la riga 1 per correggere
        if isempty(val2)                       %non trova misura 2 nell'intervallo, uso SOLO la info della misura 1 per la correzione
            H = H(1,:);
            M = M(1,:);
            e_1 =( y_d1(i)) - (D-l1*sin(x_hat(1))-l2*sin(x_hat(2))); e=[e_1;0]; %devo mettere la seconda componente a 0 sennò si salvano assieme scalari e vettori (quando usa entrambe le misure)
            S = H*P*H' + M*R*M';
            L = P*H'*S^(-1);
            x_hat = x_hat + L*e_1;
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';

        else                                  %ha entrambe le misure, faccio una correzione completa
            e_1 = (y_d1(i)) - (D-l1*sin(x_hat(1))-l2*sin(x_hat(2)));
            e_2 = (y_alfa(j)) - (pi+x_hat(1)-x_hat(2));
            e = [e_1; atan2(sin(e_2),cos(e_2))];
            S = H*P*H' + M*R*M';
            L = P*H'*S^(-1);
            x_hat = x_hat + L*e;
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';

        end
    end
    %salvo tutti i risultati ad ogni iterazoine così dapoterloi graficare
    %successivamente
    log_results.x_hat(k,1:4) = x_hat;
    log_results.P(:,:,k) = P;
    log_results.e(k,1:2) = e;
    log_results.tr(:,:,k)=svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end

log_EKF = log_results;             %rinomino il luogo di salvataggio

%=== Faccio partire uno script apposito per graficare i risultati ottenuti ed
%eventualmente avviare una simulazione di confronto ===%
run("Subplot_confronto.m");

%% 8.EKF con outliers
%Vedi nota riga 568
k=0;            %Inizializzo la variabile k che conta le iterazioni del ciclo for (ad ogni ciclo predice e corregge con la misura attuale)
tic

for t = dt:dt:t_max
    k = k + 1; %Aggiorno k

    %============ Fase di Predizione =============%

    %Calcolo le matrici jacobiane della dinamica tc (f) (F=df/dx ;D_predict=df/dw)
    [F,D_predict] = jacobiano_f(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    %Aggiono F e D_predict per trovare le jacobiane della dinamica
    %discretizzata x_k+1= x_k+dt*f == F*x_k+ D_predict*w_k  (sono queste le
    %jacobiane che voglio in realtà)
    F = eye(4)+dt*F;            %jacobiana di x_k+1 (din.td) ripetto a x
    D_predict = dt*D_predict;   %jacobiana di x_k+1 (din.td) ripetto a w
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td ('eq_dinamica' è la dinamica tc) (è x_k|k-1)
    x_hat = x_hat + dt * eq_dinamica(x_hat,tau_m1,zeros(4,1),param,g,N,T,G);
    P = F*P*F'+D_predict*Q*D_predict'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%

    %Calcolo le matrici jacobiane del modello di osservazione H=dh/dx ;M=dh/dv
    [H,M] = jacobiano_h(x_hat,param);
    %Calcolo la innovazione e riporto la seconda componente (legata a misura angolare) tra -pi e pi
    e = [y_d1(k); y_alfa(k)] - [D-l1*sin(x_hat(1))-l2*sin(x_hat(2)); pi+x_hat(1)-x_hat(2)];
    e = [e(1); atan2(sin(e(2)),cos(e(2)))];                                %innovazione [m; rad]
    S = H*P*H' + M*R*M';                                                   %matrice di covarianza della innovazione
    L = P*H'*S^(-1);                                                       %guadagno di correzione
    D_Maha_1=sqrt(e(1)'*(S(1,1))^(-1)*e(1));
    D_Maha_2=sqrt(e(2)'*(S(2,2))^(-1)*e(2));
% commentare fino a riga 596 e decommentare fino a 609 per eseguire invece
% la condizione vettoriale invece che la doppia scalare (Mahalanobis
% diventa una norma euclidea scalata)
    if D_Maha_1<= soglia_Maha       %no outliers in sensore 1
        if D_Maha_2<=soglia_Maha    %no outliers in sensore 2
            x_hat = x_hat + L*e;                                           %calcolo della stima corretta k|k con stima BLUE
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';            %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
        else                        %outliers in sensore 2
            e_1=[1 0;0 0]*e; e=e_1;                                        %rinominazione necessaria per mettere a 0 (escludo anche dai grafici) la componente della innovazione affetta da outliers
            x_hat = x_hat + L*(e_1);                                       %calcolo della stima corretta k|k con stima BLUE ma scarto la innovazione del secondo sensore
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';            %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
            num_out_s2=num_out_s2+1;                                       %aggiorno il numero di outliers
        end
    else                           %outliers in sensore 1
        if D_Maha_2<=soglia_Maha   %no outliers in sensore 2
            e_2=[0 0;0 1]*e; e=e_2;                                        %rinominazione necessaria per mettere a 0 (escludo anche dai grafici) la componente della innovazione affetta da outliers
            x_hat = x_hat + L*(e_2);                                       %calcolo della stima corretta k|k con stima BLUE ma scarto la innovazione del primo sensore
            P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';            %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
            num_out_s1=num_out_s1+1;                                       %aggiorno il numero di outliers
        else                       %outliers per entrambi i sensori, non correggo
            x_hat=x_hat;
            P=P;
            e=[0;0];                                                       %rinominazione necessaria per mettere a 0 (escludo anche dai grafici) la componente della innovazione affetta da outliers
            L=zeros(4,2);                                                  %stesso discorso della riga sopra
            num_out_s1=num_out_s1+1;                                       %aggiorno il numero di outliers
            num_out_s2=num_out_s2+1;                                       %aggiorno il numero di outliers
            num_out=num_out+1;                                             %conta quante volte non correggo proprio per niente
        end
    end

%     %Come in teoria, usando la condizione vettoriale:
%     D_Maha=sqrt(e'*(S)^(-1)*e);
%     if D_Maha<=soglia_Maha %no outliers nei sensori
%         x_hat = x_hat + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
%         P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
%     else %misura vista come outliers, non corregge
%         x_hat=x_hat;
%         P=P;
%         e=[0;0];
%         L=zeros(4,2);
%         num_out=num_out+1;
%     end

    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    log_results.x_hat(k,1:4) = x_hat;
    log_results.P(:,:,k) = P;
    log_results.e(k,1:2) = e;
    log_results.tr(:,:,k)=svd(L);                          %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
    log_results.S(:,:,k)=S;
end
num_out_s2                                                 %mostra quanti outliers del sensore 2 sono stati trovati
num_out_s1                                                 %mostra quanti outliers del sensore 1 sono stati trovati
num_out                                                    %mostra, per il caso "Come in teoria" quanti outliers sono stati trovati   
toc
log_EKF = log_results;                                     %rinomino il luogo di salvataggio per evitare fraintendimenti futuri se dovessi salvare altro

%=== Faccio partire uno script apposito per graficare i risultati ottenuti ed
%eventualmente avviare una simulazione di confronto ===%
run("Subplot_confronto.m");


%% 9.Index Funzioni
%Raccolta e descrizione delle funzioni usate nello script

%=== Implementazione della dinamica del sistema ===%

function [ddteta1,ddteta2]=dynamics(tau_m1,tau_d,dteta1,dteta2,teta1,teta2)
%raccolo le velocità angolari per fare i conti con le matrici
dteta=[dteta1;dteta2];
%inserisco le matrici della dinamica
M=[I1+m1*h1^2+m2*l1^2 , l1*m2*h2*cos(teta1-teta2); l1*m2*h2*cos(teta1-teta2) , I2+m2*h2^2];
N=[c1+c2 , -c2; -c2 , c2];
q=[l1*m2*h2*sin(teta1-teta2)*dteta2^2-(m1*h1+m2*l1)*g*sin(teta1) ; -l1*m2*h2*sin(teta1-teta2)*dteta1^2-m2*h2*g*sin(teta2)];
%inversione della dinamica
ddteta=pinv(M)*(G*tau_m1+T*tau_d-N*dteta-q);
%separazione della accelerazione angolare nelle sue due componenti
ddteta1=ddteta(1);
ddteta2=ddteta(2);
end

%=== Equazione della cinematica per calcolare posizione dei giunti ===%

function [p1,p2]=cinematica(p0, q1, q2,l1,l2)
% q1=teta1;
% q2=teta2;
p1 = [l1*sin(q1); l1*cos(q1)];            %posizione primo giunto
p2 = p0 + p1 + [l2*sin(q2); l2*cos(q2)];  %posizione secondo giunto
end

%=== Funzione per graficare posizione stimata dei giunti ===%

function [o1,o2]=grafico_ekf(or, x_results1, x_results2, l1, l2)
p_1 = x_results1;
p_2 = x_results2;
o1 = [l1*sin(p_1); l1*cos(p_1)];
o2 = or + o1 + [l2*sin(p_2); l2*cos(p_2)];
end

%=== Calcolo del jacobiano in x e w della dinamica tempo continua ===%

function [F,D] = jacobiano_f(x,u,d,param,g,N,T,G)
%estraggo i valori dei parametri che mi servono dalla struttura 'param'
M1 = param.M(1);
M2 = param.M(2);
l1 = param.l(1);
h1 = param.h(1);
h2 = param.h(2);
I1 = param.I(1);
I2 = param.I(2);
C1 = N(1,1)-N(2,2);
C2 = N(2,2);
G1 = G(1,1);
G2 = G(2,1);
T1_1 = T(1,1);
T1_2 = T(1,2);
T2_1 = T(2,1);
T2_2 = T(2,2);
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
d1 = d(1);
d2 = d(2);
d3 = d(3);
d4 = d(4);

%Calcolo le F=df/dx e D=df/dw con il simbolico a parte e riporto i
%risultati dato che l'uso del simbolico rallenterebbe l'algoritmo
F = ...
    [                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                                                                                        1,                                                                                                                                                                                                                                                                                                                                   0;
    0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                   1;
    ((M2*h2^2 + I2)*(- M2*h2*l1*cos(x1 - x2)*x4^2 + g*cos(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2^2*h2^2*l1^2*x3^2*cos(x1 - x2)^2)/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (M2*h2*l1*sin(x1 - x2)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (2*M2^3*h2^3*l1^3*cos(x1 - x2)^2*sin(x1 - x2)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2 - (2*M2^2*h2^2*l1^2*cos(x1 - x2)*sin(x1 - x2)*(M2*h2^2 + I2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2,                          (M2*h2*l1*x4^2*cos(x1 - x2)*(M2*h2^2 + I2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*h2*l1*cos(x1 - x2)*(- M2*h2*l1*cos(x1 - x2)*x3^2 + M2*g*h2*cos(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (2*M2^3*h2^3*l1^3*cos(x1 - x2)^2*sin(x1 - x2)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2 - (M2*h2*l1*sin(x1 - x2)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (2*M2^2*h2^2*l1^2*cos(x1 - x2)*sin(x1 - x2)*(M2*h2^2 + I2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2,         - ((M2*h2^2 + I2)*(C1 + C2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*h2*l1*cos(x1 - x2)*(C2 + 2*M2*h2*l1*x3*sin(x1 - x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2),             ((M2*h2^2 + I2)*(C2 - 2*M2*h2*l1*x4*sin(x1 - x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (C2*M2*h2*l1*cos(x1 - x2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2);
    (M2*h2*l1*sin(x1 - x2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*h2*l1*cos(x1 - x2)*(- M2*h2*l1*cos(x1 - x2)*x4^2 + g*cos(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (M2*h2*l1*x3^2*cos(x1 - x2)*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (2*M2^3*h2^3*l1^3*cos(x1 - x2)^2*sin(x1 - x2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2 - (2*M2^2*h2^2*l1^2*cos(x1 - x2)*sin(x1 - x2)*(M1*h1^2 + M2*l1^2 + I1)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2, ((- M2*h2*l1*cos(x1 - x2)*x3^2 + M2*g*h2*cos(x2))*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2^2*h2^2*l1^2*x4^2*cos(x1 - x2)^2)/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*h2*l1*sin(x1 - x2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (2*M2^3*h2^3*l1^3*cos(x1 - x2)^2*sin(x1 - x2)*(- M2*h2*l1*sin(x1 - x2)*x4^2 + C2*x4 + T1_1*d3 + T1_2*d4 + G1*u - x3*(C1 + C2) + g*sin(x1)*(M1*h1 + M2*l1)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2 + (2*M2^2*h2^2*l1^2*cos(x1 - x2)*sin(x1 - x2)*(M1*h1^2 + M2*l1^2 + I1)*(M2*h2*l1*sin(x1 - x2)*x3^2 + C2*x3 + T2_1*d3 + T2_2*d4 - C2*x4 + G2*u + M2*g*h2*sin(x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)^2, ((C2 + 2*M2*h2*l1*x3*sin(x1 - x2))*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) + (M2*h2*l1*cos(x1 - x2)*(C1 + C2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2), - (C2*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*h2*l1*cos(x1 - x2)*(C2 - 2*M2*h2*l1*x4*sin(x1 - x2)))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)];

D = ...
    [1, 0,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0;
    0, 1,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0;
    0, 0,           (T1_1*(M2*h2^2 + I2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*T2_1*h2*l1*cos(x1 - x2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2),           (T1_2*(M2*h2^2 + I2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*T2_2*h2*l1*cos(x1 - x2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2);
    0, 0, (T2_1*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*T1_1*h2*l1*cos(x1 - x2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2), (T2_2*(M1*h1^2 + M2*l1^2 + I1))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2) - (M2*T1_2*h2*l1*cos(x1 - x2))/(- M2^2*h2^2*l1^2*cos(x1 - x2)^2 + M2^2*h2^2*l1^2 + M1*M2*h1^2*h2^2 + I1*M2*h2^2 + I2*M2*l1^2 + I2*M1*h1^2 + I1*I2)];

end

%=== Implementa la equazione dinamica in forma di stato tc ===%

function [f] = eq_dinamica(x,u,w,param,g,N,T,G)
%estraggo i valori dei parametri che mi servono dalla struttura 'param'
m1 = param.M(1);
m2 = param.M(2);
l1 = param.l(1);
h1 = param.h(1);
h2 = param.h(2);
I1 = param.I(1);
I2 = param.I(2);
%definisco delle variabili locali per maggiore chiarezza
teta   = [x(1);x(2)];
dteta  = [x(3);x(4)];
tau_d   = [w(3);w(4)];
tau_m   = u;
%Riporto le matrici dinamiche
M = [I1 + m1*h1^2 + m2*l1^2,      l1*m2*h2*cos(teta(1)-teta(2));
    l1*m2*h2*cos(teta(1)-teta(2)), I2+m2*h2^2];
q = [l1*m2*h2*sin(teta(1)-teta(2))*dteta(2)^2 - (m1*h1+m2*l1)*g*sin(teta(1));
    -l1*m2*h2*sin(teta(1)-teta(2))*dteta(1)^2 - m2*h2*g*sin(teta(2))];
%Creo la rappresentazione in forma di stato della dinamica (f è un vettore
%colonna)
f(1,1) = dteta(1) + w(1);
f(2,1) = dteta(2) + w(2);
f(3:4,1) = M\(-N*dteta - q + G*tau_m + T*tau_d);

end

%=== Calcolo del jacobiano in x e v del modello di osservazione h=[y_d1; y_alfa]=[D-l1*sen(teta1)-l2*sen(teta2)+v1 ;pi+teta1-teta2+v2] ===%

function [H,M] = jacobiano_h(x,param)
%definisco delle variabili locali per semplicità
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
%estraggo i parametri fisici necessari
l1 = param.l(1);
l2 = param.l(2);
%Riporto i jacobiani calcolati manualmente H=dh/dx e M=dh/dv
H = ...
    [-l1*cos(x1),   -l2*cos(x2),   0,   0;
    1,             -1,          0,   0];
M = ...
    [1, 0;
    0, 1];
end

%=== Dinamica discretizzata per propagare i sigma point di UKf ===%

function [f_td] = eq_td (x1, x2, x3, x4, u, w1, w2, w3, w4, param, g, N, T, G, dt)
%estraggo i parametri fisici necessari
m1 = param.M(1);
m2 = param.M(2);
l1 = param.l(1);
h1 = param.h(1);
h2 = param.h(2);
I1 = param.I(1);
I2 = param.I(2);
%definisco variabili locali per maggior chiarezza
teta   = [x1; x2];
dteta  = [x3; x4];
tau_d   = [w3; w4];
tau_m   = u;
%riporto le matrici della dinamica
M = [I1 + m1*h1^2 + m2*l1^2,      l1*m2*h2*cos(teta(1)-teta(2));
    l1*m2*h2*cos(teta(1)-teta(2)), I2+m2*h2^2];
q = [l1*m2*h2*sin(teta(1)-teta(2))*dteta(2)^2 - (m1*h1+m2*l1)*g*sin(teta(1));
    -l1*m2*h2*sin(teta(1)-teta(2))*dteta(1)^2 - m2*h2*g*sin(teta(2))];
%scrivo la dinamica in forma di stato in tempo continuo
f(1,1) = dteta(1) + w1;
f(2,1) = dteta(2) + w2;
f(3:4,1) = M\(-N*dteta - q + G*tau_m + T*tau_d);
%La dinamica discretizzata in forma di stato sarà dunque:
f_td = [x1+dt*(f(1,1));
    x2+dt*(f(2,1));
    x3+dt*f(3,1);
    x4+dt*f(4,1)];

end