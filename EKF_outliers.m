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
%    x_hat=[-pi/2;-pi/2;0;0];                              %prova per filtro che parte da condizione diversa rispetto a sistema
dev_std_0=0.05;                                 %deviazione standard iniziale, legata ad affidabilità della info a priori
P=dev_std_0^2*eye(4);                           %matrice di MSE iniziale per avviare l'algoritmo 

dt = 0.001;                                     %[s] tempo campionamento del filtro
dt_sens_y_d1 = 0.001;                           %[s] tempo campionamento del sensore 2
dt_sens_y_alfa = 0.001;                         %[s] tempo campionamento del sensore 1
StopTime=10;                                    %[s] tempo di simulazione
t_max = StopTime;                               %leggero cambio di notazione 
%=== Parametri per caso con outliers ===%
amp_out_yd=2;                                   %[m]
amp_out_yalfa= pi/2;                            %[rad]
period_out_yd= 4;                               %[s] ogni quanto si manifesta outliers su sensore 1
period_out_yalfa= 3;                            %[s] ogni quanto si manifesta outliers su sensore 2
pulse_width_out_yd= 5;                          %[%period] per quanto dura l'impulso dato come outliers
pulse_width_out_yalfa= 5;                       %[%period] per quanto dura l'impulso dato come outliers
num_out=0;                                      %Inizializza il conteggio degli outliers
soglia_Maha=5;                                  %soglia per identificare un outliers

open_system('ISI_Simulink'); %apre il file dello schema Simulink così lo si può eseguire

%% 2.Estrazione variabili da Simulink
teta1 = out.teta1.data;   %estraggo i valori assunti da teta1 durante la simulazione (tutte le variabili salvate come TimeSeries)
teta2 = out.teta2.data;   %estraggo i valori assunti da teta2

y_d1 = out.y_d1.data;     %estraggo i valori letti dal primo sensore
y_alfa = out.y_alfa.data; %estraggo i valori letti dal secondo sensore

dteta1 = out.dteta1.data; %estraggo i valori assunti da dteta1
dteta2 = out.dteta2.data; %estraggo i valori assunti da dteta2

time=out.time;            %estraggo anche il vettore dei tempi dal clock della simulazione campionato ogni dt
time1=time(2:end);        %stesso vettore dei tempi dove escludo il tempo 0 iniziale (serve per graficare  le stime ed errori che non sono definiti per t=0 non andando a salvare le condizioni iniziali dello stimatore)

%% 4.EKF

k=0; %Inizializzo la variabile k che conta le iterazioni del ciclo for (ad ogni ciclo predice e corregge con la misura attuale)
tic
  num_out_s1=0;
  num_out_s2=0;
  soglia_Maha=5;
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
    D_Maha_1=sqrt(e(1)'*(S(1,1))^(-1)*e(1));
    D_Maha_2=sqrt(e(2)'*(S(2,2))^(-1)*e(2));
%     if D_Maha_1<= soglia_Maha %no outliers in sensore 1
%         if D_Maha_2<=soglia_Maha %no outliers in sensore 2
%              x_hat = x_hat + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
%              P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
%         else %outliers in sensore 2
%             e_1=[1 0;0 0]*e; e=e_1;
%             x_hat = x_hat + L*(e_1);                      %calcolo della stima corretta k|k con stima BLUE ma scarto la innovazione del secondo sensore
%             P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
%        num_out_s2=num_out_s2+1;
%         end
%     else %outliers in sensore 1
%         if D_Maha_2<=soglia_Maha %no outliers in sensore 2
%             e_2=[0 0;0 1]*e; e=e_2;
%             x_hat = x_hat + L*(e_2);                      %calcolo della stima corretta k|k con stima BLUE ma scarto la innovazione del secondo sensore
%             P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
%         num_out_s1=num_out_s1+1;
%         else %outliers per entrambi i sensori
%             x_hat=x_hat;
%             P=P;
%             e=[0;0];
%             L=zeros(4,2);
%             num_out_s1=num_out_s1+1;
%             num_out_s2=num_out_s2+1;
%         end
%     end
    
    %Come in teoria:
    D_Maha=sqrt(e'*(S)^(-1)*e);
    if D_Maha<=soglia_Maha %no outliers nei sensori
             x_hat = x_hat + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
             P = (eye(4) - L*H)*P*(eye(4) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
    else %misura vista come outliers
         x_hat=x_hat;
            P=P;
            e=[0;0];
            L=zeros(4,2);
            num_out=num_out+1;
        end

    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    log_results.x_hat(k,1:4) = x_hat;
    log_results.P(:,:,k) = P;
    log_results.e(k,1:2) = e;
    log_results.tr(:,:,k)=svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end
num_out_s2
num_out_s1
num_out
toc
log_EKF = log_results;                                     %rinomino il luogo di salvataggio per evitare fraintendimenti futuri se dovessi salvare altro

%=== Faccio partire uno script apposito per graficare i risultati ottenuti ed
%eventualmente avviare una simulazione di confronto ===%
run("Subplot_confronto.m");
%% 8.Index Funzioni
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