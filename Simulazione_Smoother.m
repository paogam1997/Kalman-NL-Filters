%% Simulazione Smoother

%Script va a sovrapporre le simulazioni del sistema vero, stimato e regolarizzato per
%renderci conto meglio della bontà della stima e dello smoother

%Creo figura vuota, inserisco la griglia, il titolo e pongo la dimensione
%degli assi
figure(); hold on; xlim([-1,1]);xlabel('orizzontale [m]'); ylim([-1,1]);ylabel('verticale [m]') ;grid on; axis equal; title('blu=vero  rosso=stima  verde=ricostruito')
or = [0; 0];                                                               %Definizione della origine
plot_origine = plot(or(1), or(2), 'ob', 'MarkerSize',7);                   %Plot della origine
%Rinomino le variabili stimate e regolarizzate per brevità
x_results1 = squeeze(log_EKF_s.x_hat_correction(1,:,:));
x_results2 = squeeze(log_EKF_s.x_hat_correction(2,:,:));
x_smoothed1=squeeze(log_EKF_s.x_hat_smoothed(1,:,:));
x_smoothed2=squeeze(log_EKF_s.x_hat_smoothed(2,:,:));
N_ = size(log_EKF_s.x_hat_correction, 3);                                  %Quanti campioni plotta
a1=[-1;0];a2=[1;0];a3=[0;-1];a4=[0;1]; plot_lim=plot(a1(1),a1(2),'k',a2(1),a2(2),'k',a3(1),a3(2),'k',a4(1),a4(2),'k'); %servono per tenere fisso il grafico

for k = 1 : N_
    % Grafico la traiettoria vera
    q = [teta1(k), teta2(k)];                                              %definizione delle variabili di giunto
    [p1,p2] = cinematica(or, q(1), q(2), l1, l2);                          %calcolo posizione vera dei giunti
    %grafico i due link ed i due giunti
    plot_link1 = plot([or(1),p1(1)],[or(2),p1(2)],'b');
    plot_giunto2 = plot(p1(1), p1(2), 'ob', 'MarkerSize',7);
    plot_link2 = plot([p1(1), p2(1)], [p1(2), p2(2)], 'b');
    plot_ee = plot(p2(1), p2(2), '*b', 'MarkerSize',10);

    %Inserisco la traiettoria stimata  
    p = [x_results1(k), x_results2(k)];                                    %variabili di giunto stimate
    [o1, o2] = grafico_ekf(or, p(1), p(2), l1, l2);                        %calcolo posizione stimata dei giunti
    %sovrappongo i due link e giunti stimati
    plot_link1_h = plot([or(1), o1(1)],[or(2), o1(2)], 'r--');
    plot_giunto2_h = plot(o1(1), o1(2), 'hr', 'MarkerSize',6);
    plot_link2_h = plot([o1(1), o2(1)], [o1(2), o2(2)], 'r--');
    plot_ee_h = plot(o2(1), o2(2), 'hr', 'MarkerSize',6);

    %Inserrisco la traiettoria regolarizzata
    p = [x_smoothed1(k), x_smoothed2(k)];                                  %variabili di giunto regolarizzate
    [o1, o2] = grafico_ekf(or, p(1), p(2), l1, l2);                        %calcolo posizione regolarizzata dei giunti
    %sovrappongo i due link e giunti regolarizzati
    plot_link1_s = plot([or(1), o1(1)],[or(2), o1(2)], 'g-.');
    plot_giunto2_s= plot(o1(1), o1(2), 'pg', 'MarkerSize',5);
    plot_link2_s = plot([o1(1), o2(1)], [o1(2), o2(2)], 'g-.');
    plot_ee_s = plot(o2(1), o2(2), 'pg', 'MarkerSize',5);
    timer=text(-0.9,0.9,"Timer: "+num2str(time(k),2));                     %aggiungo il timer

    pause (0.001);                                                         %dopo 0.001 secondi cancella i grafici attivi e prosrgui sostituendoli con quelli della iterazione successiva
    %cancella i grafici seguenti
    delete(plot_link1);
    delete(plot_giunto2);
    delete(plot_link2);
    delete(plot_ee);
    delete(plot_link1_h);
    delete(plot_giunto2_h);
    delete(plot_link2_h);
    delete(plot_ee_h);
    delete(plot_link1_s);
    delete(plot_giunto2_s);
    delete(plot_link2_s);
    delete(plot_ee_s);
    delete(timer);
end

%===== Index funzioni usate nello script =====%
%Stesse funzioni con identico nome sono state già presentate e commentate
%nello script principale

function [p1,p2]=cinematica(p0, q1, q2,l1,l2)
% q1=teta1;
% q2=teta2;
p1 = [l1*sin(q1); l1*cos(q1)];
p2 = p0 + p1 + [l2*sin(q2); l2*cos(q2)];
end

function [o1,o2]=grafico_ekf(or, x_results1, x_results2, l1, l2)
p_1 = x_results1;
p_2 = x_results2;
o1 = [l1*sin(p_1); l1*cos(p_1)];
o2 = or + o1 + [l2*sin(p_2); l2*cos(p_2)];
end