%% Calcolo matrici linearizzato

syms u w1 w2 w3 w4  x1 x2 x3 x4 v1 v2;                     % u=tau_m1 w1,w2= dist.integraz w3=tau_d1 w4=tau_d2 x1=teta1 x2=teta2 x3=dteta1 x4=dteta4 v1=v_d v2=v_alfa

dx = [x3; x4];
w = [w1; w2; w3; w4];

G = [1; 0];
T = [1 -1; 0 1];

Ma = [        I1+m1*h1^2+m2*l1^2, l1*m2*h2*cos(x1-x2); 
      l1*m2*h2*cos(x1-x2),               I2+m2*h2^2];

N = [c1+c2, -c2;
       -c2,  c2];

q = [l1*m2*h2*sin(x1-x2)*x4^2-(m1*h1+m2*l1)*g*sin(x1);
           -l1*m2*h2*sin(x1-x2)*x3^2-m2*h2*g*sin(x2)];

a = I1+m1*h1^2+m2*l1^2;
b = l1*m2*h2;
c = I2+m2*h2^2;

detM = 1/(a*c-b^2*cos(x1-x2)^2);

M1 = c*detM;
M2 = -b*cos(x1-x2)*detM;
M3 = a*detM;

N1 = c1+c2;
N2 = c2;

q1 = q(1);
q2 = q(2);

%dinamica inversa (tc)
ddteta1_tc = -M1*(N1*x3-N2*x4)-M2*(N2*x4-N2*x3)-M1*q1-M2*q2+u*M1+M1*(w3-w4)+M2*w4;
ddteta2_tc = -M2*(N1*x3-N2*x4)-M3*(N2*x4-N2*x3)-M2*q1-M3*q2+u*M2+M2*(w3-w4)+M3*w4;

f_tc = [x3+w1; x4+w2; ddteta1_tc; ddteta2_tc]

%modello sistema td

f = [x1+dt*(x3+w1); 
    x2+dt*(x4+w2); 
    x3+dt*ddteta1_tc;
    x4+dt*ddteta2_tc];

h = [D-l1*sin(x1)-l2*sin(x2)+v1;
    pi+x1-x2+v2];

F_primo = [diff(f(1),x1) diff(f(1),x2) diff(f(1),x3) diff(f(1),x4);
    diff(f(2),x1) diff(f(2),x2) diff(f(2),x3) diff(f(2),x4);
    diff(f(3),x1) diff(f(3),x2) diff(f(3),x3) diff(f(3),x4);
    diff(f(4),x1) diff(f(4),x2) diff(f(4),x3) diff(f(4),x4)];

F_x = subs(F_primo, [w1, w2, w3, w4], [0,0,0,0])

D_primo = [diff(f(1),w1) diff(f(1),w2) diff(f(1),w3) diff(f(1),w4);
    diff(f(2),w1) diff(f(2),w2) diff(f(2),w3) diff(f(2),w4);
    diff(f(3),w1) diff(f(3),w2) diff(f(3),w3) diff(f(3),w4);
    diff(f(4),w1) diff(f(4),w2) diff(f(4),w3) diff(f(4),w4)];

D_w = subs(D_primo, [w1, w2 w3 w4], [0,0,0,0])

H_primo = [diff(h(1),x1) diff(h(1),x2) diff(h(1),x3) diff(h(1),x4);
    diff(h(2),x1) diff(h(2),x2) diff(h(2),x3) diff(h(2),x4)];

H_x = subs(H_primo, [v1, v2], [0,0])

M_primo = [diff(h(1),v1) diff(h(1),v2);
    diff(h(2),v1) diff(h(2),v2)];

M_v = subs(M_primo, [v1, v2], [0,0])
