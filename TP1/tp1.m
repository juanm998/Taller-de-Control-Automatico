%%
clear; close all; clc;

% Configuración del Bode

optionss=bodeoptions; 

optionss.MagVisible='on';
optionss.PhaseMatching='on';
optionss.PhaseMatchingValue=-180;
optionss.PhaseMatchingFreq=100;
optionss.Grid='on';

%%
%Defino las funciones de estado no lineales:

syms x1 x2 u y;
x = [x1;x2];
f1 = x2;
f2 = 10-(u/(0.4*x1-0.01*x2));
y=x1;
f = [f1; f2];

%Defino los puntos de equilibrio:
x1e = 1;
x2e = 0;
ue = 4; 

ye = x1e;


% Defino los jacobianos:
A = jacobian(f,x);
B = jacobian(f,u);
C = jacobian(y,x);
D = jacobian(y,u);

% Evaluar el jacobiano en los valores de equilibrio:
A = subs(A, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
B = subs(B, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
C = subs(C, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
D = subs(D, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});

Ass = double(A);
Bss = double(B);
Css = double(C);
Dss = double(D);

% Transferencia linealizada:
P_lineal = zpk(ss(Ass, Bss, Css, Dss))



%%
%Elijo una frecuencia tal que el cero de fase no minima este lo "suficientemente"
%lejos de los polos de fase no minima para tener un rango amplio de
%frecuencia donde poder estabilizar
wgc=100;

%Con esto calculo el tiempo de muestreo definiendo 5 grados al
%digitalizador:
T = (4/wgc)*tand(5/2)

%Este parametro T me define los polos y ceros de la transferencia de el Padé
p_h=4/T;
H=zpk([p_h],[-p_h],-1); 


P = minreal(H*P_lineal);

%%
%Defino las transferencias "pasa todo" y "fase minima" de la nueva planta P: 
Pap=zpk([p_h -3.04],[3.04 -p_h],1)
Pmp=zpk([],[-3.29 -3.04],2.5)

%Grafico los diagramas de Bode de las transferencias:
figure();
hold on;
bode(Pap, optionss, {0.1, 10000});
bode(Pmp, optionss, {0.1, 10000});
hold off;
 

legend('Pap', 'Pmp');  
title('PAP PMP');
set(findall(gcf,'type','line'),'linewidth',2);

%Defino el compensador, en principio le pongo una ganancia unitaria luego
%la ajusto para lograr que la frecuencia wgc elegida sea la de cruce,
%tambien agrego un polo de fase minima de un modulo alto para que sea la
%transferencia del controlador no sea impropia.

%Le pongo ganancia -1 al controlador porque necesito que retrase fase

C=zpk([-3.04 -3.29],[0 -100000],-1)
%%
%Defino la transferenica L para poder ajustar la ganancia con el Bode
Pmpc=minreal(Pmp*C);
L=minreal(Pmpc*Pap);

%grafico del Bode de P,C y L
figure();
hold on;
bode(P_lineal, optionss, {0.1, 10000});
optionss.PhaseMatchingValue=90;
bode(C, optionss, {0.1, 10000});
optionss.PhaseMatchingValue=-180;
bode(L,  optionss, {0.1, 10000});

hold off;

legend('P', 'C', 'L');  
title('Transferencia de lazo abierto');
set(findall(gcf,'type','line'),'linewidth',2);
%%
%Ajusto la ganancia
k=db2mag(132)
L_final=k*L;

figure();
hold on;
bode(L_final, optionss, {0.1, 10000});
hold off;

legend('L_final');  
title('Transferencia de lazo abierto');
set(findall(gcf,'type','line'),'linewidth',2);
%%Con esto tengo un margen de fase igual a 90 grados y wgc=100
%%
%Grupo de las 4

S = minreal(1 / (1+L_final));

T = minreal(L*S);

cS = minreal(C*S);

pS = minreal(P*S);

%Todos los autovalores-> parte real negativa-> Grupo de las 4 estables
eig(S)
eig(T)
eig(cS)
eig(pS)

%%
% Paso el controlador a discreto
C_digital = c2d(C, T, 'tustin');
