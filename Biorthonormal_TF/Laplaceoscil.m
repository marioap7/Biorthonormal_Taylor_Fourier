syms s0 A t sigma f0 theta s f;
%Pi=sym(pi,"d");
s0=A*exp(sigma*t)*cos(2*pi*f0*t+theta);
S0=laplace(s0,t,s);
s3=A*(1+t/6)*exp(sigma*t)*cos(2*pi*f0*t+theta);
S3=laplace(s3,t,s);
%%%%%%%%%%%%%%%%%%%Evaluacion
A=20;sigma=-0.02; f0=0.20;theta=pi/6;
sc0=subs(s0);
Sc0=subs(S0);
Sc0f=subs(Sc0, s, 2i*pi*f);
A=15;sigma=-0.05; f0=0.50; theta=0;
sc1=subs(s0);
Sc1=subs(S0);
Sc1f=subs(Sc1, s, 2i*pi*f);
A=10;sigma=-0.09; f0=0.9;theta=0;
sc2=subs(s0);
Sc2=subs(S0);
Sc2f=subs(Sc2, s, 2i*pi*f);
A=5;sigma=-0.08; f0=0.70;theta=pi/6;
sc3=subs(s3);
Sc3=subs(S3);
Sc3f=subs(Sc3, s, 2i*pi*f);
%%%%%%%%%%%%%%Plots
figure(1)
fplot(sc0,[0,100]);
hold on;
fplot(sc1,[0,100]);
hold on;
fplot(sc2,[0,100]);
hold on;
fplot(sc3,[0,100]);
hold off;
figure(2)
fplot(sc0+sc1+sc2+sc3,[0,100],'LineWidth',1);
title('Oscillation')
xlabel('Time (s)') 
ylabel('Amplitude')
grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
fplot(abs(Sc0f), [0,1]);
hold on;
fplot(abs(Sc1f), [0,1]);
hold on;
fplot(abs(Sc2f), [0,1]);
hold on;
fplot(abs(Sc3f), [0,1]);
hold off;
figure(4)
fplot(abs(Sc2f)+abs(Sc1f)+abs(Sc0f)+abs(Sc3f), [0,1],'LineWidth',1);
title('Oscillation Spectrum')
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(1,2,1)
fplot(sc0+sc1+sc2+sc3,[0,100],'LineWidth',1);
title('Oscillation')
xlabel('Time (s)') 
ylabel('Amplitude')
grid
subplot(1,2,2)
fplot(abs(Sc2f)+abs(Sc1f)+abs(Sc0f)+abs(Sc3f), [0,1],'LineWidth',1);
title('Spectrum')
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
grid
%%%%%%%%%%%%%%%%%
