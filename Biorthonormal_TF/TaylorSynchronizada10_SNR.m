clear all; clear memory; close all; clc;
%Modelo de orden K=3
a1=20;a2=15;a3=10;a4=5;
sig1=-0.02;sig2=-0.05;sig3=-0.09;sig4=-0.08;
f1=0.20;f2=0.5; f3=0.9; f4=0.7;
theta1=pi/6; theta2=0;theta3=0;theta4=pi/6;
%Para f1=0.2, T1=5, and 4T1=20s, 20s de -10 a 10s
Fs=10;Ts=1/Fs; %  20*Fs=200
%Fs=30;Ts=1/Fs; %  20*Fs=600
t=(-10:Ts:10-Ts)';
T=[ones(20*Fs,1) t (t.^2)/2 (t.^3)/6];
B=[T T.*exp(2j*pi*f1*t) T.*exp(2j*pi*f2*t) T.*exp(2j*pi*f3*t) T.*exp(2j*pi*f4*t) T.*exp(-2j*pi*f4*t) T.*exp(-2j*pi*f3*t) T.*exp(-2j*pi*f2*t) T.*exp(-2j*pi*f1*t)];
Bp=(B' * B) \ B'; %or Bp=pinv(B) the full Moore-Penrose pseudoinverse  valid for any shape of X;
%Bp=inv(B'*B)*B';
Tp=(T' * T) \ T';  %Tp=inv(T'*T)*T'; the problem is that the inv function can be unstable.
% Best code for alpha in LMS solution: alpha = X \ x;
Tph=Tp';
Bph=Bp';
H0=[Bph(:,5),Bph(:,9), Bph(:,13), Bph(:,17)];
%H1=[Bph(:,6),Bph(:,10),Bph(:,14), Bph(:,18)];
FH0=fftshift(fft(H0,2^12));
%FH1=fftshift(fft(H1,2^12));
f=(-.5:2^-12:.5-2^-12)*Fs;  
%%%%%%%%%%%%%%%%%%%%%% Taylor-Fourier Filters
figure(1)
subplot(4,1,1)
plot(f(1024:3072),abs(FH0(1024:3072,3)),'LineWidth', 2)
xlim([-2.5 2.5])
ylabel('$H_1(f)$','interpreter','latex')
title('Frequency Response of the Adaptive Taylor-Fourier Filters')
subplot(4,1,2)
plot(f(1024:3072),abs(FH0(1024:3072,4)),'LineWidth', 2)
ylabel('$H_2(f)$','interpreter','latex')
xlim([-2.5 2.5])
subplot(4,1,3)
plot(f(1024:3072),abs(FH0(1024:3072,1)),'LineWidth', 2)
ylabel('$H_3(f)$','interpreter','latex')
xlim([-2.5 2.5])
subplot(4,1,4)
plot(f(1024:3072),abs(FH0(1024:3072,2)),'LineWidth', 2)
ylabel('$H_4(f)$','interpreter','latex')
xlim([-2.5 2.5])
xlabel('Frequency f (Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Signal Composition
t=0:Ts:120;
s1=a1*exp(sig1*t).*cos(2*pi*f1*t+theta1);
s2=a2*exp(sig2*t).*cos(2*pi*f2*t+theta2);
s3=a3*exp(sig3*t).*cos(2*pi*f3*t+theta3);
s4=a4*(1+t/6).*exp(sig4*t).*cos(2*pi*f4*t+theta4);
s=s1+s2+s3+s4;s=s';
%% variance
st3 =  s;
Eoc = sum(st3.^2);
st4 = st3/sqrt(Eoc);
snr = 10;
for kkk = 1:23
    variance = 10^(-snr/10);
    noise = sqrt(variance)*randn(size(st4));
    st = st4 + noise;
    ce = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimacion de la envolvente
E0=convn(H0,st);
% ss_r=deconv(E0(101:1101,1),s);
[ss_r1 rr1]=deconv(E0(:,1),H0(:,1));
[ss_r2 rr2]=deconv(E0(:,2),H0(:,2));
[ss_r3 rr3]=deconv(E0(:,3),H0(:,3));
[ss_r4 rr4]=deconv(E0(:,4),H0(:,4));

ss_rr=ss_r1+ss_r2+ss_r3+ss_r4;
ss_r=ss_rr(1:901);
ss_s=3.9289*st(1:901);

figure (100)
plot(ss_s)
hold on
plot(1:901,real(ss_r(1:901)))
%hold on

% pause 

error = ss_s - ss_r;

 %error2 = st.^2 - sk'.^2;
    Eerror = sum(error.^2);
    %Eerror2 = sum(error2.^2);
    vari(kkk) = Eerror/length(ss_s(1:901));
    %vari2(kkk) = Eerror2;
    snr = snr + 5;
end

snr = 10:5:120;
figure (200)
% b1 = semilogy(snr,vari,'k');
plot(snr,vari,'k');
title('Output variance vs SNR ', 'FontSize', 24)
ylabel('E = st - sk,sum(E2)', 'FontSize', 24)
xlabel('SNR', 'FontSize', 24)
grid on


e1=2*abs(E0(101:1101,1));
e2=2*abs(E0(101:1101,2));
e3=2*abs(E0(101:1101,3));
e4=2*abs(E0(101:1101,4));
figure(2)
subplot(1,4,1)
plot(t(1:1001),s1(1:1001),t(1:1001),e1,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 1')
ylabel('Signal and Envelope')
subplot(1,4,2)
plot(t(1:1001),s2(1:1001),t(1:1001),e2,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 2')
subplot(1,4,3)
plot(t(1:1001),s3(1:1001),t(1:1001),e3,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 3')
subplot(1,4,4)
plot(t(1:1001),s4(1:1001),t(1:1001),e4,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimacion de atenuaci??n 0.0200005, 0.0500073, 0.0901215
env=2*abs(E0(201:501,:));
est0=convn(Tp(1,:)',env);
est1=-convn(Tp(2,:)',env);
est2=convn(Tp(3,:)',env);
div=est1./est0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tao0=est0(:,4); tao1=est1(:,4); tao2=est2(:,4);
c0=tao0;
alpha=(tao1-sqrt(tao1.^2 - 2*tao0.*tao2))./tao0;
c1=tao1-tao0.*alpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(t(200:300),div(200:300,:),'LineWidth', 2)
xlabel('Time (s)')
title('Estimated Exponential Attenuation ')
ylabel('Attenuation $\sigma_n$','interpreter','latex')
xlim([20 30])
%la predicci??n es excelente.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
subplot(3,1,1)
plot(t(200:300),tao0(200:300),'LineWidth', 2)
xlim([20 30])
subplot(3,1,2)
plot(t(200:300),c1(200:300),'LineWidth', 2)
xlim([20 30])
subplot(3,1,3)
plot(t(200:300),alpha(200:300),'LineWidth', 2 )
xlim([20 30])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Approximation Root Mean Squared Error per window
nv=1;
for k=1:200:1200
    sigma(:,nv)=Bp*s(k:k+199);
    sh=B*sigma(:,nv);
    e(nv)=sqrt((s(k:k+199)-sh).'*(s(k:k+199)-sh))/sqrt(200);
    nv=nv+1;
end;
figure(5)
stem(1:nv-1,abs(e))
rmse=mean(abs(e))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrum
%figure(6)
%subplot(2,1,1)
%plot(t,s)
%subplot(2,1,2)
%s(2^12)=0;
%S=fftshift(abs(fft(s)));
%f=(-0.5:2^-12:0.5-2^-12)'*Fs;
%plot(f, S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse response Analysis
h=real(H0);
tf=(-10:Ts:10-Ts)';
%figure(8)
%plot(1:200,h(:,2),'LineWidth', 2)
%xlabel('Time (s)')
%ylabel('$h_2(t)$','interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
subplot(4,1,1)
plot(tf,h(:,1),'LineWidth', 2)
title('Filter Impulse Responses')
ylabel('$h_1(t)$','interpreter','latex')
subplot(4,1,2)
plot(tf,h(:,2),'LineWidth', 2)
ylabel('$h_2(t)$','interpreter','latex')
subplot(4,1,3)
plot(tf,h(:,3),'LineWidth', 2)
ylabel('$h_3(t)$','interpreter','latex')
subplot(4,1,4)
plot(tf,h(:,4),'LineWidth', 2)
xlabel('Time (s)')
ylabel('$h_4(t)$','interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=convn(h,s);
figure(9)
plot(t,r(1:1201,4),t,s4)
H=fftshift(fft(h,2^12));
figure(10)
subplot(4,1,1)
plot(f(1024:3072),abs(H(1024:3072,3)),'LineWidth', 2)
xlim([-2.5 2.5])
ylabel('$H_1(f)$','interpreter','latex')
title('Frequency Response of the Adaptive Taylor-Fourier Filters')
subplot(4,1,2)
plot(f(1024:3072),abs(H(1024:3072,4)),'LineWidth', 2)
ylabel('$H_2(f)$','interpreter','latex')
xlim([-2.5 2.5])
subplot(4,1,3)
plot(f(1024:3072),abs(H(1024:3072,1)),'LineWidth', 2)
ylabel('$H_3(f)$','interpreter','latex')
xlim([-2.5 2.5])
subplot(4,1,4)
plot(f(1024:3072),abs(H(1024:3072,2)),'LineWidth', 2)
ylabel('$H_4(f)$','interpreter','latex')
xlim([-2.5 2.5])
xlabel('Frequency f (Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of the Gram matrix and normalized for angle calculation
B=[s1',s2',s3',s4'];
G=B'*B
norm=sqrt(diag(G))
Gn= [1, G(1,2)/(sqrt(G(1,1))*sqrt(G(2,2))), G(1,3)/(sqrt(G(1,1))*sqrt(G(3,3))), G(1,4)/(sqrt(G(1,1))*sqrt(G(4,4)));...
    G(2,1)/(sqrt(G(2,2))*sqrt(G(1,1))), 1, G(2,3)/(sqrt(G(2,2))*sqrt(G(3,3))), G(2,4)/(sqrt(G(2,2))*sqrt(G(4,4)));...
    G(3,1)/(sqrt(G(3,3))*sqrt(G(1,1))), G(3,2)/(sqrt(G(3,3))*sqrt(G(2,2))), 1, G(3,4)/(sqrt(G(3,3))*sqrt(G(4,4)));...
    G(4,1)/(sqrt(G(4,4))*sqrt(G(1,1))), G(4,2)/(sqrt(G(4,4))*sqrt(G(2,2))), G(4,3)/(sqrt(G(4,4))*sqrt(G(3,3))),1];
angles=acos(Gn)*180/pi
%Total Phasor Error (TPE)
%s1=a1*exp(sig1*t).*cos(2*pi*f1*t+theta1);
phi1=a1*exp(sig1*t).*exp(1i*(theta1))/2; phi1=phi1.';
phi2=a2*exp(sig2*t).*exp(1i*(theta2))/2; phi2=phi2.';
phi3=a3*exp(sig3*t).*exp(1i*(theta3))/2; phi3=phi3.';
phi4=a4*(1+t/6).*exp(sig4*t).*exp(1i*(theta4))/2; phi4=phi4.';
%%%%%%%
phig1=E0(101:1301,1).*exp(-2i*pi*f1*t');
phig2=E0(101:1301,2).*exp(-2i*pi*f2*t');
phig3=E0(101:1301,3).*exp(-2i*pi*f3*t');
phig4=E0(101:1301,4).*exp(-2i*pi*f4*t');
%error
e1=phi1-phig1;
e2=phi2-phig2;
e3=phi3-phig3;
e4=phi4-phig4;
TPE1= sqrt(e1.*conj(e1)./(phi1.*conj(phi1)));
TPE2= sqrt(e2.*conj(e2)./(phi2.*conj(phi2)));
TPE3= sqrt(e3.*conj(e3)./(phi3.*conj(phi3)));
TPE4= sqrt(e4.*conj(e4)./(phi4.*conj(phi4)));
figure(11)
subplot(4,1,1)
plot(t,TPE1)
subplot(4,1,2)
plot(t,TPE2)
subplot(4,1,3)
plot(t,TPE3)
subplot(4,1,4)
plot(t,TPE4)
%%%%%%%%%%%%%%%
TPE=[TPE1(201:901), TPE2(201:901),TPE3(201:901), TPE4(201:901)];
max(TPE)*100
