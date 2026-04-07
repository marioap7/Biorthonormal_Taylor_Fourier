data = readmatrix('mex_28dic_1.csv');
ts=data(:,1);
s=data(:,2:6); %Excluyo el 5o y 6o porque espectralmente es bastante disperso
sm=mean(s);
sc=s-sm;      %Señal centrada
figure(1)
plot(ts, s,'LineWidth', 1)
grid on;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Frequency fluctuations')
% Señal centrada o de media nula
SC=abs(fft(sc, 2^12));
f=(0:1/2^12:1-1/2^12)'*10;
figure(2)
plot(f(1:2^11), SC(1:2^11,:), 'LineWidth', 1)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Oscillation Spectra and Filter FR')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0=0;f1=0.315;f2=0.60;%f3=1.9;
%fb=0.28 
Ts=0.1;Fs=1/Ts;Nh=floor(4*10/0.315/2)/10;%Nh=floor(4/0.32)/2;
t=(-Nh:Ts:Nh-Ts)';
T=[ones(2*Nh*Fs,1) t (t.^2)/2 (t.^3)/6];
B=[T T.*exp(2j*pi*f1*t) T.*exp(2j*pi*f2*t)  T.*exp(-2j*pi*f2*t) T.*exp(-2j*pi*f1*t)];
Bp=(B' * B) \ B'; %or Bp=pinv(B) the full Moore-Penrose pseudoinverse  valid for any shape of X;
                  %Bp=inv(B'*B)*B';

Tp=(T' * T) \ T';  %Tp=inv(T'*T)*T'; the problem is that the inv function can be unstable.
                   % Best code for alpha in LMS solution: alpha = X \ x;
P=B*Bp;    %Projection matrix
Tpt=Tp.';
Bpt=Bp.';
H0=[Bpt(:,1),Bpt(:,5),Bpt(:,9)];   %Amplitude and Phase
H1=[Bpt(:,2),Bpt(:,6),Bpt(:,10)];  %First Differentiators
H2=[Bpt(:,3),Bpt(:,7),Bpt(:,11)];  %Second Differentiators
FH0=fft(conj(H0),2^12);
%FH1=fftshift(fft(H1,2^12));
f=(0:2^-12:1-2^-12)*Fs;  
%%%%%%%%%%%%%%%%%%%%%% Taylor-Fourier Filters
figure(3)
subplot(3,1,1)
plot(f(1:2^11),abs(FH0(1:2^11,1)),'LineWidth', 2)
%xlim([-2.5 2.5])
ylabel('$H_1(f)$','interpreter','latex')
title('Frequency Response of the Taylor-Fourier Filters')
subplot(3,1,2)
plot(f(1:2^11),abs(FH0(1:2^11,2)),'LineWidth', 2)
ylabel('$H_2(f)$','interpreter','latex')
%xlim([-2.5 2.5])
subplot(3,1,3)
plot(f(1:2^11),abs(FH0(1:2^11,3)),'LineWidth', 2)
ylabel('$H_3(f)$','interpreter','latex')
%xlim([-2.5 2.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,1,1)
plot(ts, s(:,5),'LineWidth', 2)
grid on;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Frequency fluctuation')
subplot(2,1,2)
plot(f(1:2^11), SC(1:2^11,5),f(1:2^11), abs(FH0(1:2^11,:)),'LineWidth', 2)
xlim([0 1])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on;
title('Oscillation Spectrum and Filter Frequency Response')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse response Analysis
h=2*real(H0);
tf=(-Nh:Ts:Nh-Ts)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(3,1,1)
plot(tf,h(:,1),'LineWidth', 2); grid on;
title('Filter Impulse Responses')
ylabel('$h_1(t)$','interpreter','latex')
subplot(3,1,2)
plot(tf,h(:,2),'LineWidth', 2);grid on;
ylabel('$h_2(t)$','interpreter','latex')
subplot(3,1,3)
plot(tf,h(:,3),'LineWidth', 2);grid on;
ylabel('$h_3(t)$','interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Estimación de la envolvente de los modos
ns=5; %Escoger señal a analizar de 1 al 5 (el 2 o el 5)
x=s(:,ns)-sm(ns);
t=data(:,1);
H0=flip(H0,1);%by fliping the impulse responses the convolution becomes correlation or dot product.
E0=convn(x,H0,'full');
%
%x=s(:,ns);
%y0=real(E0(120:320,1)); 
y1=2*real(E0(63:263,2)); %120:320
y2=2*real(E0(63:263,3));
%Estimación de envolvente y fase
psi0=real(E0(63:263,1));
psi1=E0(63:263,2).*exp(-2j*pi*f1*t);    %Antirotation
psi2=E0(63:263,3).*exp(-2j*pi*f2*t);
%psi3=E0(120:320,4).*exp(-2j*pi*f3*t);
a0=psi0;         phi0=angle(psi0);     %zeros(201,1);  
a1=2*abs(psi1);  phi1=angle(psi1);
a2=2*abs(psi2);  phi2=angle(psi2);
%a3=2*abs(psi3);  phi3=angle(psi3);
figure(6)
subplot(1,3,1)
plot(t,a0,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('         0 (dc)')
ylabel('Mode Signal (blue) and its Envelope (red)')
subplot(1,3,2)
plot(t,y1,t,a1,'LineWidth', 2);grid on;
xlabel('Time (s)')
title(' M 1')
subplot(1,3,3)
plot(t,y2,t,a2,'LineWidth', 2);grid on;
xlabel('Time (s)')
title(' M 2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
subplot(2,3,1)
plot(t,a0,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('         0 (dc)')
ylabel('Mode Signal (blue) and its Envelope (red)')
subplot(2,3,2)
plot(t,y1,t,a1,'LineWidth', 2);grid on;
xlabel('Time (s)')
title(' M 1')
subplot(2,3,3)
plot(t,y2,t,a2,'LineWidth', 2);grid on;
xlabel('Time (s)')
title(' M 2')
subplot(2,3,4)
plot(t,phi0, 'LineWidth', 2);grid on;
xlabel('Time (s)')
title(' Mode 0 (dc)')
ylabel('Mode Phase (rads)')
subplot(2,3,5)
plot(t,unwrap(phi1),'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 1')
subplot(2,3,6)
plot(t,unwrap(phi2),'LineWidth', 2);grid on;
xlabel('Time (s)')
title('Mode 2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First derivative estimates
%FH1=fft(H1,2^12);
%plot(f(1:2^11),abs(FH1(1:2^11,1)),'LineWidth', 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H1=flip(H1,1);                %by flipping the impulse responses the convolution becomes correlation or dot product.
Xid1=convn(H1,x,'full');             %correlation
%Antirotación de la derivada
psid0=real(Xid1(63:263,1));                              %120:320
psid1=Xid1(63:263,2).*exp(-2j*pi*f1*t) - 1j*(2*pi*f1).*psi1;
psid2=Xid1(63:263,3).*exp(-2j*pi*f2*t) - 1j*(2*pi*f2).*psi2;
ad0=psid0; %.*exp(-1j*phi0)
ad1=real(2*psid1.*exp(-1j*phi1));
ad2=real(2*psid2.*exp(-1j*phi2));
phid0=zeros(201,1);
phid1=imag(2*psid1.*exp(-1j*phi1))./a1;
phid2=imag(2*psid2.*exp(-1j*phi2))./a2;
figure(8)
subplot(2,3,1)
plot(t,ad0,'LineWidth', 2);grid on;
%xlabel('Time (s)')
title('       M0')
%title('      0 (dc)')
ylabel('Amplitude 1st Deriv')
subplot(2,3,2)
plot(t,ad1,'LineWidth', 2);grid on;
%xlabel('Time (s)')
title('    M1')
subplot(2,3,3)
plot(t,ad2,'LineWidth', 2);grid on;
%xlabel('Time (s)')
title('    M2')
subplot(2,3,4)
plot(t,phid0/(2*pi)+59.779,'LineWidth', 2);grid on;  %%%%%%%%%%%%%%%%%%%%%%%%%% +59.779 
xlabel('Time (s)')
%
ylabel('Frequency (Hz)')
subplot(2,3,5)
plot(t,phid1/(2*pi)+0.315*2,'LineWidth', 2);grid on;  %%%%%%%%%%%%%%%%%%%%%%%%%% +0.315*2 derivatives + or - by resting 0.315*2
xlabel('Time (s)')
%title('Mode 1')
subplot(2,3,6)
plot(t,phid2/(2*pi)+0.60*2,'LineWidth', 2);grid on;  %%%%%%%%%%%%%%%%%%%%%%%%%%% +0.60*2 derivatives + or - by resting 0.60*2
xlabel('Time (s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second derivative estimates
H2=flip(H2);            %by fliping the impulse responses the convolution becomes correlation or dot product.
Xidd1=convn(H2,x);
%Xidd1=convn(H2,s(:,ns));
%Antirotación de la derivada
psidd0=real(Xidd1(63:263,1));                      %120:320
psidd1=Xidd1(63:263,2).*exp(-2j*pi*f1*t)- 2j*(2*pi*f1).*psid1 + (2*pi*f1).^2.*psi1;
psidd2=Xidd1(63:263,3).*exp(-2j*pi*f2*t)- 2j*(2*pi*f2).*psid2 + (2*pi*f2).^2.*psi2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Second Derivative Amplitudes
add0=real(psidd0); %.*exp(-1j*phi0)
add1=real(2*psidd1.*exp(-1j*phi1))+ a1.*phid1.^2;
add2=real(2*psidd2.*exp(-1j*phi2))+ a2.*phid2.^2;
%
phidd0=zeros(201,1);
phidd1=imag(2*psidd1.*exp(-1j*phi1))-2*ad1.*phid1./a1;
phidd2=imag(2*psidd2.*exp(-1j*phi2))-2*ad2.*phid2./a2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
subplot(2,3,1)
plot(t,add0,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('        M0')
ylabel('Amplitude 2nd Deriv')
subplot(2,3,2)
plot(t,add1,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('    M1')
subplot(2,3,3)
plot(t,add2,'LineWidth', 2);grid on;
xlabel('Time (s)')
title('    M2')
%
subplot(2,3,4)
plot(t,phidd0/(2*pi),'LineWidth', 2);grid on;
xlabel('Time (s)')
ylabel('Phase 2nd Deriv (ROCOF) Hz/s')
subplot(2,3,5)
plot(t,phidd1/(2*pi),'LineWidth', 2);grid on;
xlabel('Time (s)')
subplot(2,3,6)
plot(phidd2/(2*pi),'LineWidth', 2);grid on;
xlabel('Time (s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
plot(ts,s(:,5),ts,a0+y1+y2+59.7786,'LineWidth', 2);grid on;
xlabel('Time (s)')
ylabel('Frequency')
title('Real and Estimated Power Oscillation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x+sm(ns);
Kf=20*Nh-1;
for k=1:201-Kf
    sh=P*(x(k:k+Kf));
    Eh(k)=sh'*sh;
    Es(k)=x(k:k+Kf)'*x(k:k+Kf);
    Ee(k)=Es(k)-Eh(k);
end;
en=sqrt(Ee./Es)*100;
figure(13)
plot(t(1:62),en(1:62),'LineWidth', 2)
xlabel('Time (s)')
title('Normalized Orthogonal Error')
ylabel('Normalized Error (%)')

