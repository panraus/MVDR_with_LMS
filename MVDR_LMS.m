%% MVDR criterion  MVDR beamforming   LMS ALGORITHM
close all;clear all;clc;
source=1;           %  signal number
interference=1;     %  interference number
N=7;               %array number   
theta_s=-10;          %DOA of signal
theta_i=[-60];  %DOA of interference
ss=1024;            %snapshot
snr=[-10 10 50];    %  SNR  
j=sqrt(-1);
%% 
% w=[pi/6 pi/6 pi/3 pi/4]';
% for m=1:4
%     S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*[0:ss-1]);        %3*1024
% end

for m=1:(source+interference)
    S(m,:) = 10.^(snr(m)/10)*(randn(1,ss)+j*randn(1,ss));         %Signal and interference
end
%% 

A_i=exp(-j*pi*(0:N-1)'*sin(theta_i/180*pi));%8*4
A_s=exp(-j*pi*(0:N-1)'*sin(theta_s*pi/180));%8*1 40
A = [A_s A_i(:,1:interference)];
%% 
n=randn(N,ss)+j*randn(N,ss);

%% 
X=A*S+n;

%% 
% R=X*X'/ss;
% [Vec Val]=eig(R);
%% 
Wx = A_s'.*2^10;
% u=1/max(max(Val));
% Wx=[1 0 0 0 0 0 0].*2^14;
u=2^(-31)*2^16;
B0H_B0=[0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0;...
    0 0 1 0 0 0 0 ;...
    0 0 0 1 0 0 0;...
    0 0 0 0 1 0 0;...
    0 0 0 0 0 1 0;...
    0 0 0 0 0 0 1];
% B0=[-1 exp(-j*pi*sin(theta_s)) 0 0 0 0 0;...
%     0 -1 exp(-j*pi*sin(theta_s)) 0 0 0 0;...
%     0 0 -1 exp(-j*pi*sin(theta_s)) 0 0 0 ;...
%     0 0 0 -1 exp(-j*pi*sin(theta_s)) 0 0;...
%     0 0 0 0 -1 exp(-j*pi*sin(theta_s)) 0;...
%     0 0 0 0 0 -1 exp(-j*pi*sin(theta_s))];
% B0H_B0 = B0'*B0;
dataout=zeros(1,ss); 
% dataout(1,1)=Wx*X(:,1);
dataout(1,1)=Wx*X(:,1)./2^14;
for i=1:length(dataout)-1 
%     dataout(1,i)=Wx*data_FPGA(:,i)./2^14;
%         Wx=Wx-u*(data_FPGA(:,i)')*B0H_B0*dataout(1,i)
      Wx=Wx-u*(X(:,i)')*B0H_B0*dataout(1,i);
%     Wx=Wx-u*(X(:,i)')*B0HB0*dataout(1,i);
    dataout(1,i+1)=Wx*X(:,i+1)./2^15;
%     dataout(1,i+1)=Wx*X(:,i+1);
end
%% 
% P=inv(A_s'*inv(R)*A_s);

%%  
% W=P*inv(R)*A_s;
phi=-89:1:90;
a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
F=Wx*a;
%figure();
%plot(phi,F);
G=abs(F).^2./max(abs(F).^2);
% G=abs(F).^2;
G_dB=10*log10(G);
%figure();
%plot(phi,G);legend('N=8,d=lamda/2');
figure();
plot(phi,G_dB,'linewidth',2);legend('N=7,d=lamda/2');
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
grid on;
% axis equal;
% axis([-90 90 -100 10]);