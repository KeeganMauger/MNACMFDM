clearvars -except R3
% close all
% clc
%set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',10)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 0.5);

%--------------------------------------------------------------------------
% Modified Nodal Analysis - Transient Analysis
%--------------------------------------------------------------------------

global G b C
nodes = 6;
outNode = 5;
inNode = 9;
N = 1000; %num points
G = sparse(nodes,nodes);
C = sparse(nodes,nodes);
b = sparse(nodes,1);


t_vec_1 = linspace(0,1,1000);
t_vec_2 = linspace(0+0.06,1-0.06,1000-120);
t = 0;
u_t = 0;
v_t = 0;
w_t = 0;

f = 1/0.03;
A1 = 1;
w = 2*pi*f;

X = t_vec_2;
A2 = 1;
mu = 0.5;
sigma = 0.03;
delay = 0.06;


gaussian = A2*exp(-((X-mu).^2/(2*sigma.^2)));

z = 0;

for k = 1:numel(t_vec_1)
    t = t_vec_1(k);
    %input 1
    if t < 0.03
        u_t(k) = 0;
    else
        u_t(k) = 1;
    end
    %input 2
    v_t(k) = A1*sin(w*t);
    %input 3
    if t < 0.06
        w_t(k) = 0;
    elseif t >= 0.06 && t < 0.94
        w_t(k) = gaussian(k-60);
    elseif t >= 0.94
        w_t(k) = 0;
    else
        w_t(k) = 0;
    end
    

end

% plot(t_vec_1,u_t);
% hold on
% plot(t_vec_1,v_t);
% plot(t_vec_1,w_t);
% hold off






%--------------------------------------------------------------------------
% MNA setup
%--------------------------------------------------------------------------
Vin = 1;
Vprobe = 0;
R1 = 1;
R2 = 2;
%R3 = 123.346641;
% R3 = 50;
R4 = 0.1;
R5 = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;
In = 0.001;
Cn = 0.00001;

cap(1,2,C1);
res(1,2,R1);
res(2,0,R2);
res(3,6,R3);
res(4,5,R4);
res(5,0,R5);
xr = vol(6,0,Vprobe);
ind(2,3,L1);
vol(1,0,Vin);
vcvs(4,0,xr,0,alpha);
cur(4,0,6);
cap(4,0,Cn);

fullC = full(C);
fprintf('\nThe updated C matrix is given by:\n');
disp(fullC);

% b1 = b*u_t;
% b2 = b*v_t;
b3 = zeros(10,1000);
b3(inNode,:) = w_t;

% r = normrnd(mu,std,1000);
for i = 1:1000
    r = In*randn();
     b3(4,i) = r;
     b3(6,i) = -r;
    %b3(4,i) = 1;
end
% make b*wt and In into 2 different things

%--------------------------------------------------------------------------
% FDM Solution: Backwards Euler with Timestep 0.001
%--------------------------------------------------------------------------
h = 0.001;
x = sparse(width(G),numel(t_vec_1));
C_h = C*(1/h);

for n=1:1000
    if n == 1000
        break;
    end
    BE_LHS = C_h + G;
    BE_RHS = C_h*x(:,n) + b3(:,n+1);
    [L, U, P, Q]= lu( sparse(BE_LHS) , 0.1 );
    Z = L \(P* sparse(BE_RHS));
    Y = U \ Z;
    x(:,n+1) = Q*Y;
end
BE3 = x;


Fs = 1/h;
dF = Fs/N;
f = -Fs/2:dF:Fs/2 - dF;

figure(8)
subplot(2,1,1);
plot(t_vec_1, w_t,'LineWidth',1);
hold on;
plot(t_vec_1, BE3(outNode,:),'LineStyle','-.', 'color','r','LineWidth',1);
title('Vin and Vout for Gaussian Function Input')
ylabel('Voltage (V)')
xlabel('Time (s)')
legend('Vin','Vout')
hold off

subplot(2,1,2);
fullBE3 = full(BE3(outNode,:));
FT3 = fft(fullBE3);
plot(f,fftshift(abs(FT3)))
title('Fourier Transform of Vout with Noise')
ylabel('Magnitude')
xlabel('Frequency (Hz)')
saveas(gcf,'Figure8')

