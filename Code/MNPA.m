clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

% DC Case
global G b C
nodes = 6;
G = sparse(nodes,nodes);
C = sparse(nodes,nodes);
b = sparse(nodes,1);

% MNA setup
Vin = 1;
Vprobe = 0;
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R5 = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;

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

w = 0;
s = j*w;

A = G + s*C;
A0 = full(A);


V0 = linspace(-10,10,21);
b0 = sparse((width(G)),width(V0));
for i = 1:width(V0)
    b0(9,i) = V0(i);
end

x = sparse((width(G)),width(V0));

for j = 1:width(V0)
    x(:,j) = (G + s*C) \ b0(:,j);
end

figure(1)
plot(V0,x(5,:))
hold on
plot(V0,x(3,:))
grid on
axis([-10 10 -100 100])
title('Vout vs Vin')
xlabel('Vin')
ylabel('Vout')
legend('Vout','V3')
saveas(gcf,'Figure1')
hold off




clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
% AC Case
global G b C
nodes = 6;
G = sparse(nodes,nodes);
C = sparse(nodes,nodes);
b = sparse(nodes,1);

% MNA setup
Vin = 1;
Vprobe = 0;
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R5 = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;

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

w = linspace(0,100,100);
s = j*w;

% A = G + s.*C;
% A0 = full(A);

%V0 = linspace(-10,10,21);
% b0 = sparse((width(G)),width(w));
% for i = 1:width(w)
%     b0(6,i) = V0(i);
% end

% x = (G + s*C) \ b;
x = sparse((width(G)),width(w));

for i = 1:width(w)
    x(:,i) = (G + s(i)*C) \ b;
end
    
figure(2)
plot(w,abs(x(5,:)))
hold on
plot(w,abs(x(3,:)))
grid on
title('|Vout| vs Frequency')
xlabel('Frequency')
ylabel('Vout')
saveas(gcf,'Figure2')
hold off

mag = abs(x(5,:)./x(1,:));
gain = mag2db(mag);
figure(3)
plot(w,gain)
grid on
title('Gain vs Frequency')
ylabel('Gain (dB)')
xlabel('Frequency')
saveas(gcf,'Figure3')

% AC Case - Normal Perturbations on C

clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
global G b C
nodes = 6;
G = sparse(nodes,nodes);
C = sparse(nodes,nodes);
b = sparse(nodes,1);

% MNA setup
Vin = 1;
Vprobe = 0;
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R5 = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;

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

w = pi;
s = j*w;
std = 0.05;
mu = C1;
numR = 1000;
r = normrnd(mu,std,1,numR);

C0 = zeros(width(G),length(G),numR);
for i = 1:width(r)
    C(1:2,1:2) = 0;
    cap(1,2,r(i));
    C0(:,:,i) = C;
end

x = sparse((width(G)),numR);
for i = 1:numR
    x(:,i) = (G + s*C0(:,:,i)) \ b;
end

mag = abs(x(5,:)./x(1,:));
gain = mag2db(mag);

figure(4)
subplot(1,2,1);
hist(r)
title('Histogram of C, Normal Distribution')
xlabel('C')
ylabel('Number')
subplot(1,2,2);
hist(mag,10)
title('Histogram of Gain, Normal Distribution')
xlabel('Gain (V/V)')
ylabel('Number')
saveas(gcf,'Figure4')




    
