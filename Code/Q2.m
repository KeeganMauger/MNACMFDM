clear all
close all
clc
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
R3 = 123.346641;
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

b1 = b*u_t;
b2 = b*v_t;
b3 = b*w_t;

%--------------------------------------------------------------------------
% FDM Solution: Backwards Euler with Timestep 0.001
%--------------------------------------------------------------------------
h = 0.001;
x = sparse(width(G),numel(t_vec_1));
C_h = C*(1/h);

for j = 1:3;
    for n=1:1000
        if n == 1000
            break;
        end
        BE_LHS = C_h + G;
        BE_RHS = C_h*x(:,n) + b1(:,n+1);
        [L, U, P, Q]= lu( sparse(BE_LHS) , 0.1 );
        Z = L \(P* sparse(BE_RHS));
        Y = U \ Z;
        x(:,n+1) = Q*Y;
    end
    BE1 = x;
    for n=1:1000
        if n == 1000
            break;
        end
        BE_LHS = C_h + G;
        BE_RHS = C_h*x(:,n) + b2(:,n+1);
        [L, U, P, Q]= lu( sparse(BE_LHS) , 0.1 );
        Z = L \(P* sparse(BE_RHS));
        Y = U \ Z;
        x(:,n+1) = Q*Y;
    end
    BE2 = x;
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
end










% w = 0;
% s = j*w;
% 
% A = G + s*C;
% A0 = full(A);
% 
% 
% V0 = linspace(-10,10,21);
% b0 = sparse((width(G)),width(V0));
% for i = 1:width(V0)
%     b0(9,i) = V0(i);
% end
% 
% x = sparse((width(G)),width(V0));
% 
% for j = 1:width(V0)
%     x(:,j) = (G + s*C) \ b0(:,j);
% end
% 
% figure(2)
% plot(V0,x(5,:))
% hold on
% plot(V0,x(3,:))
% grid on
% axis([-10 10 -100 100])
% title('Vout vs Vin')
% xlabel('Vin')
% ylabel('Vout')
% legend('Vout','V3')
% saveas(gcf,'Figure2')
% hold off