clear all
close all
clc
%set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',10)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 0.5);

%--------------------------------------------------------------------------
% Finite Difference Method Simulation
%--------------------------------------------------------------------------

Vin = linspace(0.1,10,30);
for J=1:length(Vin)
fprintf('\nJ = %d',J);
% Fixed bottleneck of 0.2x10-7m
% Solving V=V0 @ x=0 and V=0 @ x=L in region LxW
% Implement funtion 'pbaspect' to fix Z aspect ratio
clearvars -except Vin J Ix_total Iavg

Conc = 1e19;
L = 200e-9;
W = 100e-9;
V0 = Vin(J);

fMesh = 1;
nx = fMesh*200;
ny = fMesh*100;
La = linspace(0,L,nx);
Wa = linspace(0,W,ny);
G = sparse(nx*ny);
%V = sparse(nx,ny);
F = sparse(1,nx*ny);


Acond = 1;              % background conductivity of region, low resistance
Bcond = 1e-2;           % Conductivity of boxes, highly resistive
cMap = zeros(nx,ny);
Lb = 40e-9;
Wb = 40e-9;
for u = 1:nx
    for v = 1:ny
        if (u >= 80 && u <= 120)
            if v >= 0 && v <= 40-1+J
                cMap(u,v) = Bcond;
            elseif v >= 60+1-J && v <= 100
                cMap(u,v) = Bcond;
            else
                cMap(u,v) = Acond;
            end
        else
            cMap(u,v) = Acond;
        end
    end
end





for i = 1:nx                %Iteration through length
    for j = 1:ny            %Iteration through width
        n = j + (i-1)*ny;

        if i == 1          % x=0 BCs
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = V0;
        elseif i == nx     % x=1 BCs
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;       %F(n)=0 sets z at final width to 0
            
% COMMENT BELOW FOR 1a
            %F(n) = 1;       %F(n)=1 sets z at final width to 1
            
        elseif j == 1                 % y=0 BCs
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;

        elseif j == ny                % y=1 BCs
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            
            
% COMMENT ABOVE FOR 1a
            
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end
            
    end
end
% figure(1)
% spy(G)

V = G\F';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end


for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap .* Ex;        %Jx
eFlowy = cMap .* Ey;        %Jy
C0 = sum(eFlowx(1, :));
Cnx = sum(eFlowx(nx, :));
Curr = (C0 + Cnx) * 0.5;

Ex = Ex';
Ey = Ey';

% figure(4)
% subplot(1, 2, 2), quiver(Ex, Ey);
% axis([0 nx 0 ny]);
% title('Electric Field Map')
% xlabel('Region Length')
% ylabel('Region Width')
% pbaspect([1 1 0.5])
% subplot(1, 2, 1), H = surf(La,Wa,Vmap');
% set(H, 'linestyle', 'none');
% %view(90, 270)
% title('Voltage Map')
% xlabel('Region Length')
% ylabel('Region Width')
% pbaspect([1 1 0.5])
% saveas(gcf,'Figure4')

% clearvars -except Ex Ey