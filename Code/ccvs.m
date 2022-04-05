function ccvs(nd1,nd2,xrE,val)
% Adds the stamp of a dependent current-controlled 
% voltage-source(CCVS)to the matrices in circuit 
% representation.
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
%   ni1 O-------|                  |----------o nd1
%               |                  |
%               |                 /+\
%               |              | /   \    Vnd1-Vnd2 = val*(Ij)
%           Ij  V           Ik | \   /
%               |              V  \-/ 
%               |                  |
%   ni2 O-------|                  |----------o nd2
%
%  (1) "nd1 & nd2" are the nodes across the dependent
%                  voltage source.
%  (2) "ni1 & ni2" are the nodes corresponding to the 
%                  controller current
%
%   nd1: (+) node   \
%   nd2: (-) node   |----->  Vnd1-Vnd2 = val*(Ij)
%   ni1: (+) node   |
%   ni2: (-) node   /
%   F: alpha, ohms
%---------------------------------------------------------------
global G b C 

no1 = nd1;
no2 = nd2;

d = size(G,1);
xr1 = d+1;
xr2 = d+2;

b(xr1) = 0;
b(xr2) = 0;
G(xr1,xr1) = 0;
G(xr2,xr2) = 0;
C(xr1,xr1) = 0;
C(xr2,xr2) = 0;

if (no1 ~= 0)
    G(xr2,no1) = 1;
    G(no1,xr2) = 1;
end
if (no2 ~= 0)
    G(xr2,no2) = -1;
    G(no2,xr2) = -1;
end

if (ni1 ~= 0)
    G(xr1,ni1) = 1;
    G(ni1,xr1) = 1;
end
if (ni2 ~= 0)
    G(xr1,ni2) = -1;
    G(ni2,xr1) = -1;
end
if (ni1 ~= 0 || ni2 ~= 0)
    G(xr2,xr1) = G(xr2,xr1) - val;
end


end

