
function vcvs(nd1,nd2,ni1,ni2,val)
% Adds the stamp of a dependent voltage-controlled 
% voltage-source(VCVS)to the matrices in circuit 
% representation.
%
%   ni1 O-------(+)        |----------o nd1
%                          |
%                         /+\
%                      | /   \    Vnd1-Vnd2 = val*(Vni1-Vni2)
%                Ivcvs | \   /
%                      V  \-/ 
%                          |
%   ni2 O-------(-)        |----------o nd2
%
%  (1) "nd1 & nd2" are the nodes across the dependent
%                  voltage source.
%  (2) "ni1 & ni2" are the nodes corresponding to the 
%                  controller voltage
%
%   nd1: (+) node   \
%   nd2: (-) node   |----->  Vnd1-Vnd2 = val*(Vni1-Vni2)
%   ni1: (+) node   |
%   ni2: (-) node   /
%   A: gain, V/V
%---------------------------------------------------------------

% The body of the function will go here!
global G b C

no1 = nd1;
no2 = nd2;

d = size(G,1);
xr = d+1;

b(xr) = 0;
G(xr,xr) = 0;
C(xr,xr) = 0;


if (no1 ~= 0)
    G(xr,no1) = 1;
    G(no1,xr) = 1;
end
if (no2 ~= 0)
    G(xr,no2) = -1;
    G(no2,xr) = -1;
end
if (ni1 ~= 0 || ni2 ~= 0)
    if (ni1 ~= 0)
        G(xr,ni1) = G(xr,ni1) - val;
    end
    if (ni2 ~= 0)
        G(xr,ni2) = G(xr,ni2) + val;
    end
end


end
