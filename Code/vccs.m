function vccs(no1, no2, ni1, ni2, val)
% Adds the stamp of a dependent voltage-controlled 
% current-source (VCCS)to the G matrix in circuit 
% representation.
%
%   ni1 O-------(+)   |----------o no1
%                     |
%                    / \
%                   / | \  val*(Vni1 - Vni2)
%                   \ V /
%                    \ /    
%                     |
%   ni2 O-------(-)   |----------o no2
%
%  (1) "no1 & no2" are the nodes across the dependent
%                  current source.
%  (2) "ni1 & ni2" are the nodes corresponding to the 
%                  controller voltage
%
%   no1:  The node at the tail of the current arrow
%   no2:  "   "   "   "  head  "    "   "      "   
%   ni1: (+) node
%   ni2: (-) node
%   val: conductance, seimans
%---------------------------------------------------------------

% The body of the function will go here!

global G

if (nil ~= 0)
    if(no1 ~= 0)
        G(n01,ni1) = val;
    end
    if(no2 ~= 0)
        G(no2,ni1) = -val;
    end
end
if (ni2 ~= 0)
    if(no1 ~= 0)
        G(no2,ni2) = -val;
    end
    if(no2 ~= 0)
        G(no2,ni1) = val;
    end
end


end