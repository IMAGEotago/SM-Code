function [gx, dgdx, dgdphi] = g_sigmoid(x, phi, ~, ~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sigmoid Observation Function                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% Author: Sophie Morris                                                   %
% Created: 14/11/2025                                                     %
% ----------------------------------------------------------------------- %
% Description: This script contains the sigmoid observation function used %
% in Sophie M's Masters' project. A timestamped anaysis plan is available %
% on (https://github.com/IMAGEotago).                                     %
% ----------------------------------------------------------------------- %
% Inputs:                                                                 %
%   x:      hidden state                                                  %
%   phi:    inverse temperature parameter                                 %
% Outputs:                                                                %
%   gx:     predicted response given x                                    %
%   dgdx:   derivative of gx w.r.t x                                      %
%   dgdphi: derivative of gx w.r.t phi                                    %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zeta = exp(phi(1));
    
    gx = 1 ./ (1 + exp(-zeta * (x - 0.5)));

    if nargout > 1
        s = gx;
        dgdx    = zeta * s .* (1 - s);
        dgdphi  = (x-0.5) .* s .* (1 - s) * zeta; 
    end
end