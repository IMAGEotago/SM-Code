function [xnext, dfdx, dfdtheta] = f_rw(x, theta, u, ~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rescorla-Wagner Evolution Function                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% Author: Sophie Morris                                                   %
% Created: 14/11/2025                                                     %
% ----------------------------------------------------------------------- %
% Description: This script contains the Rescorla-Wagner evolution         %
% used in Sophie M's Masters' project. A timestamped anaysis plan is      %
% available on (https://github.com/IMAGEotago).                           %
% ----------------------------------------------------------------------- %
% Inputs:                                                                 %
%   x:        hidden state                                                %
%   theta:    evolution parameter                                         %
%   u:        stimulus outcome                                            %
% Outputs:                                                                %
%   xnext:    predicted x                                                 %
%   dfdx:     derivative of xnext w.r.t x                                 %
%   dfdtheta: derivative of xnext w.r.t theta                             %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha = 1/(1+exp(-theta(1)));
    prediction_error = u(1) - x; 
    xnext = x + alpha * prediction_error; 

    if nargout > 1 
        dfdx = 1 - alpha; 
        dfdtheta = prediction_error * alpha * (1 - alpha);
    end
end