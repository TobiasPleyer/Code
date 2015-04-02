function [p2_out]=compressor_tiltCurve(w,p1_in,p2_in,varargin)
    % Function tiltCurve
    % This function aims to add a linear term to an existing phase curve in
    % order to be better comparable to another one.
    % This is no optimization method but simply uses a mere geometrical
    % idea. The returned phase p2_out will have the same GDD and higher
    % dispersion coefficients as p2_in, but with linear and constant term
    % closer to those of p1_in.
    % It is expected that w is a wavelength range common to both input
    % phase curves!
    if nargin > 3
        interval = varargin{1};
    else
        interval = 25;
    end
    D = p2_in-p1_in;
    idx = auxiliary_find_closest_idx(D,0);
    % Make sure we are not too far off the borders of the array -> index
    % errors
    if idx <= interval
        idx1 = 1;
    else
        idx1 = idx - interval;
    end
    if length(w) <= interval+idx
        idx2 = length(w);
    else
        idx2 = idx + interval;
    end
    x1 = w(idx1);
    x2 = w(idx2);
    y1 = p2_in(idx1) - p1_in(idx1);
    y2 = p2_in(idx2) - p1_in(idx2);
    %p2_out = p2_in + ((y2 - y1)/(x2-x1)).*w + (y1*x2 - y2*x1)/(x2-x1); 
    p2_out = p2_in - ((y2 - y1)/(x2-x1)).*(w-w(idx));
end