function [p2_out]=compressor_alignCurve(w,p1_in,p2_in,varargin)
    % Function alignCurve
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
    idx = auxiliary_find_closest_idx(w,1030);
    p1_val = p1_in(idx);
    p2_val = p2_in(idx);
    gap = p1_val - p2_val;
    p2_in = p2_in + gap;
    p2_out = compressor_tiltCurve(w,p1_in,p2_in,interval);
end