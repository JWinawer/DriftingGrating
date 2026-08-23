function H = harmonic_decompose(Y, thetaV, expCfg, b)
% HARMONIC_DECOMPOSE  Per-vertex harmonic amplitudes, and what the model predicts them to be.
%
%   H = harmonic_decompose(Y, thetaV, expCfg, b)
%
% Four orientations spaced 45 deg give each vertex's demeaned response vector exactly
% THREE degrees of freedom, spanned by cos(2a), sin(2a), cos(4a) for a sampling angle
% a (the fourth basis function, sin(4a), is identically zero at 45-deg spacing). Call
% the per-vertex amplitudes A, B, C. The four-term model is then *exactly equivalent*
% to three scalar regressions of A, B, C on the vertex's polar angle -- which is what
% makes it interpretable rather than a black box, and is the basis of the main
% diagnostic figure.
%
% The frame is chosen so the four sampling angles are the same at every vertex:
%   Cartesian (isPolar==false): ABSOLUTE frame,        a = oriAngle
%       A = b1 + b3*cos(2*thetaV)
%       B =      b3*sin(2*thetaV)
%       C = b2 + b4*cos(4*thetaV)
%   Polar (isPolar==true):      RADIAL-RELATIVE frame, a = oriAngle - 90
%       A = b3 + b1*cos(2*thetaV)
%       B =    - b1*sin(2*thetaV)
%       C = b4 + b2*cos(4*thetaV)
% i.e. the two experiments are exact mirror images, with the absolute and polar terms
% swapping roles. Read off the first pair: b1 and b3 are separated ONLY by how A is
% modulated by thetaV, plus the B channel. That partition is the whole question of how
% much apparent horizontal-vs-vertical is radial-vs-tangential leaking through geometry.
%
% Inputs
%   Y      : nVertex x nOri demeaned responses (from HARMONIC_VERTEX_DATA).
%   thetaV : nVertex x 1 polar angles, degrees.
%   expCfg : cfg.dg or cfg.da.
%   b      : optional 1 x 4 coefficients; if given, the predicted curves are returned.
%
% Output H
%   .A, .B, .C   nVertex x 1 observed amplitudes
%   .Ahat, .Bhat, .Chat  nVertex x 1 model curves (only if b supplied)
%   .a           1 x nOri sampling angles, degrees
%   .frame       'absolute' or 'radial'
%   .reconErr    max abs error of the 3-term reconstruction (should be ~1e-14;
%                a nonzero value means the sampling is not 45-deg spaced)

    ori = double(expCfg.oriAngle(:)).';
    if expCfg.isPolar
        H.frame = 'radial';
        a = mod(ori - 90, 180);
    else
        H.frame = 'absolute';
        a = mod(ori, 180);
    end
    H.a = a;

    nO = numel(a);
    c2 = cosd(2*a);  s2 = sind(2*a);  c4 = cosd(4*a);

    H.A = (2/nO) * (Y * c2.');
    H.B = (2/nO) * (Y * s2.');
    H.C = (1/nO) * (Y * c4.');

    % The three amplitudes must reconstruct the demeaned data exactly.
    recon = H.A*c2 + H.B*s2 + H.C*c4;
    H.reconErr = max(abs(recon(:) - Y(:)));

    if nargin >= 4 && ~isempty(b)
        bb = b(:); bb(~isfinite(bb)) = 0;
        tv = double(thetaV(:));
        if expCfg.isPolar
            H.Ahat = bb(3) + bb(1)*cosd(2*tv);
            H.Bhat =       - bb(1)*sind(2*tv);
            H.Chat = bb(4) + bb(2)*cosd(4*tv);
        else
            H.Ahat = bb(1) + bb(3)*cosd(2*tv);
            H.Bhat =         bb(3)*sind(2*tv);
            H.Chat = bb(2) + bb(4)*cosd(4*tv);
        end
    end
end
