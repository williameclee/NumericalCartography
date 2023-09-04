function Z_dp = shadow(Z, varargin)
    % SHADOW calculate depth of shadow casted on a DEM
    % Syntax
    %   HS = hillshade(Z)
    %   HS = hillshade(Z, "Dx", Dx, ...
    %       "azimuth", az, "altitude", al, ...
    %       "z_factor", ex_fac)
    %   HS = hillshade(Z, "Dx", Dx, ...
    %       "light_direction", light_dir, ...
    %       "z_factor", ex_fac)
    % Input
    %   Z: DEM matrix in metres. The matrix index should be in the (j,i) format, as of mdgrid
    % Parameters
    %   "Dx": Grid spacing in metres (default = 1)
    %   "azimuth": Azimuth angle in degrees (default = 45)
    %   "altitude": Altitude angle in degrees (default = 45)
    %   "light_direction": Light direction vector (default = [0,0,0])
    %   "z_factor": Vertical exaggeration factor (default = 1)
    % Either "azimuth" and "altitude" or "light_direction" should be specified.
    % Output
    %   Z_dp: Depth of shadow casted on the DEM

    % Assigning default parameters
    Dx = 1;
    az = 45;
    al = 45;
    light_dir = [0, 0, 0];
    ex_fac = 1;

    % Input parameters
    p = inputParser;
    addRequired(p, "Z"); % DEM matrix
    addOptional(p, "Dx", Dx);
    addOptional(p, "azimuth", az);
    addOptional(p, "altitude", al);
    addOptional(p, "light_direction", light_dir);
    addOptional(p, "z_factor", ex_fac);
    parse(p, Z, varargin{:});
    Z = p.Results.Z;
    Dx = p.Results.Dx;
    az = p.Results.azimuth;
    al = p.Results.altitude;
    light_dir = p.Results.light_direction;
    ex_fac = p.Results.z_factor;

    % Parameters cleaning
    Z_ex = Z * ex_fac / Dx;

    if norm(light_dir(:)) ~= 0
        al = asind(light_dir(3) / norm(light_dir));
        az = atan2d(light_dir(2), light_dir(1));
    end

    % Prepare and rotate DEM
    [Z_exp, I_exp, J_exp] = shadow_rotatefwd(Z_ex, az);

    Z_slp = Z_exp - tand(al) * J_exp;
    Z_max = cummax(Z_slp, 1, "reverse", "omitnan");
    Z_dpth = Z_slp - Z_max;

    Z_dp = shadow_rotatebwd(I_exp, J_exp, Z_dpth, az, size(Z, 2) + 2, size(Z, 1) + 2);
end


function [Z_exp, I_exp, J_exp] = shadow_rotatefwd(Z, az)
    Z_pad = padarray(Z, [1, 1], 'replicate');
    [I_pad, J_pad] = meshgrid(1:size(Z_pad, 2), 1:size(Z_pad, 1));
    N = size(Z_pad, 2);
    M = size(Z_pad, 1);
    RM = [cosd(az), -sind(az); sind(az), cosd(az)];
    C1 = RM * [1; 1];
    C2 = RM * [1; M];
    C3 = RM * [N; 1];
    C4 = RM * [N; M];

    Cx = [C1(1), C2(1), C3(1), C4(1)];
    Cy = [C1(2), C2(2), C3(2), C4(2)];
    [I_exp, J_exp] = meshgrid(floor(min(Cx)):ceil(max(Cx)), ... 
        floor(min(Cy)):ceil(max(Cy)));

    IJ_rot = RM * [I_pad(:)'; J_pad(:)'];
    I_rot = IJ_rot(1, :).';
    J_rot = IJ_rot(2, :).';

    I_rot = reshape(I_rot, size(Z_pad));
    J_rot = reshape(J_rot, size(Z_pad));

    % Interpolate DEM
    Z_exp = zeros(size(I_exp));
    W = zeros(size(Z_exp));

    iofst = 1 - floor(min(Cx));
    jofst = 1 - floor(min(Cy));

    Ifl = floor(I_rot);
    Jfl = floor(J_rot);
    T = I_rot - Ifl;
    S = J_rot - Jfl;
    Ifl = Ifl + iofst;
    Jfl = Jfl + jofst;
    Z_exp(sub2ind(size(Z_exp), Jfl, Ifl)) ...
        = Z_exp(sub2ind(size(Z_exp), Jfl, Ifl)) ...
        + (1 - T) .* (1 - S) .* Z_pad;
    Z_exp(sub2ind(size(Z_exp), Jfl + 1, Ifl)) ...
        = Z_exp(sub2ind(size(Z_exp), Jfl + 1, Ifl)) ...
        + (1 - T) .* S .* Z_pad;
    Z_exp(sub2ind(size(Z_exp), Jfl, Ifl + 1)) ...
        = Z_exp(sub2ind(size(Z_exp), Jfl, Ifl + 1)) ...
        + T .* (1 - S) .* Z_pad;
    Z_exp(sub2ind(size(Z_exp), Jfl + 1, Ifl + 1)) ...
        = Z_exp(sub2ind(size(Z_exp), Jfl + 1, Ifl + 1)) ...
        + T .* S .* Z_pad;
    W(sub2ind(size(Z_exp), Jfl, Ifl)) ...
        = W(sub2ind(size(Z_exp), Jfl, Ifl)) ...
        + (1 - T) .* (1 - S);
    W(sub2ind(size(Z_exp), Jfl + 1, Ifl)) ...
        = W(sub2ind(size(Z_exp), Jfl + 1, Ifl)) ...
        + (1 - T) .* S;
    W(sub2ind(size(Z_exp), Jfl, Ifl + 1)) ...
        = W(sub2ind(size(Z_exp), Jfl, Ifl + 1)) ...
        + T .* (1 - S);
    W(sub2ind(size(Z_exp), Jfl + 1, Ifl + 1)) ...
        = W(sub2ind(size(Z_exp), Jfl + 1, Ifl + 1)) ...
        + T .* S;

    Z_exp(W ~= 0) = Z_exp(W ~= 0) ./ W(W ~= 0);
    Z_exp(W == 0) = nan;
end

function Z_dp = shadow_rotatebwd(I_exp, J_exp, Z_dpth, az, N, M)
    RM_b = [cosd(-az), -sind(-az); sind(-az), cosd(-az)];

    I_exp = I_exp(:);
    J_exp = J_exp(:);
    Z_dpth = Z_dpth(:);
    I_exp = I_exp(~isnan(Z_dpth));
    J_exp = J_exp(~isnan(Z_dpth));
    Z_dpth = Z_dpth(~isnan(Z_dpth));

    IJ_rotb = RM_b * [I_exp'; J_exp'];
    I_exp = IJ_rotb(1, :).';
    J_exp = IJ_rotb(2, :).';

    Z_dp = zeros([M, N]);
    W = zeros([M, N]);

    Ifl = floor(I_exp);
    Jfl = floor(J_exp);
    T = I_exp - Ifl;
    S = J_exp - Jfl;

    BL = (Ifl >= 1) & (Jfl >= 1) & (Ifl <= N) & (Jfl <= M);
    BR = (Ifl + 1 >= 1) & (Jfl >= 1) & (Ifl + 1 <= N) & (Jfl <= M);
    TL = (Ifl >= 1) & (Jfl + 1 >= 1) & (Ifl <= N) & (Jfl + 1 <= M);
    TR = (Ifl + 1 >= 1) & (Jfl + 1 >= 1) & (Ifl + 1 <= N) & (Jfl + 1 <= M);

    Z_dp(sub2ind(size(Z_dp), Jfl(BL), Ifl(BL))) ...
        = Z_dp(sub2ind(size(Z_dp), Jfl(BL), Ifl(BL))) ...
        + (1 - T(BL)) .* (1 - S(BL)) .* Z_dpth(BL);
    Z_dp(sub2ind(size(Z_dp), Jfl(TL) + 1, Ifl(TL))) ...
        = Z_dp(sub2ind(size(Z_dp), Jfl(TL) + 1, Ifl(TL))) ...
        + (1 - T(TL)) .* S(TL) .* Z_dpth(TL);
    Z_dp(sub2ind(size(Z_dp), Jfl(BR), Ifl(BR) + 1)) ...
        = Z_dp(sub2ind(size(Z_dp), Jfl(BR), Ifl(BR) + 1)) ...
        + T(BR) .* (1 - S(BR)) .* Z_dpth(BR);
    Z_dp(sub2ind(size(Z_dp), Jfl(TR) + 1, Ifl(TR) + 1)) ...
        = Z_dp(sub2ind(size(Z_dp), Jfl(TR) + 1, Ifl(TR) + 1)) ...
        + T(TR) .* S(TR) .* Z_dpth(TR);

    W(sub2ind(size(Z_dp), Jfl(BL), Ifl(BL))) ...
        = W(sub2ind(size(Z_dp), Jfl(BL), Ifl(BL))) ...
        + (1 - T(BL)) .* (1 - S(BL));
    W(sub2ind(size(Z_dp), Jfl(TL) + 1, Ifl(TL))) ...
        = W(sub2ind(size(Z_dp), Jfl(TL) + 1, Ifl(TL))) ...
        + (1 - T(TL)) .* S(TL);
    W(sub2ind(size(Z_dp), Jfl(BR), Ifl(BR) + 1)) ...
        = W(sub2ind(size(Z_dp), Jfl(BR), Ifl(BR) + 1)) ...
        + T(BR) .* (1 - S(BR));
    W(sub2ind(size(Z_dp), Jfl(TR) + 1, Ifl(TR) + 1)) ...
        = W(sub2ind(size(Z_dp), Jfl(TR) + 1, Ifl(TR) + 1)) ...
        + T(TR) .* S(TR);

    Z_dp(W ~= 0) = Z_dp(W ~= 0) ./ W(W ~= 0);
    Z_dp(W == 0) = nan;

    Z_dp = Z_dp(2:end - 1, 2:end - 1);
end
