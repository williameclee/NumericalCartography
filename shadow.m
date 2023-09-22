function Z_dp = shadow(Z, varargin)
    % SHADOW calculate depth of shadow casted on a DEM
    % Syntax
    %   HS = hillshade(Z)
    %   HS = hillshade(Z, "Dx", Dx, ...
    %       "azimuth", az, "altitude", al, ...
    %       "z_factor", ex_fac)
    %   HS = hillshade(Z, "Dx", Dx, ...
    %       "light_vector", lvector, ...
    %       "z_factor", ex_fac)
    % Input
    %   Z: DEM matrix in metres. The matrix index should be in the (j,i) format, as of mdgrid
    % Parameters
    %   "Dx": Grid spacing in metres (default = 1)
    %   "azimuth": Azimuth angle in degrees (default = 45)
    %   "altitude": Altitude angle in degrees (default = 45)
    %   "light_vector": Light direction vector (default = [0,0,0])
    %   "z_factor": Vertical exaggeration factor (default = 1)
    % Either "azimuth" and "altitude" or "light_vector" should be specified.
    % Output
    %   Z_dp: Depth of shadow casted on the DEM

    % Assigning default parameters
    Dx = 1;
    az = 45;
    al = 45;
    lvector = [0, 0, 0];
    ex_fac = 1;
    algorithm = "dirty";

    % Input parameters
    p = inputParser;
    addRequired(p, "Z"); % DEM matrix
    addOptional(p, "Dx", Dx);
    addOptional(p, "azimuth", az);
    addOptional(p, "altitude", al);
    addOptional(p, "light_vector", lvector);
    addOptional(p, "z_factor", ex_fac);
    addOptional(p, "quality", algorithm);
    parse(p, Z, varargin{:});
    Z = p.Results.Z;
    Dx = p.Results.Dx;
    az = p.Results.azimuth;
    al = p.Results.altitude;
    lvector = p.Results.light_vector;
    ex_fac = p.Results.z_factor;
    algorithm = p.Results.quality;

    % Parameters cleaning
    Z_ex = Z * ex_fac / Dx;

    if norm(lvector(:)) ~= 0
        al = asind(lvector(3) / norm(lvector));
        az = atan2d(lvector(2), lvector(1));
    end

    % Prepare and rotate DEM
    [Z_exp, I_exp, J_exp] = shadow_rotatefwd(Z_ex, az, algorithm);

    Z_slp = Z_exp - tand(al) * J_exp;
    Z_max = cummax(Z_slp, 1, "reverse", "omitnan");
    Z_dpth = Z_slp - Z_max;

    Z_dp = shadow_rotatebwd(I_exp, J_exp, Z_dpth, az, size(Z), algorithm);
end

%% Supplementary functions
% Rotate DEM
function [Z_exp, I_exp, J_exp] = shadow_rotatefwd(Z, az, algorithm)
    Z_pad = padarray(Z, [1, 1], 'replicate');
    [I_pad, J_pad] = meshgrid(1:size(Z_pad, 2), 1:size(Z_pad, 1));
    N_pad = size(Z_pad, 2);
    M_pad = size(Z_pad, 1);
    RM = [cosd(az), -sind(az); sind(az), cosd(az)];
    C1_rot = RM * [1; 1];
    C2_rot = RM * [1; M_pad];
    C3_rot = RM * [N_pad; 1];
    C4_rot = RM * [N_pad; M_pad];

    Cx_rot = [C1_rot(1), C2_rot(1), C3_rot(1), C4_rot(1)];
    Cy_rot = [C1_rot(2), C2_rot(2), C3_rot(2), C4_rot(2)];
    [I_exp, J_exp] = meshgrid(floor(min(Cx_rot)):ceil(max(Cx_rot)), ...
        floor(min(Cy_rot)):ceil(max(Cy_rot)));

    IJ_rot = RM * [I_pad(:)'; J_pad(:)'];
    I_rot = IJ_rot(1, :).';
    J_rot = IJ_rot(2, :).';

    clear IJ_rot I_pad J_pad C1_rot C2_rot C3_rot C4_rot

    % Interpolate DEM
    switch algorithm
        case "high"
            Z_exp = griddata(I_rot, J_rot, Z_pad(:), I_exp, J_exp);
        case "dirty"
            I_rot = reshape(I_rot, size(Z_pad));
            J_rot = reshape(J_rot, size(Z_pad));
            Z_exp = zeros(size(I_exp));
            W = zeros(size(Z_exp));
            
            % offset between rotated and expended grids
            ifst = 1 - floor(min(Cx_rot));
            jfst = 1 - floor(min(Cy_rot));

            Ifl = floor(I_rot);
            Jfl = floor(J_rot);
            T = I_rot - Ifl;
            S = J_rot - Jfl;
            Ifl = Ifl + ifst;
            Jfl = Jfl + jfst;
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

end

% Rotate DEM back
function Z = shadow_rotatebwd(I_exp, J_exp, Z_exp, az, MN, algorithm)
    N = MN(2);
    M = MN(1);
    RM_b = [cosd(-az), -sind(-az); sind(-az), cosd(-az)];

    I_exp = I_exp(:);
    J_exp = J_exp(:);
    Z_exp = Z_exp(:);
    I_exp = I_exp(~isnan(Z_exp));
    J_exp = J_exp(~isnan(Z_exp));
    Z_exp = Z_exp(~isnan(Z_exp));

    IJ_rotb = RM_b * [I_exp'; J_exp'];
    I_exp = IJ_rotb(1, :).';
    J_exp = IJ_rotb(2, :).';

    clear IJ_rotb

    % Interpolate DEM
    switch algorithm
        case "high"
            [I, J] = meshgrid(1:N, 1:M);
            Z = griddata(I_exp - 1, J_exp - 1, Z_exp, I, J);
        case "dirty"
            Z = zeros(MN);
            W = zeros(MN);

            Ifl = floor(I_exp) - 1;
            Jfl = floor(J_exp) - 1;
            T = (I_exp - 1) - Ifl;
            S = (J_exp - 1) - Jfl;

            BL = (Ifl >= 1) & (Jfl >= 1) ...
                & (Ifl <= N) & (Jfl <= M);
            BR = (Ifl + 1 >= 1) & (Jfl >= 1) ...
                & (Ifl + 1 <= N) & (Jfl <= M);
            TL = (Ifl >= 1) & (Jfl + 1 >= 1) ...
                & (Ifl <= N) & (Jfl + 1 <= M);
            TR = (Ifl + 1 >= 1) & (Jfl + 1 >= 1) ...
                & (Ifl + 1 <= N) & (Jfl + 1 <= M);

            Z(sub2ind(size(Z), Jfl(BL), Ifl(BL))) ...
                = Z(sub2ind(size(Z), Jfl(BL), Ifl(BL))) ...
                + (1 - T(BL)) .* (1 - S(BL)) .* Z_exp(BL);
            Z(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                = Z(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                + (1 - T(TL)) .* S(TL) .* Z_exp(TL);
            Z(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                = Z(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                + T(BR) .* (1 - S(BR)) .* Z_exp(BR);
            Z(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                = Z(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                + T(TR) .* S(TR) .* Z_exp(TR);

            W(sub2ind(size(Z), Jfl(BL), Ifl(BL))) ...
                = W(sub2ind(size(Z), Jfl(BL), Ifl(BL))) ...
                + (1 - T(BL)) .* (1 - S(BL));
            W(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                = W(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                + (1 - T(TL)) .* S(TL);
            W(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                = W(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                + T(BR) .* (1 - S(BR));
            W(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                = W(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                + T(TR) .* S(TR);

            Z(W ~= 0) = Z(W ~= 0) ./ W(W ~= 0);
            Z(W == 0) = nan;
    end

end
