%% SIMPLESHADOW
% Compute shadow depth map from a Digital Elevation Model (DEM).
% The algorithm is naive and should no longer be used.
%
% See also: SHADOW
%
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-28

function shadowDepthMap = simpleshadow(Z, varargin)
    %% Initialisation
    % Input parameters
    p = inputParser;
    addRequired(p, 'Z'); % DEM matrix
    addOptional(p, 'Dx', 1);
    addOptional(p, 'azimuth', 45);
    addOptional(p, 'altitude', 45);
    addOptional(p, 'LightVector', [0, 0, 0]);
    addOptional(p, 'ZFactor', 1);
    addOptional(p, 'Quality', 'dirty');
    parse(p, Z, varargin{:});
    Z = p.Results.Z;
    Dx = p.Results.Dx;
    az = p.Results.azimuth;
    al = p.Results.altitude;
    lvector = p.Results.LightVector;
    zFactor = p.Results.ZFactor;
    algorithm = p.Results.Quality;

    % Parameters cleaning
    Z_ex = Z * zFactor / Dx;

    if norm(lvector(:)) ~= 0
        al = asind(lvector(3) / norm(lvector));
        az = atan2d(lvector(2), lvector(1));
    end

    % Prepare and rotate DEM
    [Z_exp, I_exp, J_exp] = rotatefwd(Z_ex, az, algorithm);

    Z_slp = Z_exp - tand(al) * J_exp;
    Z_max = cummax(Z_slp, 1, 'reverse', 'omitnan');
    Z_dpth = Z_slp - Z_max;

    shadowDepthMap = ...
        rotatebwd(I_exp, J_exp, Z_dpth, az, size(Z), algorithm);
end

%% Subfunctions
% Rotate the DEM
function [Zrotated, Inew, Jnew] = rotatefwd(Z, az, algorithm)
    Zpadded = padarray(Z, [1, 1], 'replicate');
    [Ipadded, Jpadded] = meshgrid(1:size(Zpadded, 2), 1:size(Zpadded, 1));
    N_pad = size(Zpadded, 2);
    M_pad = size(Zpadded, 1);
    RotationMatrix = [cosd(az), -sind(az); sind(az), cosd(az)];
    corner1 = RotationMatrix * [1; 1];
    corner2 = RotationMatrix * [1; M_pad];
    corner3 = RotationMatrix * [N_pad; 1];
    corner4 = RotationMatrix * [N_pad; M_pad];

    cornerX = [corner1(1), corner2(1), corner3(1), corner4(1)];
    cornerY = [corner1(2), corner2(2), corner3(2), corner4(2)];
    [Inew, Jnew] = meshgrid(floor(min(cornerX)):ceil(max(cornerX)), ...
        floor(min(cornerY)):ceil(max(cornerY)));

    IJrotated = RotationMatrix * [Ipadded(:)'; Jpadded(:)'];
    Irotated = IJrotated(1, :).';
    Jrotated = IJrotated(2, :).';

    clear IJrotated Ipadded Jpadded corner1 corner2 corner3 corner4

    % Interpolate DEM
    switch algorithm
        case 'high'
            Zrotated = ...
                griddata(Irotated, Jrotated, Zpadded(:), Inew, Jnew);
        case 'dirty'
            Irotated = reshape(Irotated, size(Zpadded));
            Jrotated = reshape(Jrotated, size(Zpadded));
            Zrotated = zeros(size(Inew));
            W = zeros(size(Zrotated));
            
            % offset between rotated and expended grids
            ifst = 1 - floor(min(cornerX));
            jfst = 1 - floor(min(cornerY));

            Ifl = floor(Irotated);
            Jfl = floor(Jrotated);
            T = Irotated - Ifl;
            S = Jrotated - Jfl;
            Ifl = Ifl + ifst;
            Jfl = Jfl + jfst;
            Zrotated(sub2ind(size(Zrotated), Jfl, Ifl)) ...
                = Zrotated(sub2ind(size(Zrotated), Jfl, Ifl)) ...
                + (1 - T) .* (1 - S) .* Zpadded;
            Zrotated(sub2ind(size(Zrotated), Jfl + 1, Ifl)) ...
                = Zrotated(sub2ind(size(Zrotated), Jfl + 1, Ifl)) ...
                + (1 - T) .* S .* Zpadded;
            Zrotated(sub2ind(size(Zrotated), Jfl, Ifl + 1)) ...
                = Zrotated(sub2ind(size(Zrotated), Jfl, Ifl + 1)) ...
                + T .* (1 - S) .* Zpadded;
            Zrotated(sub2ind(size(Zrotated), Jfl + 1, Ifl + 1)) ...
                = Zrotated(sub2ind(size(Zrotated), Jfl + 1, Ifl + 1)) ...
                + T .* S .* Zpadded;
            W(sub2ind(size(Zrotated), Jfl, Ifl)) ...
                = W(sub2ind(size(Zrotated), Jfl, Ifl)) ...
                + (1 - T) .* (1 - S);
            W(sub2ind(size(Zrotated), Jfl + 1, Ifl)) ...
                = W(sub2ind(size(Zrotated), Jfl + 1, Ifl)) ...
                + (1 - T) .* S;
            W(sub2ind(size(Zrotated), Jfl, Ifl + 1)) ...
                = W(sub2ind(size(Zrotated), Jfl, Ifl + 1)) ...
                + T .* (1 - S);
            W(sub2ind(size(Zrotated), Jfl + 1, Ifl + 1)) ...
                = W(sub2ind(size(Zrotated), Jfl + 1, Ifl + 1)) ...
                + T .* S;

            Zrotated(W ~= 0) = Zrotated(W ~= 0) ./ W(W ~= 0);
            Zrotated(W == 0) = nan;
    end

end

% Rotate the DEM back
function Z = rotatebwd(I, J, Z, az, mapSize, algorithm)
    N = mapSize(2);
    M = mapSize(1);
    RotationMatrix = [cosd(-az), -sind(-az); sind(-az), cosd(-az)];

    I = I(:);
    J = J(:);
    Z = Z(:);
    I = I(~isnan(Z));
    J = J(~isnan(Z));
    Z = Z(~isnan(Z));

    IJ_rotb = RotationMatrix * [I'; J'];
    I = IJ_rotb(1, :).';
    J = IJ_rotb(2, :).';

    clear IJ_rotb

    % Interpolate DEM
    switch algorithm
        case 'high'
            [I, J] = meshgrid(1:N, 1:M);
            Z = griddata(I - 1, J - 1, Z, I, J);
        case 'dirty'
            Z = zeros(mapSize);
            W = zeros(mapSize);

            Ifl = floor(I) - 1;
            Jfl = floor(J) - 1;
            T = (I - 1) - Ifl;
            S = (J - 1) - Jfl;

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
                + (1 - T(BL)) .* (1 - S(BL)) .* Z(BL);
            Z(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                = Z(sub2ind(size(Z), Jfl(TL) + 1, Ifl(TL))) ...
                + (1 - T(TL)) .* S(TL) .* Z(TL);
            Z(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                = Z(sub2ind(size(Z), Jfl(BR), Ifl(BR) + 1)) ...
                + T(BR) .* (1 - S(BR)) .* Z(BR);
            Z(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                = Z(sub2ind(size(Z), Jfl(TR) + 1, Ifl(TR) + 1)) ...
                + T(TR) .* S(TR) .* Z(TR);

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
