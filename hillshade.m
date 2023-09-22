function [HS, lvector, flat_reflection] = hillshade(Z, varargin)
    % HILLSHADE calculate hillshade of a DEM
    % Syntax
    %   HS = hillshade(Z)
    %   HS = hillshade(Z, "Dx", Dx, ...
    %       "azimuth", az_deg, "altitude", al_deg, ...
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
    %   HS: Hillshade matrix, with the same size as Z
    %   lvector: Light direction vector used in the calculation
    %   flat_reflection: Hillshade value of a flat surface

    % Assigning default parameters
    Dx = 1;
    az = 45;
    al = 45;
    lvector = [0, 0, 0];
    ex_fac = 1;

    % Input parameters
    p = inputParser;
    addRequired(p, "Z"); % DEM matrix
    addOptional(p, "Dx", Dx);
    addOptional(p, "azimuth", az);
    addOptional(p, "altitude", al);
    addOptional(p, "light_vector", lvector);
    addOptional(p, "z_factor", ex_fac);
    parse(p, Z, varargin{:});
    Z = p.Results.Z;
    Dx = p.Results.Dx;
    az = p.Results.azimuth;
    al = p.Results.altitude;
    lvector = p.Results.light_vector;
    ex_fac = p.Results.z_factor;

    % Parameters cleaning
    Z_ex = Z * ex_fac / Dx;

    if norm(lvector(:)) == 0
        lvector = lvector_calc(az, al);
    else
        lvector = lvector / norm(lvector(:));
    end

    lvector = lvector * sign(lvector(3));
    lvector = reshape(lvector, [1, 1, 3]);
    flat_reflection = lvector(3);

    % Calculate hillshade
    N = zeros([size(Z_ex), 3]);
    [N(:, :, 1), N(:, :, 2), N(:, :, 3)] = surfnorm(Z_ex);
    Light_dir = repmat(lvector, [size(Z_ex)]);
    HS = sum(Light_dir .* N, 3);
    HS = max(HS, 0);
end

function lvector = lvector_calc(az, al)
    lvector = [sind(az) * cosd(al), cosd(az) * cosd(al), sind(al)];
    lvector = lvector / norm(lvector(:));
end