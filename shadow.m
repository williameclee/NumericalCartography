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
