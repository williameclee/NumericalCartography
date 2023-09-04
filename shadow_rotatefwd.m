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
