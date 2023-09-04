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
