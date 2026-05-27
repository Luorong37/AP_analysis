function Y_hp = create_temp_highpass(movie)
%CREATE_TEMP_HIGHPASS Build a temporary high-pass filtered movie.
% This is used during motion estimation to improve shift estimation quality.
gSig = 7;
gSiz = 17;
psf = fspecial('gaussian', round(2 * gSiz), gSig);
ind = (psf >= max(psf(:, 1)));
psf = psf - mean(psf(ind));
psf(~ind) = 0;
Y_hp = imfilter(movie, psf, 'symmetric');
end
