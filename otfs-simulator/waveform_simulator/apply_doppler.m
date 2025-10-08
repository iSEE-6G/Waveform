function s_new = apply_doppler(time_vector, direct_freq_doppler, fd_iT, fd_iR, s)

num_ru = size(fd_iT, 2);
num_ue = size(fd_iR, 2);
num_scatterers = size(fd_iT, 1);

s_new = zeros(num_ru, num_ue, num_scatterers+1, length(time_vector));
for ii = 1 : num_ru
    for ll = 1 : num_ue
        % TO DIRECT:
        tmp_new = squeeze(s(1, ii, ll)) * exp(2j*pi*direct_freq_doppler(ii,ll)*time_vector);

        % TA SCATTERERS:
        tmp_new_all = squeeze(s(2:end, ii, ll)) .* exp(2j*pi*fd_iT(:,ii)*time_vector) .* exp(2j*pi*fd_iR(:,ll)*time_vector);
        
        s_new(ii, ll, 1, :) = tmp_new;
        s_new(ii, ll, 2:end, :) = tmp_new_all;
    end
end
