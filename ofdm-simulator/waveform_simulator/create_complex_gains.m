function [s, P_all, faseis, fi_rand] = create_complex_gains(lambda, d, d_from_origin, d_from_destination, subpaths_per_cluster)

num_scatterers = size(d_from_origin,1);
P_all = calculate_fsl_paths(lambda, d, d_from_origin, d_from_destination);
faseis = calculate_phasechange_due_to_path(lambda, d, d_from_origin , d_from_destination);

% From path power to complex gain:
% First take the square root for gain:
a = sqrt(P_all); 

% Then add the phase shift due to distance:
s = a .* exp(1j*faseis); 

% Create a random phase shift per scatterer:
fi_rand = 0 + (2 * pi) * rand([num_scatterers, subpaths_per_cluster, size(faseis,2),  size(faseis,3)]);
fi_rand = squeeze(sum(fi_rand, 2));

% Then add the random phase shift:
s(2:end,:,:) = s(2:end,:,:) .* exp(1j*fi_rand); 