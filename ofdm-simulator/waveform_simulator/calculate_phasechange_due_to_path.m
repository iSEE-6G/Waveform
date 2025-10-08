function faseis = calculate_phasechange_due_to_path(lambda, d, d_from_origin , d_from_destination)

num_scatterers = size(d_from_origin, 1);
num_rus = size(d_from_origin, 2);
num_ues = size(d_from_destination, 2);

% Phase due to distance for LoS:
fash = 2*pi/lambda.*d;

faseis = zeros(num_scatterers+1, num_rus, num_ues);
% Phase due to distance for multipath:
for rr = 1 : num_rus
    for uu = 1 : num_ues
        faseis(2:end, rr, uu) = 2*pi/lambda.*(d_from_origin(:, rr) + d_from_destination(:,uu));
    end
end
faseis(1,:,:) = fash;

