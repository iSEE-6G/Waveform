function [fd_iT, fd_iR, direct_freq_doppler] = doppler_creation(lambda, num_scatterers, u_MaxTR, u_Max, num_ru, num_ue)

if nargin<3
    % velocity for Tx and Rx (max 50 m/s)
    u_MaxTR = 0;
    % velocity for scatterers
    u_Max = 70;
end

%% Doppler and movement:
u_T = rand(num_ru,1) * u_MaxTR;
u_R = rand(num_ue,1) * u_MaxTR;

% Direction of movement:
f_T = 0 + (2 * pi) * rand(num_ru,1);
f_R = 0 + (2 * pi) * rand(num_ue,1);

% Velocity in each axis
u_Tx = u_T.* cos(f_T);
u_Ty = u_T.* sin(f_T);
u_Rx = u_R.* cos(f_R);
u_Ry = u_R.* sin(f_R);

% Relative direction of Rx-Tx movement
f_relative_TR = zeros(num_ru, num_ue);
relative_velocity = zeros(num_ru, num_ue);
for ii = 1 : num_ru
    for ll = 1 : num_ue
        f_relative_TR(ii, ll) = f_R(ll) - f_T(ii);
        relative_velocity_x = u_Rx(ll) - u_Tx(ii);
        relative_velocity_y = u_Ry(ll) - u_Ty(ii);
        relative_velocity(ii, ll) = sqrt((relative_velocity_x ^ 2) + (relative_velocity_y ^ 2)); % // sxetikh taxythta

    end
end

%% H dieuthinsi tis sxetikhs taxutitas einai: f_relative_TR opote:
direct_freq_doppler = (relative_velocity / lambda).* cos(f_relative_TR);

% Movement of the scatterers:
u_scatterers = -u_Max + rand(num_scatterers,1) * 2*u_Max;
f_scat = 2*pi*rand(num_scatterers,1);
u_scatterers_x = u_scatterers.*cos(f_scat);
u_scatterers_y = u_scatterers.*sin(f_scat);

fd_iT = zeros(num_scatterers, num_ru);
fd_iR = zeros(num_scatterers, num_ue);
for ii = 1 : num_ru
    % TWRA BRISKW TIS SXETIKES GWNIES KINISIS APO POMPO SE SKEDASTI:
    f_relative_Tscat = f_scat - f_T(ii);
    relative_scatterers_Tx = u_scatterers_x - u_Tx(ii); % // sxetikh taxythta apo pompo se skedasth
    relative_scatterers_Ty = u_scatterers_y - u_Ty(ii);

    relative_velocity_scatterers_tx = sqrt((relative_scatterers_Tx.^ 2) + (relative_scatterers_Ty.^ 2));
    % KAI ANTISTOIXA DYO DOPPLER SHIFTS!: PROSOXI BAZOYME THN SXETIKI GWNIA (STO ALLAKSA):
    fd_iT(:,ii) = (relative_velocity_scatterers_tx / lambda).* cos(f_relative_Tscat);
end
for ll = 1 : num_ue
    % TWRA BRISKW TIS SXETIKES GWNIES KINISIS APO SKEDASTI SE POMPO:
    f_relative_scatR = f_R(ll) - f_scat;
    relative_scatterers_Rx = u_Rx(ll) - u_scatterers_x; % // sxetikh taxythta apo skedasth se dekth
    relative_scatterers_Ry = u_Ry(ll) - u_scatterers_y;

    relative_velocity_scatterers_rx = sqrt((relative_scatterers_Rx.^ 2) + (relative_scatterers_Ry.^ 2));

    % KAI ANTISTOIXA DYO DOPPLER SHIFTS!: PROSOXI BAZOYME THN SXETIKI GWNIA (STO ALLAKSA):
    fd_iR(:,ll) = (relative_velocity_scatterers_rx / lambda).* cos(f_relative_scatR);
end
