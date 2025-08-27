% LG OAM multilevel grating + dual-lens multiplexing simulation (MATLAB)
% ------------------------------------------------------------------------
% Based on the original single-level grating + single FFT setup, adds:
%  - Lens 1: perform a second Fourier transform on the first far-field distribution
%  - Place the same forked grating at the back focal plane of lens 1
%  - Lens 2: perform a Fourier transform on the distribution after the second grating to obtain the screen-plane intensity/phase

clear; close all; clc;

%% 1. Basic parameters adjusted for 35 mm film
lambda = 650e-9;       % wavelength [m]
l1      = 1;            % OAM topological charge
l2      = 1;           % OAM topological charge for second grating
w0      = 5e-3;         % Gaussian beam waist radius [m]

% Physical dimensions of 35 mm film plane
film_width  = 36e-3;     % 36 mm width [m]
film_height = 24e-3;     % 24 mm height [m]

% Pixel resolution provided by film service
Nx = 4096;               % number of horizontal pixels
Ny = 2732;               % number of vertical pixels

dx = film_width / Nx;    % horizontal sampling interval [m]
dy = film_height / Ny;   % vertical sampling interval [m]

% Coordinates on film plane
x = (-Nx/2:Nx/2-1) * dx; 
y = (-Ny/2:Ny/2-1) * dy; 
[X, Y] = meshgrid(x, y);
[phi, r] = cart2pol(X, Y);

% grating period and spatial frequency
Lambda  = 10e-6;        
kx      = 2*pi/Lambda;  

%% 2. Compute diffraction angles for ±1…±4 orders
m_orders  = [-4:-1, 1:4];
theta     = asin(m_orders * lambda / Lambda);
theta_deg = theta*180/pi;
for idx = 1:length(m_orders)
    fprintf('m = %+d diffraction angle: %.9f°\n', m_orders(idx), theta_deg(idx));
end

%% 3. Construct forked grayscale grating (Fourier series method)
M            = 2;                     % maximum order
n_list       = -M:M;
user_weights = [2.0, 2.0, 1.0, 2.0, 2.0];
Tn           = user_weights/sum(user_weights);

T1_complex = zeros(size(X));
for idx = 1:length(n_list)
    n = n_list(idx);
    T1_complex = T1_complex + Tn(idx)*exp(1i * n * (kx*X + l1*phi));
end
T1_cont = abs(T1_complex);
T1_cont = (T1_cont - min(T1_cont(:))) / (max(T1_cont(:))-min(T1_cont(:)));

T2_complex = zeros(size(X));
for idx = 1:length(n_list)
    n = n_list(idx);
    T2_complex = T2_complex + Tn(idx)*exp(1i * n * (kx*X + l2*phi));
end
T2_cont = abs(T2_complex);
T2_cont = (T2_cont - min(T2_cont(:))) / (max(T2_cont(:))-min(T2_cont(:)));

% first quantize into 16 discrete levels (0–15)
gray_levels = 256;    % discrete gray levels for grating quantization
mask16_1 = uint8(round((gray_levels-1)*T1_cont));
mask16_2 = uint8(round((gray_levels-1)*T2_cont));

% then expand to full 8-bit range (0–255) for visualization
mask1 = uint8(round(double(mask16_1) * (255/(gray_levels-1))));    % 16-to-256 mapping
mask2 = uint8(round(double(mask16_2) * (255/(gray_levels-1))));    

% Resize masks to match film pixel dimensions (width x height = 4096 x 2732)
film_width_px  = 4096;
film_height_px = 2732;
mask1_film = imresize(mask1, [film_height_px, film_width_px], 'nearest');
mask2_film = imresize(mask2, [film_height_px, film_width_px], 'nearest');

% Save as 8-bit PNG at film resolution
imwrite(mask1_film, 'forked_grating1_fourier_levels_4096x2732.png');
imwrite(mask2_film, 'forked_grating2_fourier_levels_4096x2732.png');

%% 4. Superposition of Gaussian beams at multiple incidence angles
angles_deg     = [1.489858157, 0.744866115, -0.744866115, -1.489858157]; % incident angles for diffraction orders -2, -1, 1, 2
beam_amplitudes = [0.0, 1.0, 1.0, 0.0];     % choose which diffraction orders are enabled
theta_tilt     = angles_deg*pi/180;

E0 = zeros(size(X));
for idx = 1:length(theta_tilt)
    kx_tilt = (2*pi/lambda)*sin(theta_tilt(idx));
    E_temp  = exp(- (X.^2 + Y.^2)/w0^2) .* exp(1i * kx_tilt * X);
    E0 = E0 + beam_amplitudes(idx)*E_temp;
end

%% 5. First grating modulation
mask_norm1 = double(mask1)/(gray_levels-1);
mask_norm2 = double(mask2)/(gray_levels-1);
E1 = E0 .* mask_norm1;

%% 6. First FFT (far-field)
fx = (-Nx/2:Nx/2-1) / film_width;  % spatial frequencies in x
fy = (-Ny/2:Ny/2-1) / film_height; % spatial frequencies in y
[FX, FY] = meshgrid(fx, fy);
E_far = fftshift(fft2(ifftshift(E1))) * dx * dy;(fft2(ifftshift(E1))) * dx * dy;
x_far = lambda * 1 * FX;
y_far = lambda * 1 * FY;

%% 7. Central aperture filtering
fx_cut         = 0.2/Lambda;
H_center       = sqrt(FX.^2 + FY.^2) < fx_cut;
E_far_filtered = E_far .* H_center;

%% 8. OAM purity calculation (truncation radius 0.8*Rmax)
%R_max     = max(r(:));
%r_limit   = 0.8 * R_max;
%mask_r    = r <= r_limit;
%C_l       = sum(sum(E_far_filtered .* exp(-1i*l1*phi) .* mask_r .* r)) * dx^2;
%P_total   = sum(sum(abs(E_far_filtered).^2 .* mask_r .* r)) * dx^2;
%purity    = abs(C_l)^2 / P_total;
%fprintf('OAM mode l=%d purity (radius %.2f mm): %.4f\n', l1, r_limit*1e3, purity);

%% 9. Plot results of first propagation
figure; imshow(mask1, [0 gray_levels-1]); title('First grating mask'); axis off;
figure; imshow(mask2, [0 gray_levels-1]); title('Second grating mask'); axis off;
figure; imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, angle(E_far_filtered));
 axis image; axis xy; colorbar; title('OAM beam Far-field phase [rad]'); xlabel('x (mm)'); ylabel('y (mm)');
figure; imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, abs(E_far).^2);colormap('hot');
 axis image; axis xy; colorbar; title('Far-field real amplitude'); xlabel('x (mm)'); ylabel('y (mm)');
figure; imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, abs(E_far_filtered).^2);
 axis image; axis xy; colorbar; title('OAM beam intensity'); xlabel('x (mm)'); ylabel('y (mm)');
 colormap('hot');

%% 10. Lens 1: perform another Fourier transform on the far-field → near-field
mag = 5;                      % 扩束比
X_exp = X / mag;               % 缩放坐标
Y_exp = Y / mag;
% 对远场场分布做插值重采样，模拟 10× 扩束
E_expanded = interp2(X, Y, E_far_filtered, X_exp, Y_exp, 'linear', 0);

E_back1 = fftshift( fft2( ifftshift(E_expanded) ) ) * dx^2;
%E_back1 = fftshift( fft2( ifftshift(E_far_filtered) ) ) * dx^2;
%E_back1 = E_far_filtered;

% Display E_back1 intensity and phase
% figure('Name','Field after Lens 1 phase');
% imagesc(x_back1(1,:)*1e3, y_back1(:,1)*1e3, angle(E_back1));
% axis image; axis xy; colorbar; title('Lens 1 output phase [rad]'); xlabel('x (mm)'); ylabel('y (mm)');
% figure('Name','Field after Lens 1 intensity');
% imagesc(x_back1(1,:)*1e3, y_back1(:,1)*1e3, abs(E_back1).^2);
% axis image; axis xy; colorbar; title('Lens 1 output intensity'); xlabel('x (mm)'); ylabel('y (mm)'); colormap('hot');

%% 11. Second grating modulation (using the pseudo-inverse matrix of the first grating)
% Define the incidence tilt angle (degrees) for the second grating
% angle2_deg = 5;              
% theta2     = angle2_deg * pi/180;      
% Calculate the corresponding transverse spatial frequency
% kx2 = (2*pi/lambda) * sin(theta2);    

% Apply tilted incidence phase before the grating
% E_back1_tilt = E_back1 .* exp(1i * kx2 * X);

% Then apply the grating modulation
E_after_grating2 = E_back1 .* mask_norm2;

%% 12. Lens 2: another Fourier transform → screen plane
E_screen = fftshift(fft2(ifftshift(E_after_grating2))) * dx^2;
%E_screen = E_after_grating2;

% screen plane coordinates (focal length f2)
f2    = 0.2;  % focal length of Lens 2 [m]
x_scr = lambda * f2 * FX;
y_scr = lambda * f2 * FY;

% original screen-plane field
phase_screen     = angle(E_screen);
intensity_screen = abs(E_screen).^2;

%% 13. Insert central dark-field filter between lens 2 and screen
% Define dark-field radius (m), e.g., 0.5 mm
r_block = 1e-3;
% Construct polar coordinate mask on the screen plane
[~, R_scr] = cart2pol(x_scr, y_scr);
% Region with radius greater than r_block is transmitted
H_block = R_scr > r_block;

% Apply filtering
E_screen_filt      = E_screen .* H_block;
phase_screen_filt  = angle(E_screen_filt);
intensity_screen_filt = abs(E_screen_filt).^2;

%% 14. Plot screen plane distributions (before and after filtering)
% Original phase
figure('Name','Screen Plane Phase');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, phase_screen);
 axis image; axis xy; colorbar;
 title('Screen plane phase [rad] (original)');
 xlabel('x (mm)'); ylabel('y (mm)');

% Original intensity
figure('Name','Screen Plane Intensity');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, intensity_screen);
 axis image; axis xy; colorbar;
 title('Screen plane intensity');
 xlabel('x (mm)'); ylabel('y (mm)');
 colormap('hot');

% Filtered phase
figure('Name','Screen Plane Phase (Filtered)');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, phase_screen_filt);
axis image; axis xy; colorbar;
title(sprintf('Screen plane phase [rad] (r > %.1f mm)', r_block*1e3));
xlabel('x (mm)'); ylabel('y (mm)');

% Filtered intensity
figure('Name','Screen Plane Intensity (Filtered)');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, intensity_screen_filt);
axis image; axis xy; colorbar;
title(sprintf('Screen plane intensity (r > %.1f mm)', r_block*1e3));
xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot');
