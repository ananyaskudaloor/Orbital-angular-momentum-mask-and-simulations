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
l2      = -1;           % OAM topological charge for second grating
w0      = 5e-4;         % Gaussian beam waist radius [m]

% Unified focal length (50 mm) used for all lens-based Fourier transforms
f_lens = 100e-3;

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
Lambda  = 80e-6;        
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
n_list       = [12, 8, 0, 8, 12];
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

% first quantize into 256 discrete levels (0–255)
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
angles_deg     = [0.465533331, -0.465533331];
beam_amplitudes = [1.0, 1.0];     
theta_tilt     = angles_deg*pi/180;

E0 = zeros(size(X));
for idx = 1:length(theta_tilt)
    kx_tilt = (2*pi/lambda)*sin(theta_tilt(idx));
    E_temp  = exp(- (X.^2 + Y.^2)/w0^2) .* exp(1i * kx_tilt * X);
    E0 = E0 + beam_amplitudes(idx)*E_temp;
end

%% 5. First grating modulation
mask_norm1 = double(mask1)/(gray_levels-1);
mask_norm2 = double(mask2)/(gray_levels-1);  % added definition for second grating

E1 = E0 .* mask_norm1;

%% 6. First FFT (far-field) with aperture effect
% Define aperture for the FFT lens
aperture_dia0 = 1 * 25.4e-3;    % 1" diameter aperture [m]
radius0       = aperture_dia0/2; % aperture radius [m]
H_aperture0   = (X.^2 + Y.^2) <= radius0^2;

% Apply aperture effect on field after first grating
E1_ap         = E1 .* H_aperture0;

% Compute spatial frequencies for FFT
fx = (-Nx/2:Nx/2-1) / film_width;   % spatial frequencies in x
fy = (-Ny/2:Ny/2-1) / film_height;  % spatial frequencies in y
[FX, FY] = meshgrid(fx, fy);

% Perform FFT to simulate far-field through finite aperture lens
E_far = fftshift(fft2(ifftshift(E1_ap))) * dx * dy;

% Coordinates in first Fourier plane (using focal length f_lens)
x_far = lambda * f_lens * FX;
y_far = lambda * f_lens * FY;

% sampling in first Fourier plane
dx_far = abs(x_far(1,2)-x_far(1,1));
dy_far = abs(y_far(2,1)-y_far(1,1));

% diagnostics
fprintf('Δx0 = %.3f µm, Δy0 = %.3f µm', dx*1e6, dy*1e6);
fprintf('Δx_F1 = %.3f µm, Δy_F1 = %.3f µm (expected: λf/L)', dx_far*1e6, dy_far*1e6);

%% 7. Central aperture filtering – pinhole (with annotated export)
fx_cut         = 0.6/Lambda;
H_center       = sqrt(FX.^2 + FY.^2) < fx_cut;
E_far_filtered = E_far .* H_center;

% Export pinhole mask at 35 mm 4K resolution
mask_pinhole       = uint8(H_center);  
mask_pinhole_film  = imresize(mask_pinhole, [film_height_px, film_width_px], 'nearest');

% Annotate diameter in top left corner (white text)
diam_pinhole = 2 * lambda * f_lens * fx_cut * 1e3;  % mm
figure('Visible','off'); 
imshow(mask_pinhole_film * 255); 
hold on;
text(10, 10, sprintf('D=%.2f mm', diam_pinhole), ...
     'Color','white','FontSize',18,'FontWeight','bold','Units','pixels');
frame = getframe(gca);
imwrite(frame.cdata, 'pinhole_4096x2732_annotated.jpg');
close;

%% 8. Plot results of first propagation
figure('Name','First grating mask'); imshow(mask1, [0 gray_levels-1]); title('First grating mask'); axis off; set(gcf,'color','w');
figure('Name','Second grating mask'); imshow(mask2, [0 gray_levels-1]); title('Second grating mask'); axis off; set(gcf,'color','w');
figure('Name','OAM beam Far-field phase [rad]');
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, angle(E_far_filtered)); set(gcf,'color','w');
axis image; axis xy; colorbar; title('OAM beam Far-field phase [rad]'); xlabel('x (mm)'); ylabel('y (mm)');
figure('Name','OAM beam intensity');
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, abs(E_far_filtered).^2); set(gcf,'color','w');
axis image; axis xy; colorbar; title('OAM beam intensity'); xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot'); set(gcf,'color','w');
figure('Name','Gaussian beam intensity');
imagesc(x*1e3, y*1e3, abs(E0).^2);
axis image; axis xy; colorbar; title('Gaussian beam intensity'); xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot'); set(gcf,'color','w');
figure('Name','First Grating Modulated Intensity');
imagesc(x*1e3, y*1e3, abs(E1).^2);
axis image; axis xy; colorbar; title('First Grating Modulated'); xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot'); set(gcf,'color','w');
figure('Name','First Grating Far Field Intensity');
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, abs(E_far).^2);
axis image; axis xy; colorbar; title('First Grating Far Field'); xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot'); set(gcf,'color','w');
figure('Name','First Grating Far Field Phase');
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, angle(E_far));
axis image; axis xy; colorbar; title('First Grating Far Field'); xlabel('x (mm)'); ylabel('y (mm)');
colormap('gray'); set(gcf,'color','w');

%% 9. Lens 1: perform Fourier transform with aperture effect (no expansion)
f1 = f_lens;                 % focal length of Lens 1 [m]
aperture_dia1 = 1 * 25.4e-3;  % 1" diameter aperture [m]
radius1 = aperture_dia1/2;   % aperture radius [m]

% Create circular aperture mask in the near-field plane
H_aperture1 = (x_far.^2 + y_far.^2) <= radius1^2;

% Apply aperture effect (no beam expansion)
E_apertured = E_far_filtered .* H_aperture1;

% Perform Fourier transform through Lens 1
E_back1 = fftshift(fft2(ifftshift(E_apertured))) * dx_far * dy_far;

%% 10. Second grating modulation (using the pseudo-inverse matrix of the first grating). Second grating modulation (using the pseudo-inverse matrix of the first grating)
E_after_grating2 = E_back1 .* mask_norm2;

% Lens 2 aperture (optional)
aperture_dia2 = 1 * 25.4e-3;  % 1" diameter aperture [m]
radius2 = aperture_dia2/2;
H_aperture2 = (X.^2 + Y.^2) <= radius2^2; % plane 2 shares sampling with input plane
E_after_grating2 = E_after_grating2 .* H_aperture2;

%% 11. Lens 2: another Fourier transform → screen plane
E_screen = fftshift(fft2(ifftshift(E_after_grating2))) * dx * dy;

% screen plane coordinates (focal length f2)
f2    = f_lens;  % focal length of Lens 2 [m]
x_scr = lambda * f2 * FX;
y_scr = lambda * f2 * FY;

% sampling in screen plane (third Fourier plane)
dx_scr = abs(x_scr(1,2)-x_scr(1,1));
dy_scr = abs(y_scr(2,1)-y_scr(1,1));
fprintf('Δx_F2 = %.3f µm, Δy_F2 = %.3f µm (should equal Δx_F1)', dx_scr*1e6, dy_scr*1e6);

% original screen-plane field
phase_screen     = angle(E_screen);
intensity_screen = abs(E_screen).^2;

%% 12. Dark-field filter – darkfield (with annotated export)
r_block = 1e-4;
[~, R_scr]        = cart2pol(x_scr, y_scr);
H_block           = R_scr > r_block;
E_screen_filt     = E_screen .* H_block;
phase_screen_filt = angle(E_screen_filt);
intensity_screen_filt = abs(E_screen_filt).^2;

% Export dark field filter mask at 35 mm 4K resolution
mask_filter      = uint8(H_block);
mask_filter_film = imresize(mask_filter, [film_height_px, film_width_px], 'nearest');

% Annotate diameter in top left corner (black text)
%diam_dark = 2 * r_block * 1e3;  % mm
%figure('Visible','off');
%imshow(mask_filter_film * 255);
%hold on;
%text(10, 10, sprintf('\u76f4\u5f84: %.2f mm', diam_dark), ...
%     'Color','black','FontSize',18,'FontWeight','bold','Units','pixels');
%frame = getframe(gca);
%imwrite(frame.cdata, 'darkfield_filter_4096x2732_annotated.jpg');
%close;

%% 13. Plot screen plane distributions (before and after filtering)
% Original intensity
figure('Name','Second Grating Modulated Intensity');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, abs(E_after_grating2).^2);
 axis image; axis xy; colorbar;
 title('Second Grating Modulated Intensity');
 xlabel('x (mm)'); ylabel('y (mm)');
 colormap('hot'); set(gcf,'color','w');

figure('Name','Screen Plane Intensity');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, intensity_screen);
 axis image; axis xy; colorbar;
 title('Screen plane intensity');
 xlabel('x (mm)'); ylabel('y (mm)');
 colormap('hot'); set(gcf,'color','w');

 figure('Name','Screen Plane Phase');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, phase_screen);
 axis image; axis xy; colorbar;
 title('Screen plane phase');
 xlabel('x (mm)'); ylabel('y (mm)');
 colormap('gray'); set(gcf,'color','w');

% Filtered intensity
figure('Name','Screen Plane Intensity (Filtered)');
imagesc(x_scr(1,:)*1e3, y_scr(:,1)*1e3, intensity_screen_filt);
axis image; axis xy; colorbar;
title(sprintf('Screen plane intensity (r > %.1f mm)', r_block*1e3));
xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot'); set(gcf,'color','w');
