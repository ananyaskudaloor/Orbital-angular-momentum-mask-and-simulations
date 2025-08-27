clear; close all; clc;

%% 1. Basic parameters
lambda = 650e-9;         % wavelength [m]
k = 2*pi/lambda;         % wavenumber
w0_sim = 0.5e-3;         % LG beam waist [m]
l = 2;                   % OAM topological charge
x_offset = 5e-3;       % beam shift in x-direction [m]
z = 1;                   % propagation distance (initial plane)

% Film / grating size
Nx = 4096; Ny = 2732;
film_width  = 36e-3;
film_height = 24e-3;
dx = film_width / Nx;
dy = film_height / Ny;

x = (-Nx/2:Nx/2-1) * dx;
y = (-Ny/2:Ny/2-1) * dy;
[X, Y] = meshgrid(x, y);

%% === Coordinates for beam and grating ===
Xb = X - x_offset;        % Shifted vortex beam
Yb = Y;
R_beam = sqrt(Xb.^2 + Yb.^2);
Phi_beam = atan2(Yb, Xb);

R_grating = sqrt(X.^2 + Y.^2);   % Grating stays centered

%% === Laguerre-Gaussian beam (p = 0) ===
LG_beam = (sqrt(2) * R_beam / w0_sim).^abs(l) ...
    .* exp(1i * l * Phi_beam) ...
    .* exp(-R_beam.^2 / w0_sim^2) ...
    .* exp(-1i * k * z);

%% === Construct annular grating ===
Lambda = 50e-6;
k_radial = 2*pi / Lambda;

annular_grating = 0.5*(1 + cos(k_radial * R_grating));
annular_grating = (annular_grating - min(annular_grating(:))) / (max(annular_grating(:)) - min(annular_grating(:)));

%% === Display grating and beam before mask ===
figure;
imshow(uint8(annular_grating*255));
title('Annular Grating Mask');
axis on;

figure;
imagesc(x*1e3, y*1e3, abs(LG_beam));
axis image; axis xy; colorbar;
title('Amplitude of shifted LG beam (l=2) before mask');
xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot');

%% === Apply grating ===
E_after_grating = LG_beam .* annular_grating;

%% === Far-field computation via FFT ===
fx = (-Nx/2:Nx/2-1)/film_width;
fy = (-Ny/2:Ny/2-1)/film_height;
[FX, FY] = meshgrid(fx, fy);

E_far = fftshift(fft2(ifftshift(E_after_grating))) * dx * dy;

x_far = lambda * FX;
y_far = lambda * FY;

intensity_far = abs(E_far).^2;
phase_far = angle(E_far);

%% === Plot far-field ===
figure;
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, intensity_far);
axis image; axis xy; colorbar;
title('Far-field intensity of shifted LG l=2 through annular grating');
xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot');

figure;
imagesc(x_far(1,:)*1e3, y_far(:,1)*1e3, phase_far);
axis image; axis xy; colorbar;
title('Far-field phase of shifted LG l=2 through annular grating');
xlabel('x (mm)'); ylabel('y (mm)');

%% === Near-field visualization (just after grating) ===

% Amplitude of the near-field
figure;
imagesc(x*1e3, y*1e3, abs(E_after_grating));
axis image; axis xy; colorbar;
title('Near-field amplitude of shifted LG l=2 after annular grating');
xlabel('x (mm)'); ylabel('y (mm)');
colormap('hot');

% Phase of the near-field
figure;
imagesc(x*1e3, y*1e3, angle(E_after_grating));
axis image; axis xy; colorbar;
title('Near-field phase of shifted LG l=2 after annular grating');
xlabel('x (mm)'); ylabel('y (mm)');


%% === Save grating as JPEG ===
imwrite(uint8(annular_grating*255), 'annular_final.jpg');
disp('Annular grating mask saved as annular_final.jpg');
