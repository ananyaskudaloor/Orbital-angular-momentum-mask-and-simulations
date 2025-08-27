clc; clear; close all;

% Parameters
lambda       = 650e-9;       % wavelength (m)
w0           = 2e-3;         % Gaussian waist (2 mm)
sensor_width = 36e-3;
sensor_height= 24e-3;
Nx           = 4096;
Ny           = 2732;

% Grid
x = linspace(-sensor_width/2, sensor_width/2, Nx);
y = linspace(-sensor_height/2, sensor_height/2, Ny);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);

% Input Gaussian beam
E_in = exp(-r.^2 / w0^2);

% Parameters for fork mask
b1 = 1;
b2 = 2;
phi     = atan2(Y, X);
phi_rot = atan2(-X, Y);

% Choose grating period (simulate only one)
Lambda = 200e-6;  % simulate 50 um mask

% Construct forward amplitude mask (grayscale)
mask1 = 0.8 * (1 + cos((2*pi*Y/Lambda) + b1*phi));
mask2 = 0.8 * (1 + cos((2*pi*X/Lambda) + b2*phi_rot));
amp_mask = max(mask1,mask2);
%amp_mask = mask1+mask2;
%amp_mask = mask1*mask2;
%amp_mask = mask2/0.8;
% Modulate the beam
E_mod = E_in .* amp_mask;

% Propagate with FFT (simulate far-field)
E_far = fftshift(fft2(E_mod));
I_far = abs(E_far).^2;

% Normalize and log scale for visualization
I_far_norm = I_far / max(I_far(:));
I_far_log = log10(I_far_norm + 1e-6);  % avoid log(0)

% Display
figure;
imagesc(I_far_log);
axis image;
title('Far-Field Log Intensity after Forward Fork Mask');
xlabel('k_x'); ylabel('k_y');
colormap(hot); colorbar;

% Show the grating (amplitude mask)
figure;
imagesc(x*1e3, y*1e3, amp_mask);  % axes in mm
axis image; set(gca,'YDir','normal');
title('Amplitude Grating (mask)');
xlabel('x [mm]'); ylabel('y [mm]');
colormap(gray); colorbar;