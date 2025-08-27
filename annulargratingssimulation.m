clc; clear; close all;

%% === Parameters ===
lambda = 0.65;                 % wavelength (arbitrary units)
k = 2*pi/lambda;               % wavenumber
z = 1000;                      % propagation distance
w0 = 0.5;                      % beam waist
grating_period = 0.08;         % annular grating period (radial)
x_offset = 1.5;                % beam offset to enhance fringes

%% === Grid ===
grid_size = 1000;
x_max = 5;
x = linspace(-x_max, x_max, grid_size);
y = x;
[X,Y] = meshgrid(x,y);
dx = x(2)-x(1);

% Coordinates
R  = sqrt(X.^2 + Y.^2);        % grating centered
Xb = X - x_offset;             % beam shifted
Yb = Y;
Rb = sqrt(Xb.^2 + Yb.^2);
Phib = atan2(Yb, Xb);

% Gratings
t_amp   = double(mod(R, grating_period) < grating_period/2); % amplitude annuli (binary)
t_phase = exp(1i * 2*pi * R / grating_period);               % phase annular ramp

% Fresnel quadratic phase for FFT propagation
quad_phase = exp(1i * k/(2*z) * (X.^2 + Y.^2));

% Helper
farfield = @(U) (abs(fftshift(fft2(U .* quad_phase)) * dx^2)).^2;

ells = -3:3;

%% ---------- Amplitude grating results ----------
figure('Name','Amplitude grating l sweep');
tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');
for ii = 1:numel(ells)
    l = ells(ii);
    LG = (sqrt(2)*Rb/w0).^abs(l) .* exp(1i*l*Phib) .* exp(-(Rb.^2)/(w0^2)) .* exp(-1i*k*z);
    U  = LG .* t_amp;
    I  = farfield(U); I = I / max(I(:));

    nexttile;
    imagesc(x, y, log(1+I));
    axis image off; colormap hot;
    title(sprintf('l = %+d', l));
end
sgtitle('Far-Field Intensity (log) — Amplitude Annular Grating');

%% ---------- Phase grating results ----------
figure('Name','Phase grating l sweep');
tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');
for ii = 1:numel(ells)
    l = ells(ii);
    LG = (sqrt(2)*Rb/w0).^abs(l) .* exp(1i*l*Phib) .* exp(-(Rb.^2)/(w0^2)) .* exp(-1i*k*z);
    U  = LG .* t_phase;
    I  = farfield(U); I = I / max(I(:));

    nexttile;
    imagesc(x, y, log(1+I));
    axis image off; colormap hot;
    title(sprintf('l = %+d', l));
end
sgtitle('Far-Field Intensity (log) — Phase Annular Grating');
