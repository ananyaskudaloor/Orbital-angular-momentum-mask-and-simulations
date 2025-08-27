close all;
clc;
clear;

% --- Parameters for the simulation window and mask image ---
sensor_width = 36e-3;        % 35 mm
sensor_height = 24e-3;     % 23.3 mm
pixels_x = 4096;             % Horizontal resolution
pixels_y = 2732;             % Vertical resolution

% Spatial resolution
dx = sensor_width / pixels_x;
dy = sensor_height / pixels_y;

% Coordinate grid
x = linspace(-sensor_width/2, sensor_width/2, pixels_x);
y = linspace(-sensor_height/2, sensor_height/2, pixels_y);
[X, Y] = meshgrid(x, y);

% Topological charges
b1 = -2;   % Original fork
b2 = -1;   % Rotated fork

% Phase terms
phi = atan2(Y, X);
phi_rot = atan2(-X, Y);

% --- Loop over different grating periods ---
grating_periods = [20e-6, 50e-6, 100e-6];  % in meters

for i = 1:length(grating_periods)
    Lambda = grating_periods(i);
    
    % Fork gratings (smooth grayscale modulation)
    mask1 = 0.8 * (1 + cos((2*pi*Y/Lambda) + b1*phi));
    mask2 = 0.8 * (1 + cos((2*pi*X/Lambda) + b2*phi_rot));
    
    % Combine masks using max (smooth amplitude combination)
    mask_combined = max(mask1, mask2);
    
    % --- Positive version ---
    mask_img_positive = uint8(mask_combined * 255);
    filename_positive = sprintf('new_tilted_fork_mask_positive_%dum.jpg', round(Lambda*1e6));
    imwrite(mask_img_positive, filename_positive);
    
    figure;
    imshow(mask_img_positive);
    title(sprintf('Positive Mask - Grating %d µm', round(Lambda*1e6)));

    % --- Negative version ---
    mask_combined_inverted = 1 - mask_combined;
    mask_img_negative = uint8(mask_combined_inverted * 255);
    filename_negative = sprintf('new_tilted_fork_mask_negative_%dum.jpg', round(Lambda*1e6));
    imwrite(mask_img_negative, filename_negative);
    
    figure;
    imshow(mask_img_negative);
    title(sprintf('Negative Mask - Grating %d µm', round(Lambda*1e6)));
end
