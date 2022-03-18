%% Colours

genotype_colour = {'#9e9e9e', '#fdd835', '#f44336', '#ff9800', '#2196f3', '#4caf50', '#673ab7', '#212121'};
phenotype_colour = {'#2196f3', '#f44336', '#fdd835'};

%  And now make a vector of what range each color should be at (i.e. this vector defines the spacing of the colours, they need not be regularly/equally spaced):
X = [0 0.1 0.25 0.5 0.75 0.9 1];

% Create a matrix where each row is a color triple
% V = [
%     hex2rgb('#e0e0e0') % grey 300
%     hex2rgb('#eeeeee') % grey 200
%     hex2rgb('#f5f5f5') % grey 100
%     hex2rgb('#fafafa') % grey 50
%     hex2rgb('#2196f3') % blue 500
%     hex2rgb('#1e88e5') % blue 600
%     hex2rgb('#1976d2') % blue 700
%     ];

V = [
    hex2rgb('#e0e0e0') % grey 300
    hex2rgb('#eeeeee') % grey 200
    hex2rgb('#f5f5f5') % grey 100
    hex2rgb('#fafafa') % grey 50
    hex2rgb('#bbdefb') % blue 100
    hex2rgb('#64b5f6') % blue 300
    hex2rgb('#2196f3') % blue 500
    ];

%  And finally you can create the entire map with one simple interpolation:
Xq = linspace(0, 1, 255);
Vq1 = interp1(X, V, Xq, 'pchip');


% Create a matrix where each row is a color triple
% V = [
%     hex2rgb('#e0e0e0') % grey 300
%     hex2rgb('#eeeeee') % grey 200
%     hex2rgb('#f5f5f5') % grey 100
%     hex2rgb('#fafafa') % grey 50
%     hex2rgb('#f44336') % red 500
%     hex2rgb('#e53935') % red 600
%     hex2rgb('#d32f2f') % red 700
%     ];

V = [
    hex2rgb('#e0e0e0') % grey 300
    hex2rgb('#eeeeee') % grey 200
    hex2rgb('#f5f5f5') % grey 100
    hex2rgb('#fafafa') % grey 50
    hex2rgb('#ffcdd2') % red 100
    hex2rgb('#e57373') % red 300
    hex2rgb('#f44336') % red 500
    ];

%  And finally you can create the entire map with one simple interpolation:
Xq = linspace(0, 1, 255);
Vq2 = interp1(X, V, Xq, 'pchip');

% Create a matrix where each row is a color triple
% V = [
%     hex2rgb('#e0e0e0') % grey 300
%     hex2rgb('#eeeeee') % grey 200
%     hex2rgb('#f5f5f5') % grey 100
%     hex2rgb('#fafafa') % grey 50
%     hex2rgb('#fdd835') % red 500
%     hex2rgb('#fbc02d') % red 600
%     hex2rgb('#fbc02d') % red 700
%     ];

V = [
    hex2rgb('#e0e0e0') % grey 300
    hex2rgb('#eeeeee') % grey 200
    hex2rgb('#f5f5f5') % grey 100
    hex2rgb('#fafafa') % grey 50
    hex2rgb('#fff59d') % red 200
    hex2rgb('#ffee58') % red 400
    hex2rgb('#fdd835') % red 600
    ];

%  And now make a vector of what range each color should be at (i.e. this vector defines the spacing of the colours, they need not be regularly/equally spaced):
X = [0 0.1 0.25 0.5 0.75 0.9 1];

%  And finally you can create the entire map with one simple interpolation:
Xq = linspace(0, 1, 255);
Vq3 = interp1(X, V, Xq, 'pchip');


% % Create a matrix where each row is a color triple
% V = [
%     hex2rgb('#e0e0e0') % grey 300
%     hex2rgb('#eeeeee') % grey 200
%     hex2rgb('#f5f5f5') % grey 100
%     hex2rgb('#fafafa') % grey 50
%     hex2rgb('#673ab7') % red 500
%     hex2rgb('#5e35b1') % red 600
%     hex2rgb('#512da8') % red 700
%     ];

% Create a matrix where each row is a color triple
V = [
    hex2rgb('#e0e0e0') % grey 300
    hex2rgb('#eeeeee') % grey 200
    hex2rgb('#f5f5f5') % grey 100
    hex2rgb('#fafafa') % grey 50
    hex2rgb('#d1c4e9') % red 100
    hex2rgb('#9575cd') % red 300
    hex2rgb('#673ab7') % red 500
    ];

%  And finally you can create the entire map with one simple interpolation:
Xq = linspace(0, 1, 255);
Vq4 = interp1(X, V, Xq, 'pchip');


V = zeros(3, size(Vq3, 1), size(Vq3, 2));
V(1, :, :) = Vq1;
V(2, :, :) = Vq2;
V(3, :, :) = Vq3;
V(4, :, :) = Vq4;
