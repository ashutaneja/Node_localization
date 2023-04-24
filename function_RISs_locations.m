function [ R, R1_sqrt, R2_sqrt ] = function_RISs_locations(sqrtN, d, lambda, betaH1A, betaH2A)

%% RIS Elements Coordinates

gridPoints = (0:sqrtN-1)*d;

[X,Y] = meshgrid(gridPoints,gridPoints); % coordinates expressed as matrices

locations = X(:)+1i*Y(:); % coordinates expressed as complex vector


%% RIS correlation matrix
R = sinc(2*abs(locations - transpose(locations))/lambda);
R_sqrt = sqrtm(R);


R1_sqrt = sqrt(betaH1A)*R_sqrt; %sqrtm(R1);


R2_sqrt = sqrt(betaH2A)*R_sqrt; %sqrtm(R2);
