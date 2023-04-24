close all 
clear;

% This code utilises the optimization of phase shift matrix for obtaining the  values
%% in the indoor scenario%%%


lambda = 0.1;


kappa = 2*pi/lambda;

%The width and height of an RIS element
d = lambda/4;

% Maximum mumber of elements per dimension
MaxNumOfElements = 25;

%Number of channel realizations
numOfChan = 1000;

% Vector of number of elements per dimension
NumOfElements = 1:2:MaxNumOfElements;


% Bandidth
Bandwidth = 20e3;

% RIS element area
A = d.^2;

% Total Radiated Power in dBm
PowerdBm = 23;

% Total Radiated Power in mWatt
Power = db2pow(PowerdBm);

% Thermal noise in dBm
SigmaW2dBm = pow2db(Bandwidth)-174; % -94 dBm;

% Thermal noise in Watt (sigma2_w)
SigmaW2 = db2pow(SigmaW2dBm);

% SNR
Psigma2 = Power/SigmaW2;

% SNR in dB
Psigma2dB = pow2db(Psigma2);

% Channel gain h1
betaH1A = db2pow(-48)*A;

% Channel gain h2
betaH2A = db2pow(-38)*A;

% Channel gain of direct link
betaHd = db2pow(-inf);


rho = 2;


Sigma2dBm = PowerdBm + pow2db(betaH1A/A) - rho;


Sigma2A = db2pow(Sigma2dBm)*A;



% Preallocation
x_coordinate = zeros(1,MaxNumOfElements^2);
y_coordinate = zeros(1,MaxNumOfElements^2);

x_coordinate(1) = 0;
y_coordinate(1) = 0;

for n = 2:MaxNumOfElements^2
    x_coordinate(n) = x_coordinate(n-1) + sin(mod(floor(sqrt(4*(n-2)+1)),4)*pi/2);
    y_coordinate(n) = y_coordinate(n-1) + cos(mod(floor(sqrt(4*(n-2)+1)),4)*pi/2);
end


x_coordinate = x_coordinate-min(x_coordinate);
y_coordinate = y_coordinate-min(y_coordinate);


x_coordinate=x_coordinate*d;
y_coordinate=y_coordinate*d;


figure, plot(x_coordinate,y_coordinate,'-')
xlim([-2,MaxNumOfElements+1]*d)
ylim([-2,MaxNumOfElements+1]*d)


theta_angle_mean = pi/4;
theta_angle_variance_degrees = 20;
theta_angle_variance = (theta_angle_variance_degrees*pi/180)^2;
phi_angle_mean = 0;
phi_angle_variance = theta_angle_variance;


R1 = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
    phi_angle_mean,phi_angle_variance,kappa);


theta_angle_mean = 0;
theta_angle_variance_degrees = 20;
theta_angle_variance = (theta_angle_variance_degrees*pi/180)^2;
phi_angle_mean = pi/4;
phi_angle_variance = theta_angle_variance;


R2 = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
    phi_angle_mean,phi_angle_variance,kappa);


variances_vector = ([45 Inf]*pi/180).^2;


meanSNR = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_algorithm = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_rayleigh = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_noEMI = zeros(1,numel(NumOfElements));

for variances_index = 1:numel(variances_vector)
    
    variances = variances_vector(variances_index);
    
   
    theta_angle_mean = 0;
    theta_angle_variance = variances;
    
    
    phi_angle_mean = -pi/4;
    phi_angle_variance = variances;
    
    
    RN = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
        phi_angle_mean,phi_angle_variance,kappa);

    [ Rn, R1_sqrt, R2_sqrt ] = function_RISs_locations(sqrtN, d, lambda, betaH1A, betaH2A);
     
    
    % Prepare to save results
    SNR = zeros(numOfChan,1);
    SNR_noEMI = zeros(numOfChan,1);
    SNR_algorithm = zeros(numOfChan,1);
    SNR_rayleigh = zeros(numOfChan,1);
    
    % Loop over number of RIS elements per dimension
    for ii = 1:numel(NumOfElements)
        
        N = NumOfElements(ii)^2;
        
        disp(['N: ',num2str(N)])
        
        % Reduced-size R1 sqrt matrix
        R1sq_i = sqrtm(R1(1:N,1:N));
        
        % Reduced-size R2 sqrt matrix
        R2sq_i = sqrtm(R2(1:N,1:N));
        
        % Reduced-size Rn matrix
        RN_i = RN(1:N,1:N);
        
        % Loop over channel realizations
        parfor kk = 1:numOfChan
            
           
            h1 = sqrt(betaH1A)*R1sq_i*sqrt(.5)*(randn(N,1) + 1j*randn(N,1)); 
            h2 = sqrt(betaH2A)*R2sq_i*sqrt(.5)*(randn(N,1) + 1j*randn(N,1)); 
            hd = sqrt(betaHd)*sqrt(.5)*(randn(1,1) + 1j*randn(1,1));           
            
            
            theta = diag((exp(1j*(angle(conj(h2).*h1)-angle(hd)))));
            
           
            g2 = theta*h2;
            
           
                  
           [SNR_algorithm(kk),~,SNR_rayleigh(kk),~] = function_optimization(Power, SigmaW2, h2, h1, hd, Sigma2A, RN_i,0);
  
            
        end
        
        
        meanSNR(variances_index,ii) = mean(real(SNR));
        meanSNR_algorithm(variances_index,ii) = mean(real(SNR_algorithm));
        meanSNR_rayleigh(variances_index,ii) = mean(real(SNR_rayleigh));
        
        if variances_index == 1
            
            meanSNR_noEMI(1,ii) = mean(real(SNR_noEMI));
            
        end
        
    end
    
end



% Plot numerical results

hold on; box on; grid on
set(gca,'fontsize',16);
plot(NumOfElements.^2,smooth(meanSNR(2,:)),'LineWidth',2,'LineStyle',':','Color',[0 0 1]),
plot(NumOfElements.^2,smooth(meanSNR_algorithm(2,:)),'LineWidth',2,'LineStyle','-.','Color',[1 0 0]),

set(gca,'fontsize',18);
legend({'Opt. thermal noise',...
    'Iterative Algorithm',...
    ...
   },'Interpreter','latex','Location','best')


