%***************************************************************************
%** Implemenetation of the VB algorithm based on the journal article      **
%** "Random Matrix Based Extended Target Tracking with Orientation: A New **
%**  Model and Inference" DOI: 10.1109/TSP.2021.3065136                   **
%** Barkýn Tuncer and Emre Özkan                                          **
%** Further information:                                                  **
%** https://github.com/Metu-Sensor-Fusion-Lab                             **
%** http://sensorfusion.eee.metu.edu.tr/                                  **
%**                                                                       **
%** Source code written by Barkýn Tuncer(tuncer.barkin@gmail.com)         **
%***************************************************************************

close all
clc
clear
dbstop warning
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultFigureUnits','normalized')
pos = [ 0.1 0.1 0.7 0.6];
set(0, 'defaultFigureUnits', 'normalized','defaultFigurePosition', pos)
%% Experiment Parameters
Tfinal = 105; % Number of frames
T = 1; % Sampling time
rng(94); % For reproducibility
%% VB Parameter Initialization
m.H = [1 0 0 0 0; 0 1 0 0 0]; % Measurement transition matrix
m.R = diag([1,1]); % Measurement noise covariance matrix
sigmaQ = 4; % Standard deviation of the noise
m.Q_ = @(dt)kron([dt^3/3 dt^2/2; dt^2/2 dt], diag([sigmaQ^2 sigmaQ^2]));
m.Q = blkdiag(m.Q_(T),0.1); % Process noise covariance matrix
m.F = blkdiag([eye(2),T*eye(2); zeros(2,2),eye(2)],1); % State transition matrix
xk_1 = [450 245 0 0 -pi/2]'; % Initial augmented kinematic state mean [Pos_x Pos_y Vel_x Vel_y Orientation]
Pk_1 = blkdiag(eye(4),1); % Initial augmented kinematic state covariance
alphak_1 = [2 2]'; % Initial alpha for inverse Gamma distribution
betak_1 = [1000 250]'.*(alphak_1-1); % Initial beta for inverse Gamma distribution
m.forgettingfactor = 0.99; %Forgetting factor 
s = 1/4; % Scaling constant
visualize = 1; % Should we visualize or not?
%% Measurement Reading, Time & Measurement update
load('measurements.mat'); % Load measurement set
load('images.mat'); % Load image set
for t = 1:Tfinal
    %% Read Measurements
    y = Z{t}; % Read measurements
    img = I_(:,:,:,t); % Read image
    disp(['Time step:' num2str(t) ' , ' num2str(size(y,2)) ' measurements']);
    %% VB Time Update
    if t == 1
        xkk_1 = xk_1;
        Pkk_1 = Pk_1;
        betakk_1 = betak_1;
        alphakk_1 = alphak_1;
    else
        xkk_1 = m.F*xk_1;
        Pkk_1 = m.F*Pk_1*m.F' + m.Q;
        betakk_1 = m.forgettingfactor*betak_1;
        alphakk_1 = m.forgettingfactor*alphak_1;
    end
    %% VB Measurement Update
    [xk,Pk,betak,alphak] = UpdateVB(xkk_1,Pkk_1,betakk_1,alphakk_1,y,m.H,m.R,s);
    xk_1 = xk;    Pk_1 = Pk;     betak_1 = betak;     alphak_1 = alphak;
    len_VB = (betak_1./(alphak_1-1)).^0.5; % Get ellipse semi axis lengths
    ang_VB = xk_1(5); % Get the orientation
    VB_par = [m.H*xk_1; len_VB; ang_VB];
    %% Visualization
    if visualize
        cla
        imshow(img);
        hold all
        axis equal
        meas_points = plot(y(1,:), y(2,:), '*g','LineWidth',4,'MarkerSize',1);
        est_plot_VB = plot_ellipse([VB_par(1:2); VB_par(3:4); VB_par(5)],'-','k',2.5);
        drawnow;
%         pause(0.1);
    end
end
