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

% A single variational Bayes measurement update.
% Input: 
%   xkk_1i:     Predicted kinematic state mean vector
%   Pkk_1i:     Predicted kinematic state covariance matrix
%   betakk_1:   Predicted Extent scale parameter beta
%   alphakk_1:  Predicted Extent shape parameter alpha
%   yk:         Measurement set at time k
%   H:          Mesurement matrix
%   R:          Measurement covariance matrix
%   s:          Scaling constant s
% Output: 
%   xk_o:      Updated kinematic state mean vector
%   Pk_o:      Updated kinematic state covariance matrix
%   betak_o:   Updated Extent scale parameter beta
%   alphak_o:  Updated Extent shape parameter alpha


function [xk_o,Pk_o,betak_o,alphak_o] = UpdateVB(xkk_1i,Pkk_1i,betakk_1,alphakk_1,yk,H,R,s)
Nmeas = size(yk,2); % Number of measurements
Pkk_1 = Pkk_1i(1:4,1:4); % Position and velocity components
xkk_1 = xkk_1i(1:4); % Position and velocity components
SigmaThetakk_1 = Pkk_1i(5,5); % Seperate orientation from the kinematic state
Thetakk_1 = xkk_1i(5); % Seperate orientation from the kinematic state

%% Set Initials for Variational Iterations and Compute Inverses
betak = betakk_1;
alphak = alphakk_1;
Pk = Pkk_1;
iPkk_1 = Pkk_1^-1;
xk = xkk_1;
Thetak = Thetakk_1;
iSigmaThetakk_1 = SigmaThetakk_1^-1;
SigmaThetak = SigmaThetakk_1;
zk = yk;
SigmaZk = s*betakk_1./(alphakk_1-1);
SigmaZk = diag(SigmaZk);
iR = R^-1;

% Define rotation matrix and its derivative
L = @(tk)[cos(tk) -sin(tk); 
          sin(tk)  cos(tk)];  % Rotation matrix
Lprime = @(tk)[-sin(tk) -cos(tk);
               cos(tk)  -sin(tk)]; % Derivative of the rotation matrix

Hbar = H(:,1:4);
NumOfIterations = 10;
for i=1:NumOfIterations % Variational Iterations
    
    %compute expectations
    %E{x}
    xbar = xk;
    %E{theta}
    Thetabar = Thetak;
    %E{sX^-1}
    Xbar=diag(s*betak./(alphak));
    Xinvbar=diag(diag(Xbar).^-1);
    %E{z}
    zjbar = zk;
    %sum(E{z})
    zbar = sum(zjbar,2)./Nmeas;
    if Nmeas == 0
        zbar = zeros(2,1);
    end
    %E{L(theta)*sX^-1*L(theta)^T}
    LXLbar=ComputeLXL(Thetabar,SigmaThetak, Xinvbar);
    %E{(zjbar-Hxbar)(zjbar-Hxbar)'}
    for mt = 1:Nmeas
        innobar(:,:,mt) =(zjbar(:,mt)-Hbar*xbar)*(zjbar(:,mt)-Hbar*xbar)' + Hbar*Pk*Hbar' + SigmaZk;
    end
    
    %compute qx
    Pk=(iPkk_1 + Nmeas*Hbar'*LXLbar*Hbar)^-1;
    xk=Pk*(iPkk_1*xkk_1 + Nmeas*Hbar'*LXLbar*zbar);
    
    %compute qX
    alphak = alphakk_1 + Nmeas*0.5;    
    LsumLbar = ComputeLXL(-Thetabar,SigmaThetak, (1/(2*s))*sum(innobar,3));
    betak = betakk_1 + diag(LsumLbar);
    
    %compute qZ
    SigmaZk=(LXLbar + iR)^-1;
    for mt = 1:Nmeas
        zk(:,mt) = SigmaZk*(LXLbar*Hbar*xbar + iR*yk(:,mt));
    end
    
    %compute qTetha
    Sigmatemp=0;
    innotemp=0;
    for mt=1:Nmeas
        M=trace(Xinvbar*Lprime(Thetabar)'*innobar(:,:,mt)*Lprime(Thetabar));
        m=trace(Xinvbar*Lprime(Thetabar)'*innobar(:,:,mt)*Lprime(Thetabar)*Thetabar)-...
            trace(Xinvbar*L(Thetabar)'*innobar(:,:,mt)*Lprime(Thetabar));
        Sigmatemp=Sigmatemp+M;
        innotemp= innotemp+m ;
    end
    SigmaThetak = (iSigmaThetakk_1 + Sigmatemp)^-1;
    Thetak = SigmaThetak*(iSigmaThetakk_1*Thetakk_1 + innotemp);
end

% Set outputs
xk_o = [xk;Thetak];
Pk_o = blkdiag(Pk,SigmaThetak);
alphak_o = alphak;
betak_o = betak;

end

function LXLbar=ComputeLXL(Thetabar,SigmaThetak,X)

    %E{L(theta)^T*X*L(theta)}
    a = X(1,1); b = X(1,2); c = X(2,1);  d = X(2,2);
    LXLbar(1,1) = a*(1+cos(2*Thetabar)*exp(-2*SigmaThetak)) + d*(1-cos(2*Thetabar)*exp(-2*SigmaThetak))- (c+b)*(sin(2*Thetabar)*exp(-2*SigmaThetak));
    LXLbar(1,2) = b*(1+cos(2*Thetabar)*exp(-2*SigmaThetak)) - c*(1-cos(2*Thetabar)*exp(-2*SigmaThetak))+ (a-d)*(sin(2*Thetabar)*exp(-2*SigmaThetak));
    LXLbar(2,1) = c*(1+cos(2*Thetabar)*exp(-2*SigmaThetak)) - b*(1-cos(2*Thetabar)*exp(-2*SigmaThetak))+ (a-d)*(sin(2*Thetabar)*exp(-2*SigmaThetak));
    LXLbar(2,2) = d*(1+cos(2*Thetabar)*exp(-2*SigmaThetak)) + a*(1-cos(2*Thetabar)*exp(-2*SigmaThetak))+ (c+b)*(sin(2*Thetabar)*exp(-2*SigmaThetak));
    
    LXLbar=0.5.*LXLbar;

end