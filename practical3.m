%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get asset returns
%  I made a function to read data for you!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   r: (T1+T2)xN asset returns (in excess of risk-free rate)
%   rf: (T1+T2)x1 risk-free returns
%   r2: T2xN out-of-sample asset returns (in excess of risk-free rate)
%   rf2: T2x1 out-of-sample risk-free returns
T1 = 120; % Estimation window size
smonth = 198201; % out-of-sample start month
emonth = 202201; % out-of-sample end month

[r, rf, r2, rf2,rben] = GetData(T1, smonth, emonth);

[T2,N] = size(r2); % out-of-sample size and num. of assets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameter estimation using sample mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Er = zeros(T2,N); % Expected returns
Cr = zeros(N,N,T2); % Covariance matrix


for t = 1:T2
    Er(t,:) = mean(r(t:t+T1-1,:));
    Cr(:,:,t) = cov(r(t:t+T1-1,:));
    alpha (t,:)= rben(t,:);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Portfolio rebalancing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_ew = zeros(T2,N); % Equal-weight portfolio


w_p = zeros(T2,N); % Optimal portfolio
w_tb = zeros(T2,N);
I = ones(N,1); % Nx1 vector of ones
for t = 1:T2
    % Er and Cr at time t
    Er_t = Er(t,:)';
    Cr_t = Cr(:,:,t);
    alpha_t = alpha(t,:);
    
    % Equal-weight portfolio
    w_ew(t,:) = 1/N*I';
    
    % Sharpe ratio maximization (Tangent portfolio)
    % Four different ways of obtaining Sharpe ratio maximizing portfolio.
    % Of course, you need only one of these.
    n=(I'*(Cr_t\Er_t))
    if n>0
        n=n;
    else
        n=0-n;
    end
    w_p(t,:) = Cr_t\Er_t/(I'*(Cr_t\Er_t)); % Better than above two
    %%%%%%
    w_tb(t,:) = alpha_t;
    
%     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Portfolio evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_ew = zeros(T2,1); % Equal-weight portfolio returns (in excess of rf)
r_p = zeros(T2,1); % Optimal portfolio returns (in excess of rf)
r_tb = zeros(T2,1);
for t = 1:T2
    r_ew(t) = w_ew(t,:)*r2(t,:)' - rf2(t);
    % WRITE CODES TO CALCULATE TANGENT PORTFOLIO RETURN HERE
    r_p(t) = w_p(t,:)*r2(t,:)' - rf2(t);
    r_tb(t) = w_tb(t,:)*r2(t,:)'-rf2(t);
end

% Mean and Stddev of returns, Sharpe ratio
Mr_ew = mean(r_ew);
Sr_ew = std(r_ew);
SR_ew = Mr_ew/Sr_ew;

% CALCULATE YOUR PORTFOLIO PERFORMANCE MEASURES HERE
Mr_p = mean(r_p);
Sr_p = std(r_p);
SR_p = Mr_p/Sr_p;

% TB performance
Mr_tb = mean(r_tb);
Sr_tb = std(r_tb);
SR_tb = Mr_tb/Sr_tb;

fprintf('              EW    Mean-Variance   T-B Model\n');
fprintf('Mean Return   %.3f      %.3f       %.3f\n', Mr_ew, Mr_p, Mr_tb);
fprintf('Std dev.      %.3f      %.3f        %.3f\n', Sr_ew, Sr_p, Sr_tb);
fprintf('Sharpe Ratio  %.3f      %.3f        %.3f\n', SR_ew, SR_p, SR_tb);







    