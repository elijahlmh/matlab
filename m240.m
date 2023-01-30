T1 = 120; % Estimation window size
smonth = 196201; % out-of-sample start month
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

    %%%trans cost
    if t==1
        r_p2(t)=r_p(t);
        r_tb2(t)=r_tb(t);
    else
        dw=w_p(t,:)-w_p(t-1,:);
        tcost=0.004*abs(dw);
        r_p2(t)=r_p(t)-sum(tcost);

        dwtb=w_tb(t,:)-w_tb(t-1,:);
        tbcost=0.004*abs(dwtb);
        r_tb2(t)=r_tb(t)-sum(tbcost);

    end

end

% Mean and Stddev of returns, Sharpe ratio
Mr_ew = mean(r_ew);
Sr_ew = std(r_ew);
SR_ew = Mr_ew/Sr_ew;

% CALCULATE YOUR PORTFOLIO PERFORMANCE MEASURES HERE
Mr_p = mean(r_p);
Sr_p = std(r_p);
SR_p = Mr_p/Sr_p;

Mr_p2 = mean(r_p2);
Sr_p2 = std(r_p2);
SR_p2 = Mr_p2/Sr_p2;

% TB performance
Mr_tb = mean(r_tb);
Sr_tb = std(r_tb);
SR_tb = Mr_tb/Sr_tb;

Mr_tb2 = mean(r_tb2);
Sr_tb2 = std(r_tb2);
SR_tb2 = Mr_tb2/Sr_tb2;

fprintf('              EW          MV        MV-ac     TB Model      TB Model-ac\n');
fprintf('Mean Return   %.3f      %.3f      %.3f      %.3f         %.3f\n',  Mr_ew, Mr_p, Mr_p2, Mr_tb, Mr_tb2);
fprintf('Std dev.      %.3f      %.3f      %.3f       %.3f          %.3f\n', Sr_ew, Sr_p, Sr_p2, Sr_tb, Sr_tb2);
fprintf('Sharpe Ratio  %.3f      %.3f      %.3f      %.3f         %.3f\n', SR_ew, SR_p, SR_p2, SR_tb, SR_tb2);