%% The function that calculate the p-values
% Input: composition x, y; the elements of x and y must be positive
% Output: p-values based on clr(x) & clr(y),log(x) & log(y), and x & y
function [p_clr,p_log,p_x,clr_bacteria,naive_bacteria,log_bacteria,v_clr,v_log,v_x,n_clr,n_log,n_x] = pair_test(x,y)
[n,p] = size(x);
%% The difference of clr
sampleDiff = x-y;
theta_sampleDiff = var(sampleDiff,1,1);

sampleRatio = x./y;
log_sampleRatio = log(sampleRatio);    
theta_logsampleRatio = var(log_sampleRatio,1,1);
    
log_sampleRatio_prod = sum(log(sampleRatio),2);
clog_TsampleRatio = log (sampleRatio)-1/p.*log_sampleRatio_prod*ones(1,p);
theta_c_sampleRatio = var(clog_TsampleRatio,1,1);    
    
%% Calculate the Test Statistics
[Mn_clr,clr_bacteria] = max((mean(clog_TsampleRatio)).^2./(theta_c_sampleRatio./n));
[Mn_naive,naive_bacteria] = max((mean(sampleDiff)).^2./(theta_sampleDiff./n));
[Mn_log,log_bacteria] = max((mean(log_sampleRatio)).^2./(theta_logsampleRatio./n));

p_clr = 1-exp(-1/sqrt(pi)*exp(-(Mn_clr-(2*log(p)-log(log(p))))/2));
p_x = 1-exp(-1/sqrt(pi)*exp(-(Mn_naive-(2*log(p)-log(log(p))))/2));
p_log = 1-exp(-1/sqrt(pi)*exp(-(Mn_log-(2*log(p)-log(log(p))))/2));

%% Record the denominator
v_clr = theta_c_sampleRatio./n;
v_log = theta_logsampleRatio./n;
v_x = theta_sampleDiff./n;
%% Record the nominator
n_clr = (mean(clog_TsampleRatio)).^2;
n_log = (mean(log_sampleRatio)).^2;
n_x = (mean(sampleDiff)).^2;