%% The function that calculate the p-values
% Input: composition x, y; the elements of x and y must be positive
% Output: p-values based on clr(x) & clr(y),log(x) & log(y), and x & y
function [p_clr,p_log,p_x,clr_bacteria,naive_bacteria,log_bacteria,v_clr,v_log,v_x,n_clr,n_log,n_x] = two_sample_test(x,y)
p = size(x,2);
%% The first group 
theta_X = var(x,1,1);
 
log_X = log(x);    
theta_logX = var(log_X,1,1);

clog_TX = log (x)-1/p.*sum(log(x),2)*ones(1,p);
theta_c_X = var(clog_TX,1,1);    
%% The second group 
theta_Y = var(y,1,1);

log_Y = log(y);
theta_logY = var(log_Y,1,1);

clog_TY = log (y)-1/p.*sum(log(y),2)*ones(1,p);
theta_c_Y = var(clog_TY,1,1);    
%% Calculate the Test Statistics
nx = size(x,1);
ny = size(y,1);

[Mn_c_yx,clr_bacteria] = max((mean(clog_TY)-mean(clog_TX)).^2./(theta_c_Y./nx+theta_c_X./ny));
[Mn_yx,naive_bacteria] = max((mean(y)-mean(x)).^2./(theta_Y./nx+theta_X./ny));
[Mn_logyx,log_bacteria] = max((mean(log_Y)-mean(log_X)).^2./(theta_logY./nx+theta_logX./ny));

p_clr = 1-exp(-1/sqrt(pi)*exp(-(Mn_c_yx-(2*log(p)-log(log(p))))/2));
p_x = 1-exp(-1/sqrt(pi)*exp(-(Mn_yx-(2*log(p)-log(log(p))))/2));
p_log = 1-exp(-1/sqrt(pi)*exp(-(Mn_logyx-(2*log(p)-log(log(p))))/2));

%% Record the denominator
v_clr = theta_c_Y./nx+theta_c_X./ny;
v_log = theta_Y./nx+theta_X./ny;
v_x = theta_logY./nx+theta_logX./ny;

%% Record the nominator
n_clr = (mean(clog_TY)-mean(clog_TX)).^2;
n_log = (mean(log_Y)-mean(log_X)).^2;
n_x = (mean(y)-mean(x)).^2;