function y=virusobserve(x,H,obsnum,GA,gamma)
if strcmp(GA,'Gaussian')
    y=H*x+gamma.*randn(obsnum,1);
elseif strcmp(GA,'Normal')
    y=H*x+gamma.*(2*rand(obsnum,1)-1);
else
    y=H*x;
end
