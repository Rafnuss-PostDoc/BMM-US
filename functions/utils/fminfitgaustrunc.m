function [err,lambda] = fminfitgaustrunc(X,trans,x)

mu=x(1);
sigma=x(2);

f_lambda = @(lambda) ((trans.thr.^lambda-1)./lambda-norminv(trans.thr_percent,mu,sigma)).^2;

lambda = fminbnd(f_lambda,0,3);

Xt = (X.^lambda-1)/lambda;
edge = linspace((trans.thr.^lambda-1)/lambda,(3000.^lambda-1)/lambda,100);
nbins=histcounts(Xt,edge);
nbins=nbins./sum(nbins);
bin=edge(1:end-1)+diff(edge)/2;

nbinf = normpdf(bin,mu,sigma);
nbinf=nbinf/sum(nbinf);

% figure; hold on; plot(bin,nbins); plot(bin,nbinf)

err = sum((nbins-nbinf).^2);

end
