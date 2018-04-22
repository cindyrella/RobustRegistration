n = 1:10;
lambda = linspace(0,pi,100);
lambda = lambda(2:99);

[N,Lambda]= meshgrid(n,lambda);

P = [0,5,10,25,40,50,70]/100

%  for p = P
%  bound = sqrt(N).*(p/(1-p).*1./(1+cos(Lambda))+(1+sqrt(pi./(2*N))).^2)./max(1./sqrt(1+cos(Lambda)),1+sqrt(pi./(2.*N)));
%  figure()
%  contour(N,Lambda,min(bound,ones(size(bound))))
%  end
 
 
A = log([0.15,0.1,0.05,0.01]);
for alpha = A
for p = P
    sN = (log(N)-alpha)./(1-p).^2.*32.*N.*(p./(1+cos(Lambda))+(1-p).*(1+sqrt(pi./(2.*N))).^2);
    sN = log10(max(sN,ones(size(sN))));
    figure()
    contour(N,Lambda,sN)
end

end
