t=linspace(-10,0,10000)
plot(t,1+t.*exp(t.^2/2)*sqrt(pi/2).*(1+erf(t./sqrt(2))))