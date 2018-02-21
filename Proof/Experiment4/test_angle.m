function [] = test_angle(w,X,Y,XY,YcX,XmY)
%% Generate data
n = size(X,1);
tol = 1E-15;
Yw = Y*w;
Xw = X*w;
YcXw = YcX*w;
%% Genrate angles
m = 100;
thetavec = linspace(0,2*pi,2*m);

Col1 = zeros(numel(thetavec),1);
for id_t = 1:numel(thetavec)
    theta = thetavec(id_t);
    %den = sqrt(XmY - 2*sin(theta)*YcXw+2*(1-cos(theta))*(XY-Yw.*Xw));
    den = 1;
    %den = (XmY - 2*sin(theta)*YcXw+2*(1-cos(theta))*(XY-Yw.*Xw));
    
    % Fixed point solution
    Col1(id_t,1) = atan(sum(YcXw./den)/sum((XY-Yw.*Xw)./den));
    %Gradient
    %Col1(id_t,1) = -cos(theta)*sum(YcXw./den) + sin(theta)*sum((XY-Yw.*Xw)./den);
    %functionvalue
    %Col1(id_t,1) = sum(den);
end
%% Plot Norms
plot(thetavec,Col1);
plot(thetavec,Col1+pi);
end