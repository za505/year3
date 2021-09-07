function [y] = exponential(b,x)
%this function calculates y=A*(e^alpha*t)+y0
%where a=A, b=alpha, c=t, and y0=y0;
y=b(1)*exp(b(2)*x)+b(3);
end

