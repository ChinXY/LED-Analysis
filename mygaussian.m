function y = mygaussian(par,xdata)    
    y = par(1)*exp(-(4*log(2))*((xdata-par(2))/(par(3))).^2);
end