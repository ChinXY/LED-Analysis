function[par,par_se,rsq,xdat_max]=mygausianFIT(xdat,ydat)

format short e;

%ymean = mean(ydat);

ydat_max=0;
xdat_max=0;

for i=1:size(xdat,2)
if ydat(i)>ydat_max
    ydat_max=ydat(i);
    xdat_max=xdat(i);
else
end
end

x_dat1=xdat_max-150;
x_dat2=xdat_max+150;

WL1=0;WL2=0;x_int1=0;x_int2=0;
    for i=1:size(xdat,2)
        if xdat(i) < x_dat1  && WL1==0            
        elseif xdat(i) >= x_dat1 && WL1==0
            WL1=xdat(i);
            x_int1=i;
        else
        end
    end
    for i=1:size(xdat,2)
        if xdat(i) < x_dat2  && WL2==0            
        elseif xdat(i) >= x_dat2 && WL2==0
            WL2=xdat(i);
            x_int2=i;
        else
        end
    end
    if WL2==0
        WL2=xdat(size(xdat,2));
        x_int2=size(xdat,2);
    end
%if WL1 > WL2
%    WL
lb=[abs(0.5*ydat_max) WL1 10];
%ub=[2*ydat_max  WL2  (WL2-WL1)];
ub=[abs(2*ydat_max)  WL2  100];
par0(1,1)=ydat_max;
par0(1,2)=xdat_max;
par0(1,3)=30;

            [par,res,r,~,~,~,jac] = lsqcurvefit(@mygaussian,par0,xdat(1,x_int1:x_int2),ydat(1,x_int1:x_int2),lb,ub);
            ymean = mean(ydat(1,x_int1:x_int2));
            MSE = sum(r*r')/(size(xdat(x_int1:x_int2),2)-size(par0,2));  %for calculation of standard error
            Cov = inv(jac'*jac);                         %for calculation of standard error
            par_se=zeros(size(par0));                        %declaring dummy vector
            
            rsq = 1 - res/(sum((ydat-ymean).^2));
            
            for k=1:size(par0,2)                         %calculating standard error of each fitting parameter
            par_se(1,k)=sqrt(Cov(k,k)*MSE);                  %please refer to originlab website for the calculation
            end
            


end