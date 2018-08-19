function [par0,spec_Lum_max,data,par1,par2,par3,data_color,all_spectra_cell,data_peak]=LED_analysis8(folder,file,file_spec,fileM,WL_limit,Vth_def)
    javaaddpath('poi_library/poi-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('poi_library/dom4j-1.6.1.jar');
    javaaddpath('poi_library/stax-api-1.0.1.jar');

    format short e;
    filename = [folder '\' file];
    filename_spec = [folder '\' file_spec];
    FileName_spec=filename_spec;
    file_meta = [folder '\' fileM];
    fid=fopen(file_meta,'r');
    metadata=textscan(fid,'%s','HeaderLines',0);
    fclose(fid);
    A=metadata{1,1}{21,1};
    B=A(1,6:end-6);
    dev_area=str2double(B);
    clear A B metadata
    V_rev=-1.5; % Reverse current density extracted at V_rev
    
    %find integration boundaries
    x1=WL_limit(1,1);
    x2=WL_limit(1,2);
    fid = fopen(FileName_spec, 'r');
    PLdata = textscan(fid,'%f %f %f','HeaderLines',22);
    PLwavelength = PLdata{1,1};
    fclose(fid);
    WL=PLwavelength(:,1); % in nm
    WL1=0;WL2=0;x_int1=0;x_int2=0;
    for i=1:size(WL,1)
        if WL(i,1) < x1  && WL1==0            
        elseif WL(i,1) >= x1 && WL1==0
            WL1=WL(i,1);
            x_int1=i;
        else
        end
    end
    for i=1:size(WL,1)
        if WL(i,1) < x2  && WL2==0            
        elseif WL(i,1) >= x2 && WL2==0
            WL2=WL(i,1);
            x_int2=i;
        else
        end
    end
    
    disp(' ');
    disp(' ');
    display(['For data file ' file]);
    display(['and spectra file ' file_spec]);
    disp(['device area of ' num2str(dev_area) ' m2']);
    disp(['integrating from ' num2str(WL1) ' nm to ' num2str(WL2) ' nm']);
    disp(['Vth_def = ' num2str(Vth_def) ' V']);
    
    FileName_data=filename;
    fid = fopen(FileName_data, 'r');
    JVLdata = textscan(fid,'%f %f %f %f %f','HeaderLines',1);
    V = JVLdata{1,1}; % un V
    J = JVLdata{1,2}; % in A
    Lum_int = JVLdata{1,3}; % in cd
    Power_out= JVLdata{1,5}; % uW??
    fclose(fid);
    
    %searching for index of V_rev
    x_Vrev=1;
    Vrev=0;
        for i=1:size(V,1)
            if V(i,1) < V_rev  && x_Vrev==1            
            elseif V(i,1) >= V_rev && x_Vrev==1
            Vrev=V(i,1);
            x_Vrev=i;
            else
            end
        end
    disp(['Extracting Jrev at V = ' num2str(Vrev) ' V']);
    
    Luminance=Lum_int/dev_area; % cd/m2
    Jden_cm2=J*1e-4*1e3/dev_area; %in mA/cm2
    Lum_eff=Lum_int./J; % in Cd/A
    Lum_pow_eff=pi*Lum_eff./V; % in lm/W
    Power_in=V.*J; % input power in W
    n_wp=100*abs(Power_out*1e-6)./abs(Power_in); % wall plug efficiency in percentage
    Photon_count=zeros(size(J));
    EQE=zeros(size(J));
    Jrev=Jden_cm2(x_Vrev,1);
    data2=zeros(size(V,1),9);
    data3=zeros(size(V,1),8);

           itpl_WL = x1:1:x2;
           mat_itpl_abs_irr         = zeros(3+size(itpl_WL,2),1+size(V,1));
           mat_itpl_abs_irr(4:size(mat_itpl_abs_irr,1),1) = itpl_WL';
           
           mat_norm_itpl_abs_irr    = zeros(3+size(itpl_WL,2),1+size(V,1));
           mat_norm_itpl_abs_irr(4:size(mat_norm_itpl_abs_irr,1),1) = itpl_WL';
           
           mat_abs_irr              = zeros(3+size(WL(x_int1:x_int2),1),1+size(V,1));
           mat_abs_irr (4:size(mat_abs_irr,1),1) = WL(x_int1:x_int2);
           
           mat_norm_abs_irr         = zeros(3+size(WL(x_int1:x_int2),1),1+size(V,1));
           mat_norm_abs_irr(4:size(mat_norm_abs_irr,1),1) = WL(x_int1:x_int2);
           
           all_spectra_cell=cell([2,4]);

   
    % calculating EQE
    for i=1:size(V,1)
        headerlinesIn=22+(i-1)*1054;
        fid = fopen(FileName_spec, 'r');
        PLdata = textscan(fid,'%f %f %f','HeaderLines',headerlinesIn);
        PLwavelength = PLdata{1,1};
        PLIntensity = PLdata{1,3};
        fclose(fid);
        display(['calculating EQE ' num2str(i) ' / ' num2str(size(V,1)) ]);
        wavelength=PLwavelength(x_int1:x_int2,1); % in nm
        Energy_wavelength=(1./wavelength)*6.62607004e-34*299792458*1e9; % in Joule
        abs_irradiate=PLIntensity(x_int1:x_int2,1); % in uW/nm
        
        [xy31,xyz31,uv1960,cct,duv]=CIE_xy(wavelength,abs_irradiate);
        data2(i,:)=[xy31(1,1),xy31(1,2),xyz31(1,1),xyz31(1,2),xyz31(1,3),uv1960(1,1),uv1960(1,2),cct,duv];

        itpl_abs_irr = interp1(wavelength,abs_irradiate,itpl_WL','PCHIP',0);
        norm_itpl_abs_irr = itpl_abs_irr/max(itpl_abs_irr);
        norm_abs_irr = abs_irradiate/max(abs_irradiate);

        mat_itpl_abs_irr(1,i+1) = V(i,1);
        mat_itpl_abs_irr(2,i+1) = abs(Jden_cm2(i,1));
        mat_itpl_abs_irr(3,i+1) = Luminance(i,1);
        mat_itpl_abs_irr(4:4+size(itpl_WL,2)-1,i+1) = itpl_abs_irr(:,1);

        mat_abs_irr(1,i+1) = V(i,1);
        mat_abs_irr(2,i+1) = abs(Jden_cm2(i,1));
        mat_abs_irr(3,i+1) = Luminance(i,1);
        mat_abs_irr(4:4+size(wavelength,1)-1,i+1) = abs_irradiate(:,1);

        mat_norm_itpl_abs_irr(1,i+1) = V(i,1);
        mat_norm_itpl_abs_irr(2,i+1) = abs(Jden_cm2(i,1));
        mat_norm_itpl_abs_irr(3,i+1) = Luminance(i,1);
        mat_norm_itpl_abs_irr(4:4+size(itpl_WL,2)-1,i+1) = norm_itpl_abs_irr(:,1);

        mat_norm_abs_irr(1,i+1) = V(i,1);
        mat_norm_abs_irr(2,i+1) = abs(Jden_cm2(i,1));
        mat_norm_abs_irr(3,i+1) = Luminance(i,1);
        mat_norm_abs_irr(4:4+size(wavelength,1)-1,i+1) = norm_abs_irr(:,1);
            
        Photon_count(i,1)=trapz(wavelength(:,1),abs_irradiate(:,1)./Energy_wavelength(:,1));%*1e-6;
        EQE(i,1)=(1/3.152)*100*Photon_count(i,1)/(J(i,1)/1.60217662e-19); % in percentage
        
        abs_irr=abs_irradiate/(1e-3*max(abs_irradiate));
        [par,par_se,rsq,xdat_max]= mygausianFIT(wavelength',abs_irr');
        amplitude       = par(1)*(1e-3*max(abs_irradiate));
        amplitude_se    =par_se(1)*(1e-3*max(abs_irradiate));
        peak_center     =par(2);
        peak_center_se  =par_se(2);
        FWHM            =par(3);
        FWHM_se         =par_se(3);
        R2              =rsq;
        peak_center_t   =xdat_max;
        data3(i,:) = [amplitude, amplitude_se, peak_center, peak_center_se, FWHM, FWHM_se, R2, peak_center_t];
      
    end
    
        all_spectra_cell(1,1) = {'mat_itpl_abs_irr'};
        all_spectra_cell(1,2)={'mat_abs_irr'};
        all_spectra_cell(1,3)={'mat_norm_itpl_abs_irr'};
        all_spectra_cell(1,4)={'mat_norm_abs_irr'};
        all_spectra_cell{2,1}=mat_itpl_abs_irr;
        all_spectra_cell{2,2}=mat_abs_irr;
        all_spectra_cell{2,3}=mat_norm_itpl_abs_irr;
        all_spectra_cell{2,4}=mat_norm_abs_irr;
        xlwrite([folder '\' file_spec '_all_spectra.xlsx'],mat_itpl_abs_irr,1);
        xlwrite([folder '\' file_spec '_all_spectra.xlsx'],mat_abs_irr,2);
        xlwrite([folder '\' file_spec '_all_spectra.xlsx'],mat_norm_itpl_abs_irr,3);
        xlwrite([folder '\' file_spec '_all_spectra.xlsx'],mat_norm_abs_irr,4);
    % calculating turn on voltage
    V_on_1=0;
    Jth=0;
    for i=1:size(V,1)
        if V(i,1)<Vth_def
        else
            if Luminance(i,1)<1 && V_on_1==0            
            elseif Luminance(i,1)>=1 && V_on_1==0
                V_on_1=V(i,1);
                Jth=Jden_cm2(i,1);
            else
            end
        end
    end
            disp(' ');
            disp('At Luminance 1Cd/m2');
            disp(['voltage = ' num2str(V_on_1) ' V']);
            Vth=V_on_1;
    
    if Vth == 0
        Jth= 0;
        Lum_max= 0;
        V_Luminance_max= 0;
        CurEff_max= 0;
        V_Lum_eff_max= 0;
        Lum_pow_max= 0;
        V_Lum_pow_max_0= 0;
        Jrev= 0;
        EQE_max= 0;
        V_EQE_max_0 = 0;
        V_on_100= 0;
        Jden_100= 0;
        Lum_100= 0;
        EQE_100= 0;
        Lum_eff_100= 0;
        Lum_pow_eff_100= 0;
        V_on_1000= 0;
        Jden_1000= 0;
        Lum_1000= 0;
        EQE_1000= 0;
        Lum_eff_1000= 0;
        Lum_pow_eff_1000= 0;
    	xy31(1,1)= 0;
        xy31(1,2)= 0;
        xyz31(1,1)= 0;
        xyz31(1,2)= 0;
        xyz31(1,3)= 0;
        uv1960(1,1)= 0;
        uv1960(1,2)= 0;
        cct= 0;
        duv= 0;
        x_spec =1;
        headerlinesIn=22+(x_spec-1)*1054;
        fid = fopen(FileName_spec, 'r');
            MaxEL = textscan(fid,'%f %f %f','HeaderLines',headerlinesIn);
            EL_wavelength = MaxEL{1,1};
            EL_counts = MaxEL{1,2};
            EL_abs_irr = MaxEL{1,3};
            spec_Lum_max0(:,1) = EL_wavelength(x_int1:x_int2,1);
            spec_Lum_max0(:,2) = EL_counts(x_int1:x_int2,1);
            spec_Lum_max0(:,3) = EL_abs_irr(x_int1:x_int2,1);
            spec_Lum_max0(:,4) = EL_abs_irr(x_int1:x_int2,1)/max(EL_abs_irr(x_int1:x_int2,1));
        fclose(fid);
            spec_Lum_max=zeros(1+size(spec_Lum_max0,1),size(spec_Lum_max0,2));
            spec_Lum_max(2:size(spec_Lum_max,1),1) = spec_Lum_max0(:,1);
            spec_Lum_max(2:size(spec_Lum_max,1),2) = spec_Lum_max0(:,2);
            spec_Lum_max(2:size(spec_Lum_max,1),3) = spec_Lum_max0(:,3);
            spec_Lum_max(2:size(spec_Lum_max,1),4) = spec_Lum_max0(:,4);
            spec_Lum_max=num2cell(spec_Lum_max);
            spec_Lum_max{1,1} = 'Wavelength (nm)';
            spec_Lum_max{1,2} = 'Counts';
            spec_Lum_max{1,3} = 'Absolute Irradiance (W/nm)';
            spec_Lum_max{1,4} = 'Normalized Absolute Irradiance';
            par3=zeros(1,8);

    else    
    %calculating maximum wall plug efficiency
    n_wp_max=0;
    for i=1:size(V,1)
        if n_wp(i,1)> n_wp_max && V(i,1) >= V_on_1 
           n_wp_max=n_wp(i,1);
        else
        end
    end

    %calculating maximum luminance
    Luminance_max=0;    
    V_Luminance_max=0;
    x_spec=1;
    for i=1:size(V,1)
        if Luminance(i,1)> Luminance_max && V(i,1) >= V_on_1 
           Luminance_max=Luminance(i,1);
           V_Luminance_max=V(i,1);
           x_spec=i;
        else
        end
    end
    headerlinesIn=22+(x_spec-1)*1054;
    fid = fopen(FileName_spec, 'r');
        MaxEL = textscan(fid,'%f %f %f','HeaderLines',headerlinesIn);
        EL_wavelength = MaxEL{1,1};
        EL_counts = MaxEL{1,2};
        EL_abs_irr = MaxEL{1,3};
        spec_Lum_max0(:,1) = EL_wavelength(x_int1:x_int2,1);
        spec_Lum_max0(:,2) = EL_counts(x_int1:x_int2,1);
        spec_Lum_max0(:,3) = EL_abs_irr(x_int1:x_int2,1);
        spec_Lum_max0(:,4) = EL_abs_irr(x_int1:x_int2,1)/max(EL_abs_irr(x_int1:x_int2,1));
    fclose(fid);

    % Calculating the color coordinate of EL spectrum at maximum luminance
    WL_input=spec_Lum_max0(:,1);
    int_input=spec_Lum_max0(:,3);
    %[xy31,xyz31]=CIE_xy(WL_input,int_input);
    [xy31,xyz31,uv1960,cct,duv]=CIE_xy(WL_input,int_input);
    abs_irr=int_input/(1e-3*max(int_input));
    
    [par,par_se,rsq,xdat_max]= mygausianFIT(WL_input',abs_irr');
        max_Lum_amplitude=par(1)*(1e-3*max(int_input));
        max_Lum_amplitude_se=par_se(1)*(1e-3*max(int_input));
        max_Lum_peak_center=par(2);
        max_Lum_peak_center_se=par_se(2);
        max_Lum_FWHM=par(3);
        max_Lum_FWHM_se=par_se(3);
        max_Lum_R2=rsq;
        max_Lum_peak_center_t=xdat_max;
    par3=[max_Lum_amplitude, max_Lum_amplitude_se, max_Lum_peak_center, max_Lum_peak_center_se, max_Lum_FWHM, max_Lum_FWHM_se, max_Lum_R2,max_Lum_peak_center_t]; %8 element
    
    disp(' ');
    disp(['max luminance ' num2str(Luminance_max) 'Cd/m2']);
    Lum_max=Luminance_max;
    disp(['voltage @ max luminance ' num2str(V_Luminance_max) ' V']);
    disp(' ');

    spec_Lum_max=zeros(1+size(spec_Lum_max0,1),size(spec_Lum_max0,2));
    spec_Lum_max(2:size(spec_Lum_max,1),1) = spec_Lum_max0(:,1);
    spec_Lum_max(2:size(spec_Lum_max,1),2) = spec_Lum_max0(:,2);
    spec_Lum_max(2:size(spec_Lum_max,1),3) = spec_Lum_max0(:,3);
    spec_Lum_max(2:size(spec_Lum_max,1),4) = spec_Lum_max0(:,4);
    spec_Lum_max=num2cell(spec_Lum_max);
    spec_Lum_max{1,1} = 'Wavelength (nm)';
    spec_Lum_max{1,2} = 'Counts';
    spec_Lum_max{1,3} = 'Absolute Irradiance (W/nm)';
    spec_Lum_max{1,4} = 'Normalized Absolute Irradiance';
    
    %Extracting parameter at 100 cd/m2
    V_on_100=0;
    Jden_100=0;
    EQE_100=0;
    Lum_eff_100=0;
    Lum_pow_eff_100=0;
    Lum_100=0;
    for i=1:size(Luminance,1)
        if Luminance(i,1) < 100  && Lum_100==0 && V(i,1) >= V_on_1            
        elseif Luminance(i,1) >= 100 && Lum_100==0 && V(i,1) >= V_on_1
            V_on_100=V(i,1);
            Jden_100=Jden_cm2(i,1);
            Lum_100=Luminance(i,1);
            EQE_100=EQE(i,1);
            Lum_eff_100=Lum_eff(i,1);
            Lum_pow_eff_100=Lum_pow_eff(i,1);
        else
        end
    end
    
    %Extracting parameter at 1000 cd/m2
    V_on_1000=0;
    Jden_1000=0;
    EQE_1000=0;
    Lum_eff_1000=0;
    Lum_pow_eff_1000=0;
    Lum_1000=0;
    for i=1:size(Luminance,1)
        if Luminance(i,1) < 1000  && Lum_1000==0 && V(i,1) >= V_on_1           
        elseif Luminance(i,1) >= 1000 && Lum_1000==0 && V(i,1) >= V_on_1
            V_on_1000=V(i,1);
            Jden_1000=Jden_cm2(i,1);
            Lum_1000=Luminance(i,1);
            EQE_1000=EQE(i,1);
            Lum_eff_1000=Lum_eff(i,1);
            Lum_pow_eff_1000=Lum_pow_eff(i,1);
        else
        end
    end
    
    % Calculating maximum current efficiency
    Lum_eff_max=0;
    V_Lum_eff_max=0;
    for i=1:size(V,1)
        if Lum_eff(i,1)> Lum_eff_max && V(i,1) >= V_on_1 
           Lum_eff_max=Lum_eff(i,1);
           V_Lum_eff_max=V(i,1);
        else
        end
    end
    CurEff_max=Lum_eff_max;
    
    % Calculating maximum Luminous power efficiency
    Lum_pow_max_0=0;
    V_Lum_pow_max_0=0;
        for i=1:size(V,1)
            if Lum_pow_eff(i,1)> Lum_pow_max_0 && V(i,1) >= V_on_1 
            Lum_pow_max_0=Lum_pow_eff(i,1);
            V_Lum_pow_max_0=V(i,1);
            else
            end
        end
    Lum_pow_max=Lum_pow_max_0;
    
    % Calculating max EQE
    EQE_max_0=0;
    V_EQE_max_0=0;
        for i=1:size(V,1)
           if EQE(i,1)> EQE_max_0 && V(i,1) >= V_on_1 
           EQE_max_0=EQE(i,1); 
           V_EQE_max_0=V(i,1);
           end
        end
    EQE_max=EQE_max_0;
    end
    data_color=zeros(size(V,1)+1,3+9);
    data_color(2:1+size(V,1),1) = V(:,1);
    data_color(2:1+size(V,1),2) = abs(Jden_cm2(:,1));
    data_color(2:1+size(V,1),3) = Luminance(:,1);
    data_color(2:1+size(V,1),4:size(data_color,2)) = data2(:,:);
    data_color=num2cell(data_color);
    data_color{1,1} = 'Voltage (V)';
    data_color{1,2} = 'Current Density (mA/cm2)';
    data_color{1,3} = 'Luminance (cd/m2)';
    data_color(1,4:size(data_color,2)) = {'CIE31_x','CIE31_y','CIE31_X','CIE31_Y','CIE31_Z','CIE60_u','CIE60_v','color correlated temperature (K)','duv'};
    
    data=zeros(size(V,1)+1,10);
    data(2:1+size(V,1),1)=V(:,1); % in V
    data(2:1+size(V,1),2)=J(:,1); % in A
    data(2:1+size(V,1),3)=Jden_cm2(:,1); % in mA/cm2
    data(2:1+size(V,1),4)=abs(Jden_cm2(:,1));
    data(2:1+size(V,1),5)=Lum_int(:,1); % in Cd
    data(2:1+size(V,1),6)=Luminance(:,1); % in Cd/m2
    data(2:1+size(V,1),7)=Lum_eff(:,1); % in Cd/A
    data(2:1+size(V,1),8)=Lum_pow_eff(:,1); % in lm/W
    data(2:1+size(V,1),9)=EQE(:,1); % in percentage
    data(2:1+size(V,1),10)=n_wp(:,1); % in percentage
    data=num2cell(data);
    data{1,1} = 'Voltage (V)';
    data{1,2} = 'Current (A)';
    data{1,3} = 'Current Density (mA/cm2)';
    data{1,4} = 'Current Density (mA/cm2)';
    data{1,5} = 'Luminous Intensity (cd)';
    data{1,6} = 'Luminance (cd/m2)';
    data{1,7} = 'Current Efficiency (cd/A)';
    data{1,8} = 'Luminous Power Efficiency (lm/W)';
    data{1,9} = 'EQE (%)';
    data{1,10} = 'wall plug efficiency (%)';
 
    data_peak = zeros(size(V,1)+1,3+8);
    data_peak(2:1+size(V,1),1) = V(:,1);
    data_peak(2:1+size(V,1),2) = abs(Jden_cm2(:,1));
    data_peak(2:1+size(V,1),3) = Luminance(:,1);
    data_peak(2:1+size(V,1),4:size(data_peak,2)) = data3(:,:);
    data_peak=num2cell(data_peak);
    data_peak{1,1} = 'Voltage (V)';
    data_peak{1,2} = 'Current Density (mA/cm2)';
    data_peak{1,3} = 'Luminance (cd/m2)';
    data_peak(1,4:size(data_peak,2)) = {'Amp', 'Amp Std Err', 'Peak Pos (nm)', 'Peak Pos Std Err', 'FWHM', 'FWHM Std Err', 'R2', 'Traced Peak Pos'};
    
    par0=[dev_area,Vth,Jth,Lum_max,V_Luminance_max,CurEff_max,V_Lum_eff_max,Lum_pow_max,V_Lum_pow_max_0,Jrev,EQE_max,V_EQE_max_0]; % 11 element
    par1=[V_on_100,Jden_100,Lum_100,EQE_100,Lum_eff_100,Lum_pow_eff_100,V_on_1000,Jden_1000,Lum_1000,EQE_1000,Lum_eff_1000,Lum_pow_eff_1000]; % 12 element
    par2=[xy31(1,1),xy31(1,2),xyz31(1,1),xyz31(1,2),xyz31(1,3),uv1960(1,1),uv1960(1,2),cct,duv]; % 9 element
%     
%     data_out=cell(size(data));
%     data_out(1,:)=data(1,:);
%     for a=2:size(data,1)
%         for b=1:size(data,2)
%         data_out(a,b)= cellstr(num2str(data{a,b}));
%         end
%     end
%     clear a b
%     xlswrite([folder '\' file '_data.xlsx'],data_out,1);
%     clear data_out
%     
%     par0_out=cell(2,size(par0,2));
%     par0_out(1,:)={'V1 (V)','Jth1 (mA/cm2)','Max Lum (cd/m2)','V @ Max Lum (V)','Max C.E.(cd/A)',...
%             'V @ Max C.E.(V)','Max L.P.E. (lm/W)','V @ Max L.P.E.(V)','Jrev','Max EQE','V @ Max EQE'};
%         for b=1:size(par0,2)
%             
%             par0_out(2,b)= cellstr(num2str(par0(1,b)));
%         end
%     clear b
%     xlswrite([folder '\' file '_data.xlsx'],par0_out,2);
%     clear par0_out
%     
%     spec_Lum_max_out=cell(size(spec_Lum_max));
%     spec_Lum_max_out(1,:)=spec_Lum_max(1,:);
%     for a=2:size(spec_Lum_max,1)
%         for b=1:size(spec_Lum_max,2)
%         spec_Lum_max_out(a,b)= cellstr(num2str(spec_Lum_max{a,b}));
%         end
%     end
%     clear a b
%     
%     xlswrite([folder '\' file '_data.xlsx'],spec_Lum_max_out,3);
%     clear spec_Lum_max_out
%     par1_out=cell(2,size(par1,2));
%     par1_out(1,:)={'V100 (V)','Lum100 (cd/m2)','EQE100 (%)','CE100 (cd/A)','LPE100 (lm/W)','V1000 (V)','Lum1000(cd/m2)','EQE1000 (%)','CE1000 (cd/A)','LPE1000 (lm/W)'};
% 
%         for b=1:size(par1,2)
%             a=1;
%             par1_out(a,b)= cellstr(num2str(par1(a,b)));
%         end
% 
%     clear a b
%     xlswrite([folder '\' file '_data.xlsx'],par1_out,4);
%     clear par1_out
%     
%     par2_out=cell(2,size(par2,2));
%     par2_out(1,:)={'CIE31_x','CIE31_y','CIE31_X','CIE31_Y','CIE31_Z','CIE60_u','CIE60_v','color correlated temperature (K)','duv'};
%         for b=1:size(par2,2)
%             par2_out(2,b)= cellstr(num2str(par2(1,b)));
%         end
%     clear b
%     xlswrite([folder '\' file '_data.xlsx'],par2_out,5);
%     clear par2_out
%     
%     par3_out=cell(2,size(par3,2));
%     par3_out(1,:)={'Amp', 'Amp Std Err', 'Peak Pos (nm)', 'Peak Pos Std Err (nm)', 'FWHM (nm)', 'FWHM Std Err (nm)', 'R2', 'Traced Peak Pos (nm)'};    
%         for b=1:size(par3,2)
%             par3_out(2,b)= cellstr(num2str(par3(1,b)));
%         end
%     clear b
%     xlswrite([folder '\' file '_data.xlsx'],par3_out,6);
%     clear par3_out
%     
%     data_color_out=cell(size(data_color));
%     data_color_out(1,:)=data_color(1,:);
%     for a=2:size(data_color,1)
%         for b=1:size(data_color,2)
%         data_color_out(a,b)= cellstr(num2str(data_color{a,b}));
%         end
%     end
%     clear a b
%     xlswrite([folder '\' file '_data.xlsx'],data_color_out,7);
%     clear data_color_out
%         
%     data_peak_out=cell(size(data_peak));
%     data_peak_out(1,:)=data_peak(1,:);
%     for a=2:size(data_peak,1)
%         for b=1:size(data_peak,2)
%         data_peak_out(a,b)= cellstr(num2str(data_peak{a,b}));
%         end
%     end
%     clear a b
%     xlswrite([folder '\' file '_data.xlsx'],data_peak_out,8);
%     clear data_peak_out
%     
%     
 end
