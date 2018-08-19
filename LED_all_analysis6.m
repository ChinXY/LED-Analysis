function [par00,par11,par22,par33,data_all,spectra,data_color_all,data_peak_all]=LED_all_analysis6(folder,file_array,date_array,WL_limit,Vth_def)
    format short e;
    size_x=0;
    if size(file_array,2)==0
        par00=0;
        par11=0;
        par22=0;
        par33=0;
        data_all=0;
        spectra=0;
        data_color_all=0;
        data_peak_all=0;

    else
        for i=1:size(file_array,2)
        file=[num2str(date_array(1,i)) '_' file_array{1,i} 'data.ciemo'];
        filename=[folder '\' file];
        FileName_data=filename;
        fid = fopen(FileName_data, 'r');
        JVLdata = textscan(fid,'%f %f %f %f %f','HeaderLines',1);
        V = JVLdata{1,1}; % un V
        fclose(fid);
            if size(V,1)>size_x
               size_x=size(V,1);
            else
            end
        end
        clear V
        par00=cell(1+size(file_array,2),2+11+1+2);
        par11=cell(1+size(file_array,2),2+12+1);
        par22=cell(1+size(file_array,2),2+9+1);
        par33=cell(size(file_array,2),2+8+1);
        data_all=cell(size_x+3, 10*size(file_array,2));
        data_color_all=cell(size_x+3,12*size(file_array,2));
        data_peak_all=cell(size_x+3,11*size(file_array,2));
        
        for i=1:size(file_array,2)
        file=[num2str(date_array(1,i)) '_' file_array{1,i} 'data.ciemo'];
        file_spec=[num2str(date_array(1,i)) '_' file_array{1,i} 'spectrum_raw.ciemo'];
        fileM=[num2str(date_array(1,i)) '_' file_array{1,i} 'meta.ciemo'];
        [par0,spec_Lum_max,data,par1,par2,par3,data_color,all_spectra_cell,data_peak]=LED_analysis8(folder,file,file_spec,fileM,WL_limit,Vth_def);
        
        
        par00{1+i,1}                    = i;
        par00{1+i,2}                    = date_array(1,i);
        par00{1+i,3}                    = str2double(file_array{1,i});
        par00(1+i,4:size(par00,2)-1)    = num2cell(par0);
        par00{1+i,size(par00,2)}        = all_spectra_cell;
        
        par11{1+i,1}                    = i;
        par11{1+i,2}                    = date_array(1,i);
        par11{1+i,3}                    = str2double(file_array{1,i});
        par11(1+i,4:size(par11,2))      = num2cell(par1);
        
        par22{1+i,1}                    = i;
        par22{1+i,2}                    = date_array(1,i);
        par22{1+i,3}                    = str2double(file_array{1,i});
        par22(1+i,4:size(par22,2))      = num2cell(par2);
        
        par33{1+i,1}                    = i;
        par33{1+i,2}                    = date_array(1,i);
        par33{1+i,3}                    = str2double(file_array{1,i});
        par33(1+i,4:size(par33,2))       = num2cell(par3);
        
        a0(1,1:10)=i;
        a1(1,1:10)=str2double(file_array{1,i});
        data_all(1,10*(i-1)+1:10*(i-1)+10)=num2cell(a0);
        data_all(2,10*(i-1)+1:10*(i-1)+10)=num2cell(a1);
        data_all(3:2+size(data,1),10*(i-1)+1:10*(i-1)+10)=data(:,:);
        clear a0 a1
        
        a0(1,1:12)=i;
        a1(1,1:12)=str2double(file_array{1,i});
        data_color_all(1,12*(i-1)+1:12*(i-1)+12)=num2cell(a0);
        data_color_all(2,12*(i-1)+1:12*(i-1)+12)=num2cell(a1);
        data_color_all(3:2+size(data_color,1),12*(i-1)+1:12*(i-1)+12)=data_color(:,:);
        clear a0 a1
        
        a0(1,1:4) = i;
        a1(1,1:4) = str2double(file_array{1,i});
        spectra(1,4*(i-1)+1:4*i) = num2cell(a0);
        spectra(2,4*(i-1)+1:4*i) = num2cell(a1);
        spectra(3:2+size(spec_Lum_max),4*(i-1)+1:4*i) = spec_Lum_max(:,:);
        clear a0 a1
        
        a0(1,1:11) = i;
        a1(1,1:11) = str2double(file_array{1,i});
        data_peak_all(1,11*(i-1)+1:11*(i-1)+11)=num2cell(a0);
        data_peak_all(2,11*(i-1)+1:11*(i-1)+11)=num2cell(a1);
        data_peak_all(3:2+size(data_peak,1),11*(i-1)+1:11*(i-1)+11)=data_peak(:,:);
        clear a0 a1
        end
        
        A=data_all(3,:);
        B=data_all(1:2,:);
        data_all(1,:)=A;
        data_all(2:3,:)=B;
        clear A B
        
        A=data_color_all(3,:);
        B=data_color_all(1:2,:);
        data_color_all(1,:)=A;
        data_color_all(2:3,:)=B;
        clear A B
        
        A=spectra(3,:);
        B=spectra(1:2,:);
        spectra(1,:)=A;
        spectra(2:3,:)=B;
        clear A B
        
        A=data_peak_all(3,:);
        B=data_peak_all(1:2,:);
        data_peak_all(1,:)=A;
        data_peak_all(2:3,:)=B;
        clear A B

        par00{1,1} = 'number';
        par00{1,2} = 'date';
        par00{1,3} = 'time';
        par00(1,4:size(par00,2)) = {'DevArea(m2)','V1 (V)','Jth1 (mA/cm2)','Max Lum (cd/m2)','V @ Max Lum (V)','Max C.E.(cd/A)',...
            'V @ Max C.E.(V)','Max L.P.E. (lm/W)','V @ Max L.P.E.(V)','Jrev','Max EQE','V @ Max EQE','all spectra'};
        
        par11{1,1} = 'number';
        par11{1,2} = 'date';
        par11{1,3} = 'time';
        par11(1,4:size(par11,2)) = {'V100 (V)','Jden100 (mA/cm2)','Lum100 (cd/m2)','EQE100 (%)','CE100 (cd/A)','LPE100 (lm/W)','V1000 (V)','Jden100(mA/cm2)','Lum1000(cd/m2)','EQE1000 (%)','CE1000 (cd/A)','LPE1000 (lm/W)'};
        
        par22{1,1} = 'number';
        par22{1,2} = 'date';
        par22{1,3} = 'time';
        par22(1,4:size(par22,2)) = {'CIE31_x','CIE31_y','CIE31_X','CIE31_Y','CIE31_Z','CIE60_u','CIE60_v','color correlated temperature (K)','duv'};
        
        par33{1,1} = 'number';
        par33{1,2} = 'date';
        par33{1,3} = 'time';
        par33(1,4:size(par33,2)) = {'Amp', 'Amp Std Err', 'Peak Pos (nm)', 'Peak Pos Std Err (nm)', 'FWHM (nm)', 'FWHM Std Err (nm)', 'R2', 'Traced Peak Pos (nm)'};
        
    end 
    
end
