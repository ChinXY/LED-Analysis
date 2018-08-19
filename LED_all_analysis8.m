function [all_analysis]=LED_all_analysis8(data_folder,WL_limit,Vth_def)
format short e;
[folder_list]= FolderList(data_folder);

all_analysis{1,1}=0;
all_analysis{1,2}='folder';
all_analysis{1,3}='data_all';
all_analysis{1,4}='par00';
all_analysis{1,5}='spectra';
all_analysis{1,6}='par11';
all_analysis{1,7}='par22';
all_analysis{1,8}='par33';
all_analysis{1,9}='data_color_all';
all_analysis{1,10}='data_peak_all';
    

for i=1:size(folder_list,1)
    folder=[data_folder '\' folder_list{i,1}];
    [par00,par11,par22,par33,data_all,spectra,data_color_all,data_peak_all,num_file]=LED_all_analysis7(folder,WL_limit,Vth_def);
    all_analysis{i+1,1}=i;
    all_analysis{i+1,2}=folder_list{i,1};
    all_analysis{i+1,3}=data_all;
    all_analysis{i+1,4}=par00;
    all_analysis{i+1,5}=spectra;
    all_analysis{i+1,6}=par11;
    all_analysis{i+1,7}=par22;
    all_analysis{i+1,8}=par33;
    all_analysis{i+1,9}=data_color_all;
    all_analysis{i+1,10}=data_peak_all;
    
    if num_file == 0
    else
    for a=2:size(data_all,1)
        for b=1:size(data_all,2)
        data_all(a,b)= cellstr(num2str(data_all{a,b}));
        end
    end
    clear a b
    
    for a=2:size(par00,1)
        for b=1:size(par00,2)-1
        par00(a,b)= cellstr(num2str(par00{a,b}));
        end
    end
    clear a b
    
    for a=2:size(par11,1)
        for b=1:size(par11,2)
        par11(a,b)= cellstr(num2str(par11{a,b}));
        end
    end
    clear a b
    
    for a=2:size(par22,1)
        for b=1:size(par22,2)
        par22(a,b)= cellstr(num2str(par22{a,b}));
        end
    end
    clear a b
    
    for a=2:size(par33,1)
        for b=1:size(par33,2)
        par33(a,b)= cellstr(num2str(par33{a,b}));
        end
    end
    clear a b
    
    for a=2:size(spectra,1)
        for b=1:size(spectra,2)
        spectra(a,b)= cellstr(num2str(spectra{a,b}));
        end
    end
    clear a b
    
    for a=2:size(data_color_all,1)
        for b=1:size(data_color_all,2)
        data_color_all(a,b)= cellstr(num2str(data_color_all{a,b}));
        end
    end
    clear a b
    
    for a=2:size(data_peak_all,1)
        for b=1:size(data_peak_all,2)
        data_peak_all(a,b)= cellstr(num2str(data_peak_all{a,b}));
        end
    end
    clear a b
    
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],data_all,1);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],par00(:,1:size(par00,2)-1),2);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],spectra,3);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],par11,4);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],par22,5);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],par33,6);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],data_color_all,7);
    xlwrite([data_folder '\' folder_list{i,1} '_data_all.xlsx'],data_peak_all,8);
    


    end
end

save([data_folder '\' ' all analysis new.mat'], 'all_analysis','data_folder','WL_limit','Vth_def');
end
