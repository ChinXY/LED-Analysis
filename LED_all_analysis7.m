function [par00,par11,par22,par33,data_all,spectra,data_color_all,data_peak_all,num_file]=LED_all_analysis7(folder,WL_limit,Vth_def)

[file_array,num_file,date_array]=file_list2(folder);

disp(' ');
disp(' ');
disp(['For folder ' folder]);
disp(['Number of files = ' num2str(num_file)]);
[par00,par11,par22,par33,data_all,spectra,data_color_all,data_peak_all]=LED_all_analysis6(folder,file_array,date_array,WL_limit,Vth_def);
end