function[folder_list]= FolderList(data_folder)
d = dir(data_folder);
isub = [d(:).isdir];
folder_list= {d(isub).name}';
folder_list(ismember(folder_list,{'.','..'})) = [];
end
