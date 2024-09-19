function process_all_tiff_folders(root_path, save_root_path)
    % PROCESS_ALL_TIFF_FOLDERS Processes all subfolders containing TIFF files under the root path.
    % It creates TIFF stacks for each subfolder using create_tiff_stack function.
    %
    % Parameters:
    % root_path - String specifying the root directory to start searching for TIFF folders.
    % save_root_path - String specifying the root directory to save the TIFF stack files.
    %
    % Example:
    % process_all_tiff_folders('E:\1_Data\Luorong\20240709_optopatch', 'E:\1_Data\Luorong\processed\');
    
    % Get all subfolders
    subfolders = get_subfolders_with_tiffs(root_path);
    
    % Process each subfolder
    for i = 1:length(subfolders)
        fprintf('Processing folder: %s\n', subfolders{i});
        % Create corresponding save path
        relative_path = strrep(subfolders{i}, root_path, '');
        save_path = fullfile(save_root_path, relative_path);
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        % Create TIFF stack
        create_tiff_stack(subfolders{i}, save_path);
        
        % Copy non-TIFF files to the new folder
        copy_non_tiff_files(subfolders{i}, save_path);
    end
end

function subfolders = get_subfolders_with_tiffs(root_path)
    % GET_SUBFOLDERS_WITH_TIFFS Recursively finds all subfolders containing TIFF files.
    %
    % Parameters:
    % root_path - String specifying the root directory to start searching.
    %
    % Returns:
    % subfolders - Cell array of strings, each string is a path to a subfolder containing TIFF files.

    subfolders = {};
    folder_list = dir(root_path);
    
    for i = 1:length(folder_list)
        if folder_list(i).isdir && ~strcmp(folder_list(i).name, '.') && ~strcmp(folder_list(i).name, '..')
            current_folder = fullfile(root_path, folder_list(i).name);
            tif_files = dir(fullfile(current_folder, '*.tif'));
            if ~isempty(tif_files)
                subfolders{end+1} = current_folder;
            end
            % Recursively search subdirectories
            subfolders = [subfolders, get_subfolders_with_tiffs(current_folder)];
        end
    end
end

function copy_non_tiff_files(source_folder, destination_folder)
    % COPY_NON_TIFF_FILES Copies non-TIFF files from source_folder to destination_folder.
    %
    % Parameters:
    % source_folder - String specifying the folder to copy files from.
    % destination_folder - String specifying the folder to copy files to.

    % Get list of all files in the source folder
    all_files = dir(source_folder);
    
    for i = 1:length(all_files)
        if ~all_files(i).isdir && ~endsWith(all_files(i).name, '.tif')
            source_file = fullfile(source_folder, all_files(i).name);
            destination_file = fullfile(destination_folder, all_files(i).name);
            copyfile(source_file, destination_file);
        end
    end
end
