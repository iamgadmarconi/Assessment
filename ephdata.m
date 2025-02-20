classdef ephdata
    methods (Static)
        function tbl = read_data(path)
            % Read CSV file and return table
            % Usage: data = data.read_data('path/to/file.csv');
            tbl = readtable(path);
            fprintf('Loaded data with %d rows and %d columns\n', size(tbl, 1), size(tbl, 2));
        end
        
        function filtered = filter_by_eccentricity(tbl, max_e)
            % Filter data by maximum eccentricity
            % Usage: filtered = data.filter_by_eccentricity(data_table, 0.5);
            filtered = tbl(tbl.e < max_e, :);
        end
        
        function [cleaned_names, tbl] = clean_object_names(tbl)
            % Clean object names from full_name column
            raw_names = tbl.full_name;
            
            % First remove parentheses, then trim whitespace
            cleaned_names = cellfun(@(x) ...
                strtrim(...                          % Remove leading/trailing spaces
                    strrep(strrep(x, '(', ''), ')', '')...  % Remove both parentheses
                ), ...
                raw_names, 'UniformOutput', false);
            
            tbl.object_name = cleaned_names;
        end
        
        function subset = select_columns(tbl, columns)
            % Select specific columns from table
            % Usage: subset = data.select_columns(data_table, {'a', 'e', 'object_name'});
            subset = tbl(:, columns);
        end
        
        function filtered = filter_data(tbl, condition)
            % Generic filter using logical array
            % Usage: filtered = data.filter_data(data_table, data_table.H < 20);
            filtered = tbl(condition, :);
        end
        
        function stats = get_stats(tbl, columns)
            % Get basic statistics for specified columns
            % Usage: stats = data.get_stats(data_table, {'a', 'e', 'i'});
            stats = struct();
            for i = 1:length(columns)
                col = columns{i};
                stats.(col) = struct(...
                    'mean', mean(tbl.(col)), ...
                    'std', std(tbl.(col)), ...
                    'min', min(tbl.(col)), ...
                    'max', max(tbl.(col)) ...
                );
            end
        end
        
        function unique_objects = get_unique_objects(tbl)
            % Get unique object names
            [~, idx] = unique(tbl.object_name);
            unique_objects = tbl(idx, :);
        end

        function id = getAsteroidId(tbl, name)
            % getAsteroidId returns the numeric ID for the given asteroid name.
            idx = find(strcmpi(tbl.full_name, name), 1);
            
            if isempty(idx)
                error('Asteroid "%s" not found in the database.', name);
            end
            
            % Map the row index to the asteroid ID.
            id = idx + 11;
        end
    end
end 