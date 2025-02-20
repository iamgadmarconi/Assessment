clc; clear; close all;

if ~exist('lambertsolver', 'class')
    error('lambertsolver class not found. Verify file location and MATLAB path.');
end

dep_date = date2mjd2000([2033, 01, 01, 0, 0, 0])
arr_date = date2mjd2000([2037, 12, 31, 0, 0, 0])

