n_repeats = 100;  % number of times to run fit with new v0
days = [0, 3, 6, 9];  % timepoints to analyze 
data_folder = 'input_folder';

output_table = table(...
    'Size', [6 * n_repeats, 19], ...
    'VariableTypes', {...
        'string', 'string', ...
        'double', 'double', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'double', 'logical', 'logical'}, ...
    'VariableNames', {...
        'cell_line', 'conc', ... 
        'g_l0', 'd_l0', 'r_lh0', 'r_hl0', 'g_h0', 'd_h0', ...
        'g_l', 'd_l', 'r_lh', 'r_hl', 'g_h', 'd_h', ...
        'i', 'err_null', 'err', 'ineq', 'trns'});
    

%% Find parameter values for n random intializations + starting error
rowi = 1;
for cell_line = ["WT", "KO"]
    for conc = ["0", "1", "4"]
        disp("Running Tet2" + cell_line + " " + conc + "uM, " + datestr(now, 'HH:MM:SS') + "...");
        trns_mat = fit_fn(cell_line, conc, days, data_folder, n_repeats);
        output_table(rowi:(rowi+n_repeats-1), :) = trns_mat;
        rowi = rowi + n_repeats;
    end
end
writetable(output_table, "output.txt", 'Delimiter', '\t');

disp("Done.");







