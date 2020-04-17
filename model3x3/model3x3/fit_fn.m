function output_table = fit_fn( cell_line, conc, days, data_folder, n_repeats )

    start_days = days(1:end-1);
    dat = cell(length(start_days) * 2, 4);
    k = 1;
    for day=start_days
        for j={'lo', 'hi'}
            datlo = readtable(strcat(data_folder, '/KG1_Tet2', cell_line, '_', j{1}, num2str(day), '_', conc, 'uM.txt'));
            dat{k, 1} = table2array(datlo(:,1));
            dat{k, 2} = table2array(datlo(:,2));
            dat{k, 3} = table2array(datlo(:,5));
            dat{k, 4} = table2array(datlo(:,6));
            k = k + 1;
        end
    end
    
    options = optimset('MaxIter', 100000, 'MaxFunEvals', 100000, 'TolX', 1e-5);
    output_table = table(...
        'Size', [n_repeats, 19], ...
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
    
    rowi = 1;
    for iter = 1:n_repeats
        guess = rand(1,6);
        bestfit = fminsearch(@f2, guess, options);
        
        g_l0 = guess(1);
        d_l0 = guess(2);
        r_lh0 = guess(3);
        r_hl0 = guess(4);
        g_h0 = guess(5);
        d_h0 = guess(6);
        g_l = bestfit(1);
        d_l = bestfit(2);
        r_lh = bestfit(3);
        r_hl = bestfit(4);
        g_h = bestfit(5);
        d_h = bestfit(6);
        net_l = g_l - d_l;
        net_h = g_h - d_h;
        
        err_null = f2(guess);
        err = f2(bestfit);
        ineq = (net_l) > (net_h);
        trns = r_lh > r_hl;

        output_table(rowi,:) = table(...
            cell_line, conc, ...
            g_l0, d_l0, r_lh0, r_hl0, g_h0, d_h0, ...
            g_l, d_l, r_lh, r_hl, g_h, d_h, ...
            iter, err_null, err, ineq, trns); 
        rowi = rowi + 1;
    end
    
    function outval = f2(v1)
        g_lf = v1(1);
        d_lf = v1(2);
        net_lf = g_lf - d_lf;
        r_lhf = v1(3);
        r_hlf = v1(4);
        g_hf = v1(5);
        d_hf = v1(6);
        net_hf = g_hf - d_hf;
        
        fit_err = zeros(1, length(dat));
        for i=1: size(dat,1)
            x11 = cell2mat(dat(i, 1));
            x12 = cell2mat(dat(i, 2));
            x1 = [(x11*(1 + net_lf - r_lhf)) + (r_hlf * x12) (x12*(1 + net_hf - r_hlf)) + (r_lhf * x11)];
            y1 = cell2mat(dat(i, 3:4));
            fit_err(1, i) = sum(abs(x1 - y1), 'all');
        end
        
        outval = sum(fit_err);
        if r_lhf < 0 || r_hlf < 0 || g_lf < 0 || d_lf < 0 || g_hf < 0 || d_hf < 0
            outval = outval + 1e5;
        end
    end
end
    
