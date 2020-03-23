function term = fit_fn( cell_line )

    dat1 = readtable('data/' + cell_line + '_lo0.txt');
    dat2 = readtable('data/' + cell_line + '_hi0.txt');
    dat3 = readtable('data/' + cell_line + '_lo3.txt');
    dat4 = readtable('data/' + cell_line + '_hi3.txt');
    
    x1 = table2array(dat1(:,{'lo', 'hi'}));
    x2 = table2array(dat2(:,{'lo', 'hi'}));
    x3 = table2array(dat3(:,{'lo', 'hi'}));
    x4 = table2array(dat4(:,{'lo', 'hi'}));
    
    y1 = table2array(dat1(:,{'lo_1', 'hi_1'}));
    y2 = table2array(dat2(:,{'lo_1', 'hi_1'}));
    y3 = table2array(dat3(:,{'lo_1', 'hi_1'}));
    y4 = table2array(dat4(:,{'lo_1', 'hi_1'}));
    
    options = optimset('MaxIter', 100000, 'MaxFunEvals', 100000, 'TolX', 1e-7, 'PlotFcns', @optimplotfval); 
    guess_v1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    bestfit = fminsearch(@f2, guess_v1, options);
    term = [bestfit, f2(bestfit)];

    function outval = f2(v1)
        a = v1(1);
        b = v1(2);
        c = v1(3);
        d = v1(4);
        e = v1(5);
        f = v1(6);
        
        x1_t1 = [(x1(:,1)*(1 + b - c - a)) + (e * x1(:,2)) (x1(:,2)*(1 + f - e - d)) + (c * x1(:,1))];
        x2_t1 = [(x2(:,1)*(1 + b - c - a)) + (e * x2(:,2)) (x2(:,2)*(1 + f - e - d)) + (c * x2(:,1))];
        x3_t1 = [(x3(:,1)*(1 + b - c - a)) + (e * x3(:,2)) (x3(:,2)*(1 + f - e - d)) + (c * x3(:,1))];
        x4_t1 = [(x4(:,1)*(1 + b - c - a)) + (e * x4(:,2)) (x4(:,2)*(1 + f - e - d)) + (c * x4(:,1))];
        
        outval = sum(abs(x1_t1 - y1), 'all') + ...
            sum(abs(x2_t1 - y2), 'all') + ...
            sum(abs(x3_t1 - y3), 'all') + ...
            sum(abs(x4_t1 - y4), 'all');
        
        if a < 0 || b < 0 || c < 0 || d < 0 || e < 0 || f < 0
            outval = outval + 1e7;
        end
    end
end
