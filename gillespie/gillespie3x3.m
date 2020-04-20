function [tplot, splot, dplot, nplot] = gillespie3x3 ( cell_line, conc, s_0, d_0, t_end )
    %------------------------
    % equation         | rate
    %------------------|-----
    % S -> X           | dS
    % S -> 2*S         | gS
    % S -> D           | rSD
    % D -> X           | dD
    % D -> S           | rDS
    % D -> 2*D         | gD
    %------------------------
    
    % starting population
    n_S = s_0;
    n_D = d_0;
    n_X = 0.0;
    
    % time
    t = 0.0;
    
    % rates KG1WT 0uM
    if cell_line == "KG1" && conc == "0"
        d_S = 0.03;
        g_S = 1.32;
        r_SD = 1.38;
        d_D = 0.36;
        r_DS = 0.022;
        g_D = 0.37;
    elseif cell_line == "KG1" && conc == "1"
        % rates KG1WT 1uM
        d_S = 0.75;
        g_S = 0.85;
        r_SD = 0.59;
        d_D = 0.64;
        r_DS = 0.00;
        g_D = 0.38;
    elseif cell_line == "KG1T" && conc == "0"
        % rates KG1KO 0uM
        d_S = 0.04;
        g_S = 1.12;
        r_SD = 0.08;
        d_D = 0.33;
        r_DS = 0.36;
        g_D = 0.76;
    elseif cell_line == "KG1T" && conc == "1"
        % rates KG1KO 1uM
        d_S = 0.77;
        g_S = 0.80;
        r_SD = 0.13;
        d_D = 0.65;
        r_DS = 0.05;
        g_D = 0.38;
    end
    
    % In general, a vector of length M of stochastic reaction constants.
    splot = n_S;
    dplot = n_D;
    nplot = n_S + n_D;
    tplot = zeros(1);
    n = 0;

    while t < t_end
        w1 = n_S * d_S; 
        w2 = n_S * g_S; 
        w3 = n_S * r_SD;
        w4 = n_D * d_D;
        w5 = n_D * r_DS;
        w6 = n_D * g_D;
        a = [w1 w2 w3 w4 w5 w6];
        a0 = sum(a);
        if a0 == 0
            break
        end
        
        dt = -log(rand(1,1)) / a0;
        t = t + dt;
        
        % Note 0<=mu<=M. This decides which reaction occurs. 
        mu = randp(a/a0, 1,1);
        
        % Adjust population levels based on reaction formula(s).
        if mu == 1
            n_S = n_S - 1;
            n_X = n_X + 1;
        elseif mu == 2
            n_S = n_S + 1;
        elseif mu == 3
            n_S = n_S - 1;
            n_D = n_D + 1;
        elseif mu == 4
            n_D = n_D - 1;
            n_X = n_X + 1;
        elseif mu == 5
            n_D = n_D - 1;
            n_S = n_S + 1;
        elseif mu == 6
            n_D = n_D + 1;
        end
        
        if n_S < 0
            n_S = 0;
        end
        if n_D < 0
            n_D = 0;
        end
        if n_S + n_D > 500
            break
        end
        
        n = n+1;
        splot(n+1) = n_S;
        dplot(n+1) = n_D;
        nplot(n+1) = n_S + n_D;
        tplot(n+1) = t;  
    end
end