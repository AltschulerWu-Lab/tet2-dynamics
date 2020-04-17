n_cd38lo = 0;
n_cd38hi = 5;
time = 6;
n_reps = 10;
n_cell_threshold = 50;

disp("Run for " + n_cd38lo + " CD38lo & " + n_cd38hi + " CD38hi cells for " + time + " time steps.");
for cell_line = ["KG1", "KG1T"]
    for conc = "1"
        disp("Gillespie for " + cell_line + " at " + conc + "uM AraC");
        disp("% colonies grown: ");
        for j = 1:n_reps
            grown_out = zeros(1, 100);
            for i = 1:100
                [tplot, y_S, y_D, y_n] = gillespie3x3( cell_line, conc, n_cd38lo, n_cd38hi, time);
%                 grown_out(i) = y_n(end);
                if y_n(end) > n_cell_threshold
                    grown_out(i) = 1;
                end
            end
            fprintf(num2str(mean(grown_out)) + ", ");
        end
    end
    fprintf("\n\n")
end