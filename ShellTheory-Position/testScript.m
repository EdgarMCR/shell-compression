close all;
clear all;
H=0;
Np=1000;
counter=0;
errors=0;
for H=1:-0.02:0.5
        try
            fprintf('H=%d \n\n ', H);
            [X, P] = initialGuess(Np, H);
            savepath=sprintf('Solution_%d_H=%.1f_p=%.1f.txt', counter, H, P);
            problem(Np,  H, P, X, savepath)
        catch
            errors=errors+1;
        end
        counter=counter+1;
end

fprintf('\n\n Total runs = %d with %d errors, sucessful = %d.\n', counter, errors, counter-errors);