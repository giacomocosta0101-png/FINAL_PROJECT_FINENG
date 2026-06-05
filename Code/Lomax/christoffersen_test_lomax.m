function out = christoffersen_test_lomax(exceptions)
% CHRISTOFFERSEN_TEST_LOMAX  Run POF, independence and conditional-coverage tests.
%
%   out = christoffersen_test_lomax(exceptions)
%
%   INPUT
%     exceptions : VaR exceptions for one model (column 1 = 95%, column 2 =
%                  99%). Accepts the Fixed-mode M x 2 layout or the
%                  Rolling-window 1 x 2 x M array.
%
%   OUTPUT
%     out        : struct; out.level(L) holds the counts, LR
%                  statistics and p-values for level L.

% Rolling-window exceptions come as 1 x 2 x M; bring them back to M x 2.
if ndims(exceptions) == 3 && size(exceptions,1) == 1 && size(exceptions,2) == 2
    exceptions = squeeze(exceptions).';
end

%% Tests
levels    = [0.05, 0.01];
lvl_names = {'95%','99%'};
out       = struct();

for i = 1:2
    out.level(i) = single_test(double(exceptions(:,i)), ...
                               levels(i), lvl_names{i});
end

print_results(out);

end


function s = single_test(I, p, label)
% Run sthe three VaR backtests on a single hit series.
%   I     : N x 1 vector of 0/1 VaR exceptions (the "hit" sequence)
%   p     : nominal tail probability (e.g. 0.05 or 0.01)
%   label : level name for reporting (e.g. '95%')

    N = numel(I);
    x = sum(I);                         % observed number of exceptions

    % First-order Markov transition counts (state i at t-1 -> j at t)
    % We can model the consecutive set as: 2*I(t-1)+I(t): 0->00, 1->01, 2->10, 3->11

    tr  = 2*I(1:end-1) + I(2:end);
    n00 = sum(tr==0);   n01 = sum(tr==1);
    n10 = sum(tr==2);   n11 = sum(tr==3);

    pi01 = n01/(n00+n01);                   % P(hit today | no hit yesterday)
    pi11 = n11/(n10+n11);                   % P(hit today | hit yesterday)
    pi   = (n01+n11)/(n00+n01+n10+n11);     % overall hit rate

    % Christoffersen Independence: iid (H0) vs first-order Markov (H1)
    logL_iid    = (n00+n10)*log(1-pi)  + (n01+n11)*log(pi);
    logL_markov = n00*log(1-pi01) + n01*log(pi01) + ...
                  n10*log(1-pi11) + n11*log(max(pi11,eps)); % since log(0) = NaN, and here is very likely to occur, we hardcoded here
    LR_ind   = -2*(logL_iid - logL_markov);
    pval_ind = 1 - chi2cdf(LR_ind, 1);

    % Kupiec Test: unconditional coverage (hit rate == p)
    pi_hat = x/N;

    LR_kupiec = -2*((N-x)*log(1-p)      + x*log(p) ...
                    - (N-x)*log(1-pi_hat) - x*log(pi_hat) );

    pval_kupiec = 1 - chi2cdf(LR_kupiec, 1);

    % Christoffersen conditional coverage:
    LR_Christoffersen  = LR_kupiec + LR_ind;
    pval_Christoffersen = 1 - chi2cdf(LR_Christoffersen, 2);

    s = struct( ...
        'label',    label,   'N',   N,   'x', x,   'expected', p*N, ...
        'n00', n00, 'n01', n01, 'n10', n10, 'n11', n11, ...
        'pi01', pi01, 'pi11', pi11, 'pi', pi, ...
        'LR_kupiec', LR_kupiec, 'pval_kupiec', pval_kupiec, ...
        'LR_independence', LR_ind, 'pval_independence', pval_ind, ...
        'LR_christoffersen',  LR_Christoffersen,  'pval_christoffersen',  pval_Christoffersen);
end


function print_results(out)
    bar = repmat('=',1,98);
    fprintf('\n%s\n  Christoffersen / Kupiec  —  Kupiec + Independence + Conditional Coverage\n%s\n', bar, bar);
    fprintf('%-6s %-6s %4s %5s | %4s %4s %4s %4s | %8s %7s | %8s %7s | %8s %7s\n', ...
        'Model','Lev','x','exp','n00','n01','n10','n11', ...
        'LR_kup','p_kup','LR_ind','p_ind','LR_christ','p_christ');
    fprintf('%s\n', repmat('-',1,98));
    for m = 1:numel(out)
        for L = 1:2
            s = out(m).level(L);
            fprintf('%-6s %-6s %4d %5.1f | %4d %4d %4d %4d | %8.3f %7.4f | %8.3f %7.4f | %8.3f %7.4f\n', ...
                sprintf('M%d',m), s.label, s.x, s.expected, ...
                s.n00, s.n01, s.n10, s.n11, ...
                s.LR_kupiec, s.pval_kupiec, ...
                s.LR_independence, s.pval_independence, ...
                s.LR_christoffersen,  s.pval_christoffersen);
        end
    end
    fprintf('%s\n\n', bar);
end
