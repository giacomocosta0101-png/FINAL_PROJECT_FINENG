function out = christoffersen_test(exceptions)
% CHRISTOFFERSEN_TEST  Test di indipendenza + Kupiec POF + Conditional
% Coverage per ogni modello e ogni livello di confidenza.
%
%   out = christoffersen_test(exceptions)
%
%   exceptions: cell {1xM} di matrici Nx2 (col1 = 95%, col2 = 99%)
%   out:        struct con conteggi, statistiche e p-value

    levels    = [0.05, 0.01];
    lvl_names = {'95%','99%'};
    M = numel(exceptions);
    out = struct();

    for m = 1:M
        for L = 1:2
            I = double(exceptions{m}(:,L));   % serie binaria di eccezioni
            N = numel(I);
            p = levels(L);

            % ---- conteggio transizioni (Markov del primo ordine) ----
            tr  = I(1:end-1)*2 + I(2:end);   % 0=00, 1=01, 2=10, 3=11
            n00 = sum(tr==0);
            n01 = sum(tr==1);
            n10 = sum(tr==2);
            n11 = sum(tr==3);

            % ---- stime di probabilità ----
            pi01 = safe_div(n01, n00+n01);
            pi11 = safe_div(n11, n10+n11);
            pi   = safe_div(n01+n11, n00+n01+n10+n11);

            % ---- Independence (Christoffersen) ----
            logL0 = xlogy(n00+n10, 1-pi)   + xlogy(n01+n11, pi);
            logL1 = xlogy(n00, 1-pi01) + xlogy(n01, pi01) + ...
                    xlogy(n10, 1-pi11) + xlogy(n11, pi11);
            LR_ind   = -2*(logL0 - logL1);
            pval_ind = 1 - chi2cdf(LR_ind, 1);

            % ---- Kupiec POF (unconditional coverage) ----
            x      = sum(I);
            pi_hat = x/N;
            if x==0 || x==N
                LR_POF = 0;
            else
                LR_POF = -2*( xlogy(N-x,1-p) + xlogy(x,p) ...
                            - xlogy(N-x,1-pi_hat) - xlogy(x,pi_hat) );
            end
            pval_POF = 1 - chi2cdf(LR_POF, 1);

            % ---- Conditional Coverage (joint) ----
            LR_cc   = LR_POF + LR_ind;
            pval_cc = 1 - chi2cdf(LR_cc, 2);

            % ---- salvo ----
            out(m).level(L).label    = lvl_names{L};
            out(m).level(L).N        = N;
            out(m).level(L).x        = x;
            out(m).level(L).expected = p*N;
            out(m).level(L).n00 = n00;  out(m).level(L).n01 = n01;
            out(m).level(L).n10 = n10;  out(m).level(L).n11 = n11;
            out(m).level(L).pi01 = pi01;  out(m).level(L).pi11 = pi11;
            out(m).level(L).pi   = pi;
            out(m).level(L).LR_POF  = LR_POF;   out(m).level(L).pval_POF = pval_POF;
            out(m).level(L).LR_ind  = LR_ind;   out(m).level(L).pval_ind = pval_ind;
            out(m).level(L).LR_cc   = LR_cc;    out(m).level(L).pval_cc  = pval_cc;
        end
    end

    print_results(out);
end

% ============== helper ==============
function y = xlogy(x, p)
    if x==0,  y = 0;  else, y = x*log(p);  end
end

function v = safe_div(a, b)
    if b==0,  v = 0;  else, v = a/b;  end
end

function print_results(out)
    bar = repmat('=',1,98);
    fprintf('\n%s\n  Christoffersen / Kupiec  —  POF + Independence + Conditional Coverage\n%s\n', bar, bar);
    fprintf('%-6s %-6s %4s %5s | %4s %4s %4s %4s | %8s %7s | %8s %7s | %8s %7s\n', ...
        'Model','Lev','x','exp','n00','n01','n10','n11', ...
        'LR_POF','p_POF','LR_ind','p_ind','LR_cc','p_cc');
    fprintf('%s\n', repmat('-',1,98));
    for m = 1:numel(out)
        for L = 1:2
            s = out(m).level(L);
            fprintf('%-6s %-6s %4d %5.1f | %4d %4d %4d %4d | %8.3f %7.4f | %8.3f %7.4f | %8.3f %7.4f\n', ...
                sprintf('M%d',m), s.label, s.x, s.expected, ...
                s.n00, s.n01, s.n10, s.n11, ...
                s.LR_POF, s.pval_POF, ...
                s.LR_ind, s.pval_ind, ...
                s.LR_cc,  s.pval_cc);
        end
    end
    fprintf('%s\n\n', bar);
end