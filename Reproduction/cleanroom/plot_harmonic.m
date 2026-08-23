function figs = plot_harmonic(R, Ddg, Dda, cfg, variant)
% PLOT_HARMONIC  Figures for the per-vertex harmonic model.
%
%   figs = plot_harmonic(R, Ddg, Dda, cfg, variant)
%
% Writes three figures to cfg.figDir (png + pdf, as elsewhere in the clean room):
%   harmonic_decomposition_<variant>   the A/B/C amplitudes vs polar angle, with the
%                                      model's three curves. The diagnostic that shows
%                                      whether the functional form actually holds.
%   harmonic_coefficients_<variant>    Fit A vs Fit B coefficients with bootstrap CIs.
%   harmonic_crossprediction_<variant> per-wedge tuning, observed vs cross-predicted.

    if ~isfolder(cfg.figDir), mkdir(cfg.figDir); end
    figs = struct();
    figs.decomposition   = fig_decomposition(R, Ddg, Dda, cfg, variant);
    figs.coefficients    = fig_coefficients(R, cfg, variant);
    figs.crossprediction = fig_crossprediction(R, Ddg, Dda, cfg, variant);
end

% ========================================================================
function out = fig_decomposition(R, Ddg, Dda, cfg, variant)
% Four 45-deg-spaced orientations give each vertex exactly three degrees of freedom.
% In the frame where the sampling angles are vertex-independent (absolute for dg,
% radial-relative for da) the model reduces to three scalar regressions on thetaV.
    fig = figure('Color','w','Position',[100 100 1250 620],'Visible','off');
    tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Orientation asymmetries vs pRF polar angle (%s betas)', variant), ...
          'FontWeight','bold');
    subtitle(tl, ['Every panel is a pro-minus-con response difference; offsets and amplitudes ' ...
                  'are the reported 2b values.   Blue departing from red is the context effect.'], ...
             'FontSize',9);

    tvD = linspace(0,360,721).';
    exps = {'dg','Cartesian gratings'; 'da','Polar gratings'};
    for ei = 1:2
        en = exps{ei,1};
        D = Ddg; if strcmp(en,'da'), D = Dda; end
        b = R.fits.(en).B.bMean;

        H  = harmonic_decompose(D.Y, D.tvCont, cfg.(en));
        Hm = harmonic_decompose(predict_harmonic(b, tvD, cfg.(en)), tvD, cfg.(en), b);

        % CONTEXT-FREE NULL: the OTHER experiment's coefficients, pushed through THIS
        % experiment's geometry. Under no context effect the same four weights drive
        % both experiments, so this thin line would coincide with the red one. Because
        % the two experiments hold different terms constant, the null appears here with
        % offset and ripple exchanged -- e.g. in the left column it must equal the red
        % curve at the horizontal meridian and be its exact negative at the vertical
        % meridian. The dashed version rescales by the fitted cross-experiment gain,
        % separating "polar responses are weaker overall" from a reference-frame change.
        if strcmp(en,'dg')
            enOther = 'da';  gDir = mean(R.cross.gain.da2dg, 'omitnan');
        else
            enOther = 'dg';  gDir = mean(R.cross.gain.dg2da, 'omitnan');
        end
        bOther = R.fits.(enOther).B.bMean;
        Hx = harmonic_decompose(predict_harmonic(bOther, tvD, cfg.(en)), tvD, cfg.(en), bOther);

        % Plotted at TWICE the harmonic amplitude. 2A, 2B and 2C are each exactly a
        % pro-minus-con response difference between two conditions (or two pairs), so
        % the reader sees the raw quantity with no mental arithmetic, and the panel
        % offsets and amplitudes are the same 2b numbers reported everywhere else.
        %   dg:  2A = horizontal - vertical    2B = 45deg - 135deg
        %        2C = mean(cardinal) - mean(oblique)
        %   da:  2A = radial - tangential      2B = ccspiral - cspiral
        %        2C = mean(polar-cardinal) - mean(polar-oblique)
        obs   = {2*H.A, 2*H.B, 2*H.C};
        curv  = {2*Hm.Ahat, 2*Hm.Bhat, 2*Hm.Chat};
        xnull = {2*Hx.Ahat, 2*Hx.Bhat, 2*Hx.Chat};
        if cfg.(en).isPolar
            lab = {'radial - tangential  = 2b3 + 2b1cos2\theta_V', ...
                   'ccspiral - cspiral  = -2b1sin2\theta_V', ...
                   'polar-cardinal - polar-oblique  = 2b4 + 2b2cos4\theta_V'};
        else
            lab = {'horizontal - vertical  = 2b1 + 2b3cos2\theta_V', ...
                   '45\circ - 135\circ  = 2b3sin2\theta_V', ...
                   'cardinal - oblique  = 2b2 + 2b4cos4\theta_V'};
        end

        % The unit of inference here is the OBSERVER (n = 8), not the vertex: vertices
        % within a subject are far from independent, and a vertex-level SEM over the
        % ~400-500 vertices per bin understates the real uncertainty by roughly an order
        % of magnitude. So average within subject within bin first, then take the mean
        % and SEM across subjects. This also equalises each observer's contribution to
        % the plotted points, matching the per-subject fits and the observer bootstrap.
        edges = 0:15:360;  ctr = edges(1:end-1) + 7.5;
        [~,~,bin] = histcounts(mod(D.tvCont,360), edges);
        nSubj = numel(cfg.subjects);
        okV   = bin > 0 & D.subj > 0;
        subs  = [bin(okV), D.subj(okV)];
        for k = 1:3
            nexttile((ei-1)*3 + k); hold on;
            P  = accumarray(subs, obs{k}(okV), [numel(ctr) nSubj], @mean, NaN);
            nC = sum(isfinite(P), 2);                   % observers contributing per bin
            m  = mean(P, 2, 'omitnan');
            s  = std(P, 0, 2, 'omitnan') ./ sqrt(nC);
            thin = nC < 3;                              % too few observers to mean anything
            m(thin) = NaN;  s(thin) = NaN;
            cNull = [0.10 0.45 0.80];
            hN  = plot(tvD, xnull{k},      '-',  'Color', cNull, 'LineWidth', 1.0);
            hNG = plot(tvD, gDir*xnull{k}, '--', 'Color', cNull, 'LineWidth', 1.0);
            hD  = errorbar(ctr, m, s, 'o', 'Color',[.35 .35 .35], 'MarkerFaceColor',[.35 .35 .35], ...
                     'MarkerSize',3.5, 'LineWidth',0.6, 'CapSize',0);
            hM  = plot(tvD, curv{k}, 'r-', 'LineWidth', 1.8);
            yline(0,'k:');
            xlim([0 360]); xticks(0:90:360); grid on;
            xlabel('pRF polar angle \theta_V (deg, 0 = right HM)');
            if k==1, ylabel(sprintf('%s\nresponse difference', exps{ei,2})); end
            title(lab{k}, 'FontWeight','normal','FontSize',9);
            if k==1   % one legend per row: enOther and gDir differ between rows
                legend([hD hM hN hNG], {'data', 'fit to this experiment', ...
                        sprintf('null: %s coefficients', enOther), ...
                        sprintf('null \\times gain (%.2f)', gDir)}, ...
                       'Location','best', 'FontSize',7, 'Box','off');
            end
        end
    end
    out = save_fig(fig, cfg, sprintf('harmonic_decomposition_%s', variant));
end

% ========================================================================
function out = fig_coefficients(R, cfg, variant)
% Fit A (thetaV quantised to the wedge centre) vs Fit B (true per-vertex thetaV).
% Their difference is the within-ROI local-orientation artifact, isolated.
    labels = {'horiz vs vert','card vs obl','rad vs tang','polC vs polO'};
    cA = [0.55 0.60 0.68];  cB = [0.15 0.35 0.65];

    fig = figure('Color','w','Position',[100 100 1050 460],'Visible','off');
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Per-vertex harmonic coefficients, 2b with %d%% bootstrap CI (%s betas)', ...
          cfg.ciLevel, variant), 'FontWeight','bold');

    ymax = 0;
    for e = {'dg','da'}
        for f = {'A','B'}
            ymax = max(ymax, max(abs(2*R.fits.(e{1}).(f{1}).ci(:))));
        end
    end
    ymax = 1.15*ymax;

    ttl = {'A: Cartesian gratings (dg)','B: Polar gratings (da)'};
    ens = {'dg','da'};
    for ei = 1:2
        en = ens{ei};
        nexttile(ei); hold on;
        for j = 1:4
            for fi = 1:2
                f = 'AB'; f = f(fi);
                r  = R.fits.(en).(f);
                v  = 2*r.bMean(j);  ci = 2*r.ci(j,:);
                x  = j + (fi-1.5)*0.32;
                c  = cA; if fi==2, c = cB; end
                bar(x, v, 0.28, 'FaceColor', c, 'EdgeColor','none');
                errorbar(x, v, v-ci(1), ci(2)-v, 'k', 'LineWidth',1, 'CapSize',5);
            end
        end
        yline(0,'k-');
        xlim([0.5 4.5]); ylim([-ymax ymax]);
        set(gca,'XTick',1:4,'XTickLabel',labels); xtickangle(18); grid on;
        ylabel('2b  (pro minus con)');
        title(ttl{ei});
        if ei==1
            h = [bar(nan,nan,'FaceColor',cA,'EdgeColor','none'), ...
                 bar(nan,nan,'FaceColor',cB,'EdgeColor','none')];
            legend(h, {'Fit A: wedge centre','Fit B: true pRF angle'}, ...
                   'Location','southeast','FontSize',8);
        end
    end
    out = save_fig(fig, cfg, sprintf('harmonic_coefficients_%s', variant));
end

% ========================================================================
function out = fig_crossprediction(R, Ddg, Dda, cfg, variant)
% Per polar-angle wedge, the four stimuli plotted at their LOCAL orientation:
% observed, against the prediction from the other experiment's coefficients, with a
% RESIDUAL row under each. The residual row is the informative one -- the raw curves
% are dominated by the shared tuning, and the systematic part of the failure is only
% a few hundredths of a percent signal change, invisible on the raw scale.
%
% Predictions are built PER SUBJECT (that subject's coefficients from the other
% experiment, and that subject's fitted cross-experiment gain), so the residual's
% error bars are honest: they carry the uncertainty of the transferred coefficients,
% not just of the observations.
%
% Note on the x-axis: the wedge centres are all multiples of 45 deg and the four
% stimuli are 45 deg apart, so the local orientations always land on {0,45,90,135} in
% both rows. What differs is WHICH stimulus sits at each position -- fixed for dg,
% rotating with thetaV for da (at thetaV=90 the pinwheel is at 90 deg; at thetaV=0 it
% is at 0 deg, because radial on the horizontal meridian is horizontal).
    pairs = {'da','dg'; 'dg','da'};
    nP = numel(cfg.paBins);
    nS = numel(cfg.subjects);

    fig = figure('Color','w','Position',[40 40 1650 900],'Visible','off');
    tl  = tiledlayout(fig,4,nP,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Observed vs cross-predicted local orientation tuning (%s betas)', variant), ...
          'FontWeight','bold');
    subtitle(tl, ['black = observed \pm SEM across observers;   red = prediction from the OTHER ' ...
                  'experiment''s per-observer coefficients \times its fitted gain.   ' ...
                  'Residual rows are observed minus predicted, on an expanded scale.'], 'FontSize',9);

    store = cell(2,1);
    ylObs = 0;  ylRes = 0;
    for k = 1:2
        tgt = pairs{k,1}; src = pairs{k,2};
        Dt = Ddg; if strcmp(tgt,'da'), Dt = Dda; end
        if strcmp(src,'dg'), B = R.cross.bDG; else, B = R.cross.bDA; end

        % Per-observer COEFFICIENTS (so the residual error bars carry the uncertainty
        % of the transferred fit), but a SINGLE group-level gain. Per-observer gains
        % are not usable here: in the da->dg direction they run 0.09 to 4.11 across the
        % eight observers, because that direction divides by the small polar-experiment
        % prediction. Averaging ratios that unstable moves the mean gain from 2.29 to
        % 1.32 and flips the sign of the residual at the vertical meridian. One gain for
        % the whole display is both stable and the stricter null -- a single free
        % parameter rather than eight. HARMONIC_CROSSEXP still reports the per-observer
        % free-gain R2; its gain-equalised coefficient differences were checked against
        % a single group gain and move by less than 0.02.
        Yp = nan(size(Dt.Y));
        for si = 1:nS
            m = Dt.subj == si;
            if ~any(m) || ~all(isfinite(B(si,:))), continue; end
            Yp(m,:) = predict_harmonic(B(si,:), Dt.tvCont(m), cfg.(tgt));
        end
        ok = isfinite(Yp) & isfinite(Dt.Y);
        gGrp = sum(Yp(ok) .* Dt.Y(ok)) / sum(Yp(ok).^2);
        Yp   = gGrp * Yp;

        [~, Mo] = harmonic_roi_roundtrip(Dt.Y, Dt, cfg, cfg.(tgt), 'mean');
        [~, Mp] = harmonic_roi_roundtrip(Yp,   Dt, cfg, cfg.(tgt), 'mean');
        [mo, so] = msem(Mo);
        [mp, ~ ] = msem(Mp);
        [mr, sr] = msem(Mo - Mp);
        store{k} = struct('tgt',tgt,'src',src,'g',gGrp, ...
                          'mo',mo,'so',so,'mp',mp,'mr',mr,'sr',sr);
        ylObs = max([ylObs; abs(mo(:))+abs(so(:)); abs(mp(:))]);
        ylRes = max([ylRes; abs(mr(:))+abs(sr(:))]);
    end
    ylObs = 1.15*ylObs;  ylRes = 1.15*ylRes;

    for k = 1:2
        S = store{k};
        expCfg = cfg.(S.tgt);
        for p = 1:nP
            th  = cfg.paBins(p);
            loc = mod(double(expCfg.oriAngle(:)) + double(expCfg.isPolar)*(th-90), 180);
            [ls, ord] = sort(loc);
            xo   = [ls; ls(1)+180];
            wrap = @(v) [v(ord,p); v(ord(1),p)];

            % --- observed vs predicted ---
            nexttile((2*k-2)*nP + p); hold on;
            errorbar(xo, wrap(S.mo), wrap(S.so), 'ko-', 'MarkerFaceColor','k', ...
                     'MarkerSize',4, 'LineWidth',1, 'CapSize',3);
            plot(xo, wrap(S.mp), 'r.-', 'MarkerSize',12, 'LineWidth',1.4);
            yline(0,'k:');
            xlim([-5 185]); ylim([-ylObs ylObs]); xticks(0:90:180); grid on;
            set(gca,'XTickLabel',[]);
            if p==1
                ylabel(sprintf('%s observed\npredicted from %s (g=%.2f)', S.tgt, S.src, S.g));
            else
                set(gca,'YTickLabel',[]);
            end
            title(sprintf('\\theta_V=%d\\circ', cfg.paBins(p)), 'FontWeight','normal','FontSize',9);

            % --- residual. Only the four measured orientations: unlike the line plots
            % above, a bar repeated at +180 deg would read as a fifth observation. ---
            nexttile((2*k-1)*nP + p); hold on;
            hb = bar(ls, S.mr(ord,p), 0.55, 'FaceColor',[0.35 0.45 0.60], 'EdgeColor','none');
            errorbar(ls, S.mr(ord,p), S.sr(ord,p), 'k', 'LineStyle','none', ...
                     'LineWidth',0.8, 'CapSize',3);
            uistack(hb,'bottom');
            yline(0,'k-');
            xlim([-15 195]); ylim([-ylRes ylRes]); xticks(0:90:180); grid on;
            if p==1
                ylabel(sprintf('%s residual\n(observed - predicted)', S.tgt));
            else
                set(gca,'YTickLabel',[]);
            end
            if k==2, xlabel('local orientation (deg)'); else, set(gca,'XTickLabel',[]); end
        end
    end
    out = save_fig(fig, cfg, sprintf('harmonic_crossprediction_%s', variant));
end

% ------------------------------------------------------------------------
function [m, s] = msem(M)
% Mean and SEM across the subject dimension of an nOri x nPA x nSubj array.
    m = mean(M, 3, 'omitnan');
    n = sum(isfinite(M), 3);
    s = std(M, 0, 3, 'omitnan') ./ sqrt(max(n,1));
end

% ========================================================================
function out = save_fig(fig, cfg, base)
    p = fullfile(cfg.figDir, base);
    exportgraphics(fig, [p '.png'], 'Resolution', 150);
    exportgraphics(fig, [p '.pdf']);
    close(fig);
    out = [p '.png'];
    fprintf('plot_harmonic: wrote %s\n', out);
end
