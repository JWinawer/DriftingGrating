function out = plot_spec_profile(S, figDir, L, V)
% PLOT_SPEC_PROFILE  The polar-angle profile the specification actually fits.
%
%   out = plot_spec_profile(S, figDir)         % S from SPEC_PROFILES
%   out = plot_spec_profile(S, figDir, L)      % L from SPEC_AXIS_LIMITS
%   out = plot_spec_profile(S, figDir, L, V)   % V from SPEC_VARIANTS
%
% ON THE ROI ROUTE. That route fits no per-vertex model, so the curve becomes a STEP
% function: the wedge means, held constant across each 45 deg wedge. That is not a
% drawing convenience -- it is precisely what the ROI route assumes about polar-angle
% structure, and drawing it against the same 15 deg data shows what the assumption
% costs wherever the underlying profile is smooth.
%
% ALL SIX PANELS SHARE ONE SCALE. They are the same units and the same kind of
% quantity -- a pro-minus-con orientation contrast in percent signal change -- and
% which of them is large is the whole point, so they are drawn against one axis. The
% two second-harmonic panels are correspondingly small, which is the finding, not a
% drawing problem.
%
% Rows are the two experiments; the three columns are the only three orientation
% contrasts MEASURABLE at a single polar angle. Four orientations at 45 deg spacing
% give each vertex's demeaned response exactly three degrees of freedom, so a profile
% against continuous theta_V can show three curves and no more
% (../supplement/SUPPLEMENT_harmonic_model.md section S2.3). The fourth coefficient is
% identified ACROSS vertices, from the theta_V modulation of the first and third --
% which is what the tilt of these curves is.
%
% This is the figure that shows why theta_V is continuous rather than binned: binning
% would replace each curve with eight points and put cos(4*theta_V) at only two
% values, which is exactly degenerate once one ROI class is lost.
%
% Points are the across-observer mean of each observer's 15-deg bin mean, with SEM
% across observers (n = 8) -- the observer is the unit of inference, so a vertex-level
% error bar would be roughly an order of magnitude too small.

    if nargin < 3 || isempty(L), L = spec_axis_limits(S); end
    if nargin < 4 || isempty(V)
        V = spec_variants('spec');
        if ~S.hasModel, V = spec_variants('roi'); end
    end
    expn = {'dg','da'};
    expLbl = {'Cartesian gratings (dg)','Polar gratings (da)'};
    col = {[0.16 0.20 0.52],[0.55 0.33 0.12],[0.13 0.47 0.24]};

    fig = figure('Color','w','Position',[60 60 1260 700],'Visible','off');
    tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    if S.hasModel, V.curveLabel = 'fitted four-term harmonic model';
    else,          V.curveLabel = ['wedge means, held constant across each 45' char(176) ' wedge']; end
    title(tl, {sprintf(['The three contrasts measurable at one polar angle  |  %s, ' ...
               '%g-%g%s  |  n = %d'], S.area, S.eccRange(1), S.eccRange(2), char(176), ...
               numel(S.subjects)), ...
               ['\rm\fontsize{9}points = observed 15' char(176) ' bin means ' char(177) ...
               ' SEM across observers;  line = ' V.curveLabel]}, ...
               'Interpreter','tex','FontWeight','bold');

    for ei = 1:2
        E = S.(expn{ei});
        x = E.fine.centres;
        for k = 1:3
            o = E.fine.obs(:,:,k);
            n = sum(isfinite(o),1);
            mu = mean(o,1,'omitnan');
            se = std(o,0,1,'omitnan') ./ max(sqrt(n),1);
            mu(n < 3) = NaN;
            if S.hasModel
                mm = mean(E.fine.mdlDense(:,:,k),1,'omitnan');
                xd = E.fine.denseCentres;
            else
                % Wedge means as a step: constant over each 45 deg wedge, centred on
                % the wedge centre, wrapped so the 0/360 edge closes.
                wc = E.fine.wedgeCentres(:).';
                wv = mean(E.fine.wedge(:,:,k),1,'omitnan');
                edg = [wc - 22.5, wc(end) + 22.5];
                xd = reshape([edg(1:end-1); edg(2:end)], 1, []);
                mm = reshape([wv; wv], 1, []);
                keep = xd >= 0 & xd <= 360;
                xd = xd(keep); mm = mm(keep);
            end

            nexttile((ei-1)*3 + k); hold on;
            yline(0,'-','Color',[0.75 0.75 0.75],'LineWidth',0.8);
            fill([x fliplr(x)], [mu-se fliplr(mu+se)], col{k}, ...
                 'FaceAlpha', 0.16, 'EdgeColor','none');
            plot(x, mu, 'o', 'MarkerSize', 4, 'Color', col{k}, 'MarkerFaceColor', col{k});
            plot(xd, mm, '-', 'Color', col{k}, 'LineWidth', 2);
            xlim([0 360]); ylim(L.prof);
            set(gca,'XTick',0:90:360,'FontSize',9); grid on; box off;
            title(E.fineLbl{k}, 'Interpreter','tex','FontSize',10,'FontWeight','bold');
            if k == 1
                ylabel({expLbl{ei}, '% signal change'}, 'FontSize', 9, 'FontWeight','bold');
            end
            if ei == 2, xlabel(['pRF polar angle (' char(176) ')'], 'FontSize', 9); end
        end
    end

    if ~isfolder(figDir), mkdir(figDir); end
    sfx = '';
    if ~(strcmpi(S.area,'V1') && isequal(S.eccRange(:).', [4 8]))
        sfx = sprintf('_%s_%g-%g', lower(S.area), S.eccRange(1), S.eccRange(2));
    end
    base = fullfile(figDir, sprintf('Figure_5_6_%s_profile%s', V.tag, sfx));
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
    close(fig);
    out = [base '.png'];
    fprintf('plot_spec_profile: wrote %s\n', out);
end
