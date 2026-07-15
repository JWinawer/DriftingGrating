function [proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, subset)

% recently edited to add vertical horizontal. Prior versions only included
% radialvstang for da, or "main cardinal" for da and dg

% these labels are the same because the betas were organized this way

% now subset=1 is vertical for DG and radial for DA (considered main /
% non-derived subset)

    % only used for conditions NOT derived (to index medianBOLDpa)
    if strcmp(comparisonName, 'motion_minus_orientation')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = [9, 11]; % radial: in/out
            conConditions = [8, 10]; % tangential: c, cc
        elseif subset==1 && (strcmp(projectName, 'dg'))   
            proConditions = [8, 10]; % right / left -- condition is set to have horizontal as "advantaged" based on natural scene statistics, flipped pro and con labels based on this
            conConditions = [9, 11]; % up / down
        else                        % same for dg, da
            proConditions = 8:11;   % cardinal / polar cardinal
            conConditions = 12:15;  % oblique / polar oblique
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'motion_minus_baseline')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = 20:21; % radial: in / out
            conConditions = 18:19; % tangential: c, cc
        elseif subset==1 && (strcmp(projectName, 'dg'))
            proConditions = 18:19; % right / left                       
            conConditions = 20:21; % up / down
        else                       % same for dg, da
            proConditions = 18:21; % cardinal / polar cardinal
            conConditions = 22:25; % oblique / polar oblique
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'orientation_minus_baseline')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = 27; % radial
            conConditions = 26; % tangential
        elseif subset==1 && (strcmp(projectName, 'dg'))
            proConditions = 26; % horizontal
            conConditions = 27; % vertical
        else                       % same for dg, da
            proConditions = 26:27; % cardinal / polar cardinal
            conConditions = 28:29; % oblique / polar oblique
        end
        allConditions = [proConditions, conConditions];
    end

end