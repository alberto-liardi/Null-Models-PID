% study of PID atom distributions depending on MI value

%%%%%%%%%%%%%%%%%%%%%%%%%%
% color the different lines? maybe not. 
% give interpretation of results
%%%%%%%%%%%%%%%%%%%%%%%%%

% as = "Results_post/Kin2_c20/Sleep_1000/";
% aw = "Results_post/Kin2_c20/Awake_1000/";

as = "Results_post/S2_c20/Sleep_1000/";
aw = "Results_post/S2_c20/Awake_1000/";

subject = "S2";

% Syn_distr_single(aw,'W',"Syn");
% Syn_distr_single(as,'S',"Syn");

Syn_distr_double(aw,as,"Ux",subject);
Syn_distr_double(aw,as,"Uy",subject);
Syn_distr_double(aw,as,"Syn",subject);
Syn_distr_double(aw,as,"Red",subject);


function [] = Syn_distr_single(path,cond,comp,sub)

    if cond == 'W', cond = 'Awake';
    else, cond = 'Sleep'; end

    [Un_X, Un_Y, Red, Syn, MI] = MonkeyPIDLoader(path + "MMI/data.txt");
        
    MI_bins = linspace(min(MI{2}), max(MI{2}), 10);
    [MI_s,ind] = sort(MI{2});
    MI_binned = cell(1,9);
    for s = 1:length(MI_bins)-1
        MI_binned{s} = MI_s(MI_s>=MI_bins(s) & MI_s<MI_bins(s+1));
    end
    MI_binned{end}(end+1) = MI_bins(end);
    
    if comp == "Ux", atom = Un_X{2}(ind); tit = "Unique Information X";
    elseif comp == "Uy", atom = Un_Y{2}(ind); tit = "Unique Information Y";
    elseif comp == "Red", atom = Red{2}(ind); tit = "Redundancy";
    elseif comp == "Syn", atom = Syn{2}(ind); tit = "Synergy"; 
    end
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    last = 1;
    for i = 1:9
        subplot(3,3,i);
        first = last;
        last = first + length(MI_binned{i});
%         disp(first); disp(last);
        edges = linspace(min(atom(first:last-1)), max(atom(first:last-1)), 20);
        histogram(atom(first:last-1), 'BinEdges', edges, "FaceColor", "r");
        title(MI_bins(i) + "-" + MI_bins(i+1));
        mu = mean(atom(first:last-1));
        xline(mu, 'Color', 'black', 'LineWidth', 2);
        legend("\mu= "+sprintf('%.2f',mu),'Location','northwest');
    end
    sgtitle(cond+" "+tit+" "+sub); 
end


function [] = Syn_distr_double(pathw,paths,comp,sub)

    [Un_X_w, Un_Y_w, Red_w, Syn_w, MI_w] = MonkeyPIDLoader(pathw + "MMI/data.txt");
    [Un_X_s, Un_Y_s, Red_s, Syn_s, MI_s] = MonkeyPIDLoader(paths + "MMI/data.txt");

%     ViolinPlot(Un_X_w{2}, Un_Y_w{2}, Red_w{2}, Syn_w{2}, MI_w{2}, sub+" - Var(p) PID atoms Awake", ...
%            pathw + "Violins.png");
%     ViolinPlot(Un_X_s{2}, Un_Y_s{2}, Red_s{2}, Syn_s{2}, MI_s{2}, sub+" - Var(p) PID atoms Sleep", ...
%            paths + "Violins.png");
        
    m = min(min(MI_w{2}), min(MI_s{2}));
    M = max(max(MI_w{2}), max(MI_s{2}));
    MI_bins = linspace(m,M,10);
    [MI_w_sort,ind_w] = sort(MI_w{2});
    [MI_s_sort,ind_s] = sort(MI_s{2});
    
    MI_binned_s = {};
    MI_binned_w = {};
    for s = 1:length(MI_bins)-1
        MI_binned_s{end+1} = MI_s_sort(MI_s_sort>=MI_bins(s) & MI_s_sort<MI_bins(s+1));
        MI_binned_w{end+1} = MI_w_sort(MI_w_sort>=MI_bins(s) & MI_w_sort<MI_bins(s+1));
    end
    if ismember(MI_bins(end), MI_s_sort)
        MI_binned_s{end}(end+1) = MI_bins(end);
    else, MI_binned_w{end}(end+1) = MI_bins(end);
    end

    if comp == "Ux", atom_s = Un_X_s{2}(ind_s); atom_w = Un_X_w{2}(ind_w);
                     tit = "Unique Information X";
    elseif comp == "Uy", atom_s = Un_Y_s{2}(ind_s); atom_w = Un_Y_w{2}(ind_w);
                     tit = "Unique Information Y";
    elseif comp == "Red", atom_s = Red_s{2}(ind_s); atom_w = Red_w{2}(ind_w);
                     tit = "Redundancy";
    elseif comp == "Syn", atom_s = Syn_s{2}(ind_s); atom_w = Syn_w{2}(ind_w);
                     tit = "Synergy";
    end

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    last_w = 1;
    last_s = 1;
    for i = 1:9
        subplot(3,3,i);
        first_s = last_s;
        last_s = first_s + length(MI_binned_s{i});
        first_w = last_w;
        last_w = first_w + length(MI_binned_w{i});
%         disp(first_w); disp(last_w);

        if isempty(atom_s(first_s:last_s-1)) && isempty(atom_w(first_w:last_w-1))
            continue
        elseif isempty(atom_s(first_s:last_s-1))
            edges = linspace(min(atom_w(first_w:last_w-1)), max(atom_w(first_w:last_w-1)), 20);
        elseif isempty(atom_w(first_w:last_w-1))
            edges = linspace(min(atom_s(first_s:last_s-1)), max(atom_s(first_s:last_s-1)), 20);
        else
            edges = linspace(min(min(atom_s(first_s:last_s-1)),min(atom_w(first_w:last_w-1))), ... 
                         max(max(atom_s(first_s:last_s-1)),max(atom_w(first_w:last_w-1))),20);
        end

        mu_s = mean(atom_s(first_s:last_s-1));
        mu_w = mean(atom_w(first_w:last_w-1));

        histogram(atom_s(first_s:last_s-1), 'BinEdges', edges, "FaceColor", "r");
        if ~isempty(atom_s(first_s:last_s-1))
            xline(mu_s, 'Color', 'black', 'LineWidth', 2);
        end
        hold on
        histogram(atom_w(first_w:last_w-1), 'BinEdges', edges, "FaceColor", "g");
        if ~isempty(atom_w(first_w:last_w-1)) 
            xline(mu_w, 'Color', 'black', 'LineWidth', 2);
        end

        if isempty(atom_s(first_s:last_s-1))
            legend("Sedated", "Awake", "\mu_w= "+sprintf('%.2f',mu_w),'Location','northwest');
        elseif isempty(atom_w(first_w:last_w-1))
            legend("Sedated", "\mu_s= "+sprintf('%.2f',mu_s),"Awake", 'Location','northwest');
        else
            legend("Sedated", "\mu_s= "+sprintf('%.2f',mu_s), ...
                   "Awake", "\mu_w= "+sprintf('%.2f',mu_w),'Location','northwest');
        end
        title(MI_bins(i) + "-" + MI_bins(i+1));
    end
    sgtitle(tit+" "+sub); 

end