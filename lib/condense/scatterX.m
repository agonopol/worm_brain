function [X,Y] = scatterX(M, varargin)

    dimensionalityReductionMethod = 'mds';
    sizeAssignment = 20*ones(size(M,1), 1);
    colorAssignment = ones(size(M,1), 1);
    densityColoring = false;
    addJitter = true;
    transparency = 0.5;
    labels = [];
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'dimensionalityReductionMethod'))
            dimensionalityReductionMethod = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'colorAssignment'))
            colorAssignment = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'sizeAssignment'))
            sizeAssignment = sizeAssignment .* varargin{i+1};
        end
        if (strcmp(varargin{i}, 'densityColoring'))
            densityColoring = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'addJitter'))
            addJitter = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'transparency'))
            transparency = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'labels'))
            labels = varargin{i+1};
        end
    end

    rng(122585)
    if (size(M,2) > 3)
        switch (dimensionalityReductionMethod)
            case 'mds'
                D_org = calcDistanceMatrix(M, ...
                                           'k_knn', 'full', ...
                                           'type_k_knn', 'normal', ...
                                           'distfun', 'euclidean', ...
                                           'n_pca', 10);
                try
                    M_toPlot = cmdscale(D_org, 2);
                catch
                     M_toPlot = M;
                end
            case 'tsne'
                %M_toPlot = tsne(M, [], 2, 10, 15);
                %M_toPlot = tsne(M, [], 2, 10, 15);
                M_toPlot = tsne(M);
            end
    else
        M_toPlot = M;
    end
    if addJitter
        M_jitter = rand(size(M_toPlot))-0.5;
        M_jitter = M_jitter * diag((max(M_toPlot)-min(M_toPlot)));
        M_toPlot = M_toPlot + M_jitter / 70;
    end
    X = M_toPlot(:,1);
    Y = M_toPlot(:,2);
    if (size(M_toPlot, 2) == 2)
        if (densityColoring)
            dscatter(M_toPlot(:,1), M_toPlot(:,2));
        else
            s = scatter(M_toPlot(:,1), M_toPlot(:,2), sizeAssignment, colorAssignment, 'filled');
            alpha(s, transparency);
        end
        if ( not(isempty(labels)) )
            dx = 0.1; dy = 0.1;
            x = M_toPlot(:,1) + dx;
            y = M_toPlot(:,2) + dy;
            text(x+dx, y+dy, labels, 'FontSize', 6);
        end
    elseif (size(M_toPlot, 2) == 3)
        s = scatter3(M_toPlot(:,1), M_toPlot(:,2), M_toPlot(:,3), sizeAssignment, colorAssignment, 'filled');
        alpha(s, transparency);
    end
    
end
