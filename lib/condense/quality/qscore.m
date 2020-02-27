function scores = qscore(data, assigments, varargin)
   fn = 'modularity';
   for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'quality'))
            fn = varargin{i+1};
        end
   end
   
   assigments = num2cell(assigments, 2);

   if strcmp(fn, 'modularity')
        scores = cellfun(@(x) affmodularity(data, x), assigments, 'UniformOutput', false);
   else
        scores = cellfun(@(x) evalclusters(data, x', fn), assigments, 'UniformOutput', false);
        scores = cellfun(@(x) x.CriterionValues, scores, 'UniformOutput', false);
   end
   
   scores = cell2mat(scores);
end