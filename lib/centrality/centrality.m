function scores = centrality(data, clusters)
    analysis = analize(data, [1:length(clusters); clusters]');
    analysis = cellfun(@(x) x(:, [1,3])', analysis, 'UniformOutput', false);
    analysis = sort([analysis{:}],2);
    scores = analysis(1,:);
end

function Clus = analize(ADJ, A)
    %A is a list of the cluster indices of neurons with their numerical labels. 
    %The numerical labels are in the first column, the corresponding cluster index
    %in the next column.


    n=max(A(:,2)); % To find the number of clusters

    for i=1:n
        M{i}=find(A(:,2)==i);
    end

    for i=1:n
        for j=1:length(M{i})
            for k=1:length(M{i})
                N{i}(j,k)=ADJ(M{i}(j),M{i}(k)); %ADJ is the aggregated matrix
            end
        end
    end

    %N contains list of adjacency matrices for the respective clusters. eigencent 
    %can be used on each of these matrices to find centrality scores within these 
    %clusters.
    
    for i=1:length(N)
        Clus{i}=eigcent(N{i});
    end



    for i=1:length(N)
        for j=1:length(M{i})
            Clus{i}(j,3)=M{i}(Clus{i}(j,2));
        end
    end
end

function y = eigcent(X)

    %eigencentrality calculator: Takes aggregated matrix as input, outputs list
    %of centrality scores in first column, numerical labels of neurons in
    %second column. 

    %List is arranged in increasing order of centrality scores,
    %so most central neuron is placed in the end. 

    %The numerical labels are assigned using alphabetic ordering of the text
    %labels of the neurons. 

    [A B]=eigs(X,1);
    A=abs(A);
    for i=1:length(A)
        A(i,2)=i;
    end

    A=sortrows(A,1);

    y=A;

end