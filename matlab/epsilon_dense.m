function samples = epsilon_dense(data,ep)
%% subsample all points randomly
aa = randperm(size(data,1));
samples = data(aa,:);
cnt = 0;
%check the packing condition
while (cnt<size(samples,1))
    cnt = cnt + 1;
    idx = rangesearch(samples,samples(cnt,:),ep,'SortIndices',false);
    samples = setdiff(samples,samples(idx{1}(2:end),:),'rows','stable');
end

%Repulsion forces
samples = epsilon_dense_2(data,samples,ep);

end

