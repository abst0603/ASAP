function samples = epsilon_dense(data,ep)
samples = data;
%% subsample all points randomly
for i = 1:round(size(samples,1)/5)
    aa = randperm(size(samples,1),1);
    idx = rangesearch(samples,samples(aa,:),ep);
    samples = setdiff(samples,samples(idx{1}(2:end),:),'rows','stable');
end
cnt = 0;
%check the packing condition
while (cnt<size(samples,1))
    cnt = cnt + 1;
    idx = rangesearch(samples,samples(cnt,:),ep);
    samples = setdiff(samples,samples(idx{1}(2:end),:),'rows','stable');
end

%Repulsion forces
samples = epsilon_dense_2(data,samples,ep);

end

