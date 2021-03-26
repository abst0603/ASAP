function samples1 = epsilon_dense_2(data,samples1,ep)
samples = data;

%% in this parts points repel each other so they will fill the whole space 
% and more space will be provided so we can sample again in this space. The
% final thing will be more dense and the unwanted wholes will be removed
lr = 1;
counter1 = 1;
tau = 3.5;
while(lr > 0.1*ep^3)
    lr = ((ep)^3)*exp(-counter1/tau);
    counter1 = counter1 + 1;
    samples = data;
    setloop = randperm(size(samples1,1));
    idx = rangesearch(samples1,samples1,2*ep);
    for i = setloop
%         idx = rangesearch(samples1,samples1(i,:),2*ep);
        idx{i}(1)=[];
        if(isempty(idx{i})) 
            continue;
        end
        dist = pdist2(samples1(idx{i},:),samples1(i,:),'squaredeuclidean');
        revdist = lr./dist;
        nn = sum(((samples1(i,:) - samples1(idx{i},:))./vecnorm((samples1(i,:) - samples1(idx{i},:)),2,2)).*revdist);
        bb = samples1(i,:) + nn;
        cc = data(knnsearch(data,bb,'k',1),:);
        if((nnz(ismember(samples1,cc,'row'))==0) && nnz(pdist2(samples1(idx{i},:),cc)<=ep)==0)% (idx{i},:)
            samples1(i,:) = cc;
        end
    end
    %% Check to see if some empty space is added or not with respect to previous
    % step
%% version 2(faster)
    idx = rangesearch(samples,samples1,ep,'SortIndices',false);
    innd = unique(horzcat(idx{:}));
    samples(innd,:) = [];
    samples2 = samples;
    %% check the new point by themselves to make sure the epsilon-distance is
    % true in this set too
%% version 2(faster)
    cnt = 0;
    while (cnt<size(samples2,1))
        cnt = cnt + 1;
        idx = rangesearch(samples2,samples2(cnt,:),ep,'SortIndices',false);
        samples2 = setdiff(samples2,samples2(idx{1}(2:end),:),'rows','stable');
    end
    samples1 = [samples1;samples2];
end

cnt = 0;
while (cnt<size(samples1,1))
    cnt = cnt + 1;
    idx = rangesearch(samples1,samples1(cnt,:),ep,'SortIndices',false);
    samples1 = setdiff(samples1,samples1(idx{1}(2:end),:),'rows','stable');
end
end
