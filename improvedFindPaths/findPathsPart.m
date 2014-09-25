function [stages,dests,cfgs] = findPathsPart(mni,part,Len)

%L = 2*Len;
	ui = [0 0 1 -1 0];
  vi = [1 -1 0 0 0];

Cfgs = [];
for len = 1:Len
	
  cfgs = [1;2;3;4];

	for i = 1:len-1
		cfgs = [[repmat(1,size(cfgs,1),1);repmat(2,size(cfgs,1),1);repmat(3,size(cfgs,1),1);repmat(4,size(cfgs,1),1);] repmat(cfgs,4,1)];
	end

	for i = 1:len-1
		inds = ismember(cfgs(:,i:i+1),[1 2;2 1;3 4;4 3],'rows');
		cfgs = cfgs(~inds,:);
	end

	if len > 4
		%P = perms([1 2 3 4]);
    U = ui(cfgs);
    V = vi(cfgs);
    cyc = false(size(cfgs,1),1);
    for k = 3:2:2*(floor((len-1)/2))-1
		  for i = 1:len-k
        cyc = cyc | ((sum(U(:,i:i+k),2) == 0) & (sum(V(:,i:i+k),2) == 0));
		  end
    end
    cfgs = cfgs(~cyc,:);
	end

  %U = cumsum(ui(cfgs),2);
  %V = cumsum(vi(cfgs),2);
  %U = U + Len;
  %V = V + Len;
  %U = [repmat(0,size(U,1),1) U];
  %V = [repmat(0,size(V,1),1) V];
  %links = U(:,1:end-1) + V(:,1:end-1)*L + U(:,2:end)*L^2 + V(:,2:end)*L^3;
  cfgs = [cfgs repmat(5,size(cfgs,1),Len - size(cfgs,2))];
  %links = [links repmat(0,size(links,1),Len - size(links,2))];
  Cfgs = [Cfgs;cfgs];
  %Links = [Links ; links];

end

  cfgs = Cfgs;
	numCfgs = size(cfgs,1);


  dests = [sum(reshape(ui(cfgs),size(cfgs)),2) sum(reshape(vi(cfgs),size(cfgs)),2)];
  %size(dests)

	%dest = [sum(ui(cfgs),2) sum(vi(cfgs),2)];
	%dests = unique(dest,'rows');
	%numDests = size(dests,1);
	%indcs = false(size(cfgs,1),size(dests,1));
	%for i = 1:numDests
	%	indcs(:,i) = ismember(dest,dests(i,:),'rows');
	%end

	%paths = repmat(struct('paths',[]),numOfParts,numDests);
	%D = cell(numOfParts,i);
	%allPaths = zeros(numOfParts,numOfParts,numDests);

	%cycleIndex = find(ismember(dests,[0 0],'rows'));
	%cycleCfgs = find(indcs(:,cycleIndex));
	%linkStrengths = zeros(size(mni));
	%for part = 1:numOfParts
		dst = repmat(part,numCfgs,1);
		stages = zeros(size(cfgs));
		for l = 1:Len
			stages(:,l) = dst;
			dst = mni(sub2ind(size(mni),dst,cfgs(:,l)));
		end
    stages = [stages dst];
		%for i = 1:numDests
		%	[p q] = count_unique(dst(indcs(:,i)));
		%	paths(part,i).paths = [p q];
		%	%D{part,i} = [D{part,i};[p q]];
		%	allPaths(part,p,i) = q;
		%end

		%if cycleIndex
		%	cycles = find(dst == part);
		%	corCycles = intersect(cycleCfgs,cycles);
		%	parts = stages(corCycles,:);
		%	configs = cfgs(corCycles,:);
    %  cycleParts = [cycleParts;parts];
    %  cycleConfigs = [cycleConfigs;configs];
		%	links = sub2ind(size(mni),parts(:),configs(:));
		%	[u v] = count_unique(links);
		%	linkStrengths(u) = linkStrengths(u) + v;
		%else 
		%	linkStrengths = [];
		%end
	%end


