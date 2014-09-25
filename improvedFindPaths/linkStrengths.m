
function [linkStrengths,parts,configs] = findPaths2(mni,numOfParts,len)

	%load(['findPaths2_variables_' num2str(numOfParts) '_' num2str(len)]);

	if len == 1,
		allPaths = zeros(numOfParts,numOfParts,4);
		for k = 1:4, 
			for p = 1:numOfParts
				allPaths(p,mni(p,k),k) = 1; 
			end
		end
		dests = [0,+1;
		        0,-1;
		        +1,0;
		        -1,0];

  elseif len == 2,

    parts = [];
    configs = [];

    perp = [2 1 4 3];

    linkStrengths = zeros(size(mni));

    for p = 1:numOfParts
      for dr = 1:4
        if mni(mni(p,dr),perp(dr)) == p,
          linkStrengths(p,dr) = 2;
        end
      end
    end

  else

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
		%for i = 1:len-3
		%	inds = ismember(cfgs(:,i:i+3),P,'rows');
		%	cfgs = cfgs(~inds,:);
		%end
	  ui = [0 0 1 -1];
  	vi = [1 -1 0 0];
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

	numCfgs = size(cfgs,1);

	ui = [0 0 1 -1];
	vi = [1 -1 0 0];

	dest = [sum(ui(cfgs),2) sum(vi(cfgs),2)];
	dests = unique(dest,'rows');
	numDests = size(dests,1);
	indcs = false(size(cfgs,1),size(dests,1));
	for i = 1:numDests
		indcs(:,i) = ismember(dest,dests(i,:),'rows');
	end

	paths = repmat(struct('paths',[]),numOfParts,numDests);
	%D = cell(numOfParts,i);
	allPaths = zeros(numOfParts,numOfParts,numDests);

	cycleIndex = find(ismember(dests,[0 0],'rows'));
	cycleCfgs = find(indcs(:,cycleIndex));
	linkStrengths = zeros(size(mni));

	Dst = repmat(1:numOfParts,numCfgs,1);
	Dst = Dst(:); init_parts = Dst;
	Dst = uint16(Dst);
	cfgs = uint16(cfgs);
	Cfgs = repmat(cfgs,numOfParts,1);
	stages = uint16(zeros(size(Cfgs)));
	CycleCfgs = repmat(cycleCfgs,1,numOfParts) + repmat(numCfgs*(0:numOfParts-1),size(cycleCfgs,1),1);
	CycleCfgs = CycleCfgs(:);


	mni = uint16(mni);
	linkStrengths = zeros(size(mni));

	for l = 1:len
		stages(:,l) = Dst;
		Dst = mni(sub2ind(size(mni),Dst,Cfgs(:,l)));
	end

	cycles = find(Dst == init_parts);

	corCycles = intersect(CycleCfgs,cycles);
	parts = stages(corCycles,:);
	configs = Cfgs(corCycles,:);

  no_dup = true(size(parts,1),1);
  for i = 1:size(parts,1)
    if length(parts(i,:)) ~= length(unique(parts(i,:)))
      no_dup(i) = false;
    end
  end
  parts = parts(no_dup,:);
  configs = configs(no_dup,:);

	links = sub2ind(size(mni),parts(:),configs(:));
	[u v] = count_unique(links);
	linkStrengths(u) = v;

	%for part = 1:numOfParts
	%	dst = repmat(part,numCfgs,1);
	%	stages = zeros(size(cfgs));
	%	for l = 1:len
	%		stages(:,l) = dst;
	%		dst = mni(sub2ind(size(mni),dst,cfgs(:,l)));
	%	end
	%	%for i = 1:numDests
	%	%	[p q] = count_unique(dst(indcs(:,i)));
	%	%	paths(part,i).paths = [p q];
	%	%	allPaths(part,p,i) = q;
	%	%end
	%
	%	%if cycleIndex
	%		cycles = find(dst == part);
	%		corCycles = intersect(cycleCfgs,cycles);
	%		parts = stages(corCycles,:);
	%		configs = cfgs(corCycles,:);
	%		links = sub2ind(size(mni),parts(:),configs(:));
	%		[u v] = count_unique(links);
	%		linkStrengths(u) = linkStrengths(u) + v;
	%	%else 
	%	%	linkStrengths = [];
	%	%end
	%end

	end
