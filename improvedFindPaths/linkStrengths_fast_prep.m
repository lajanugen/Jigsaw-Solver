
function findPaths2(numOfParts,len)

	mni = zeros(numOfParts,4);

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
	else

	  ui = [0 0 1 -1];
  	vi = [1 -1 0 0];

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


  	dest = [sum(ui(cfgs),2) sum(vi(cfgs),2)];
  	dests = unique(dest,'rows');
  	numDests = size(dests,1);
  	indcs = false(size(cfgs,1),size(dests,1));
  	for i = 1:numDests
  		indcs(:,i) = ismember(dest,dests(i,:),'rows');
  	end

  	paths = repmat(struct('paths',[]),numOfParts,numDests);
  	%D = cell(numOfParts,i);
  	%allPaths = zeros(numOfParts,numOfParts,numDests);

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

	save(['findPaths2_variables_' num2str(numOfParts) '_' num2str(len)], 'Dst', 'init_parts', 'Cfgs', 'CycleCfgs', 'stages');

	end
