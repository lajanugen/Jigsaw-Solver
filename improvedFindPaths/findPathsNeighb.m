function [orig,Dst,stages,dest] = findPaths2(mni,numOfParts,len)

		dests = [0,+1;
		        0,-1;
		        +1,0;
		        -1,0];

	ui = [0 0 1 -1];
  vi = [1 -1 0 0];

	cfgs = [1;2;3;4];

	for i = 1:len-1
		cfgs = [[repmat(1,size(cfgs,1),1);repmat(2,size(cfgs,1),1);repmat(3,size(cfgs,1),1);repmat(4,size(cfgs,1),1);] repmat(cfgs,4,1)];
	end

  dest = [sum(ui(cfgs),2) sum(vi(cfgs),2)];
  cfgs = cfgs(ismember(dest,dests,'rows'),:);

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

  dest = [sum(ui(cfgs),2) sum(vi(cfgs),2)];
	numCfgs = size(cfgs,1);

  dests = unique(dest,'rows');
  numDests = size(dests,1);
  indcs = false(size(cfgs,1),size(dests,1));
  for i = 1:numDests
  	indcs(:,i) = ismember(dest,dests(i,:),'rows');
  end


  Dst = repmat(1:numOfParts,numCfgs,1);
  Cfgs = repmat(cfgs,numOfParts,1);
  [~,k] = ismember(dest,dests,'rows');
  Indcs = repmat(k,numOfParts,1);

  Dst = Dst(:);
  orig = Dst;
  stages = Dst;
	for l = 1:len
		Dst = mni(sub2ind(size(mni),Dst,Cfgs(:,l)));
    stages = [stages Dst];
	end
  stages = stages(:,1:end-1) + stages(:,2:end)*numOfParts;

