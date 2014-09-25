function mni = evalMet(mni);

numOfParts = size(mni,1);
tf = [1 2 3 4];
dr = [0 2 1 3];
mni = [mni (1:numOfParts)'];
%perc = numCorPairsmni(false,mni)
str4 = linkStrengths(mni,numOfParts,4);

Mni = mni;

for maxLen = 7 
%maxLen = 3;
%[cfgs,nonOverLap2,inds2,nonOverLap3,inds3] = findNonOverlap(maxLen);
%[cfgs,inds2,inds3] = findNonOverlap(maxLen);
Cfgs = [];
inds2 = [];
inds3 = [];
load('Cfgs9');
load('inds2_9');
%load('inds3_11');
cfgs = Cfgs;

fprintf('computation done \n');
numCfgs = size(cfgs,1);

trnf = [1 2 3 4;2 1 4 3;3 4 2 1;4 3 1 2];

for iter = 1:10

  %if iter == 5,
  %  [cfgs,inds2] = findNonOverlap(7);
  %  numCfgs = size(cfgs,1);
  %end
  %if iter == 7,
  %  [cfgs,inds2] = findNonOverlap(3);
  %  numCfgs = size(cfgs,1);
  %end

  for j = 1:4

    trf = trnf(j,:);
    trf = [trf 5];
    cfg = trf(cfgs);

    Dst = repmat(1:numOfParts,numCfgs,1);
    Dst = Dst(:);
    cfg = repmat(cfg,numOfParts,1);
    %size(Dst)
    %size(cfg(:,1))
    for t = 1:size(cfgs,2)
      Dst = mni(sub2ind(size(mni),Dst,cfg(:,t)));
    end
    transl_inds = repmat((0:numOfParts-1)*numCfgs,size(inds2,1),1);
    transl_inds = repmat(transl_inds(:),1,2) + repmat(inds2,numOfParts,1);
    dsts = Dst(transl_inds);
    conf_dsts = dsts(:,1) == dsts(:,2);
    src = repmat(1:numOfParts,size(inds2,1),1);
    src = src(:);

    Src = src(conf_dsts);
    Dest = dsts(conf_dsts);

    pairs = [Src Dest];
    enc_pairs = Src + numOfParts*Dest;
    [uniq,cnt] = count_unique(enc_pairs);
    a = rem(uniq,numOfParts);
    b = (uniq - a)/numOfParts;

    %if iter == 8 && j == 4
    %  fi = find(a==218);
    %  [b(fi) cnt(fi)]
    %  pause
    %end

    for i = 1:numOfParts 
      inds = find(a==i);
      if isempty(inds), continue; end
      [~,mx] = max(cnt(inds));
      B = b(inds);
      Mni(i,j) = B(mx);
    end


    %for i = 1:numOfParts 
    %  Dst = repmat(i,numCfgs,1);
   
    %  %if find(Dst == mni(i,j)), continue; end
  
    %  consider = [];
    %  consider_count = [];

    %  parfor k = 1:size(inds2,1) 
    %    ind = inds2(k,:);
    %    dsts = Dst(ind);
    %    if all(dsts == dsts(1))
    %      %Mni(i,j) = dsts(1);
    %      consider = [consider ; dsts(1)];
    %      consider_count = [consider_count sum(find(Dst == dsts(1)))];
    %    end
    %  end

    %  if ~isempty(consider) 
    %    if size(consider,1) == 1, Mni(i,j) = consider(1);
    %    else 
    %      [~,id] = max(consider_count);
    %      Mni(i,j) = consider(id);
    %    end
    %  end
  
    %end
  end

  new_str4 = linkStrengths(Mni,numOfParts,4);
  %Mni(Mni(Mni(218,2),4),1)
  %Mni(Mni(Mni(218,1),4),2)
  %Mni(217,4)
  %mni(217,4)
  %mni(mni(mni(218,2),4),1)
  %mni(mni(mni(218,1),4),2)
  %str4(218,4)
  
  worse_inds = new_str4 <= str4;
  Mni(worse_inds) = mni(worse_inds);
  
  %parfor i = 1:numOfParts 
  %  for j = 1:4
  %    %if str4(i,j) > 0, Mni(i,j) = mni(i,j); continue; end
  %    %if j <= 2
  %    %  if( mni(mni(mni(Mni(i,j),3),3-j),4) ~= i && mni(mni(mni(Mni(i,j),4),3-j),3) ~= i )
  %    %    Mni(i,j) = mni(i,j);
  %    %  end
  %    %else
  %    %  if( mni(mni(mni(Mni(i,j),1),7-j),2) ~= i && mni(mni(mni(Mni(i,j),2),7-j),1) ~= i )
  %    %    Mni(i,j) = mni(i,j);
  %    %  end
  %    %end
  %    if new_str4(i,j) <= str4(i,j)
  %      Mni(i,j) = mni(i,j);
  %    end
  %  end
  %end
  str4 = linkStrengths(Mni,numOfParts,4);

  %for u = 1:numOfParts 
  %  for v = 1:2
  %    if (str4(u,1) >= 4) && (str4(mni(u,1),2) == 0), mni(mni(u,1),2) = u; end
  %    if (str4(u,2) >= 4) && (str4(mni(u,2),1) == 0), mni(mni(u,2),1) = u; end
  %    if (str4(u,3) >= 4) && (str4(mni(u,3),4) == 0), mni(mni(u,3),4) = u; end
  %    if (str4(u,4) >= 4) && (str4(mni(u,4),3) == 0), mni(mni(u,4),3) = u; end
  %  end
  %end

  
  %perc = numCorPairsmni(false,Mni)
  %[perc,incp,incr] = numCorPairsmni(false,Mni);
  %figure
  mni = Mni;
end

%displayIncorrect(incp,incr);

end
