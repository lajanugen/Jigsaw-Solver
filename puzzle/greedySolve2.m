function [greedyOrder,s,sco,mni,time] = greedySolve(pcs3,numOfParts,cols,rows)


ui = [0 0 1 -1];
vi = [1 -1 0 0];
dirs = [0,+1;
        0,-1;
        +1,0;
        -1,0];

numOfParts = cols*rows;
partSize = size(pcs3,1);

partsExpMat = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        partsExpMat(i,j) = (i - 1) * cols + j;
    end
end

%r2 = mscore_mgc_orig(pcs3);
%N = normalizeScores(r2);
%%[mn,mni] = min(r2,[],2);
%%[~,cost_ranks] = sort(r2,2,'ascend');
%
%s = refine13(r2);
%[mn,mni] = max(s,[],2);
%%numCorPairsmni(false,mni,cols,rows,partsExpMat)
%%pause
%[~,cost_ranks] = sort(s,2,'descend');
%
%
%%mn = squeeze(mn);
%mni = squeeze(mni);
%mni = neighb_ref(mni);

start_time = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% multiple metrics %%%%%%%%%%%%%%%%%%%%%%%
mnis = [];
sco = mscore2_mgc_orig(pcs3);    s = refine13(sco); [~,mni] = max(s,[],2); mni = squeeze(mni); mni = neighb_ref(mni); mnis = cat(3,mnis,squeeze(mni));
N = normalizeScores(sco);
[~,cost_ranks] = sort(s,2,'descend');
sc = sco;
sco = mscore2_ssd(pcs3);         s = refine13(sco); [~,mni] = max(s,[],2); mni = squeeze(mni); mni = neighb_ref(mni); mnis = cat(3,mnis,squeeze(mni));
sco = mscore2_mgc_horver(pcs3);  s = refine13(sco); [~,mni] = max(s,[],2); mni = squeeze(mni); mni = neighb_ref(mni); mnis = cat(3,mnis,squeeze(mni));
sco = mscore2_laplacian(pcs3);   s = refine13(sco); [~,mni] = max(s,[],2); mni = squeeze(mni); mni = neighb_ref(mni); mnis = cat(3,mnis,squeeze(mni));
sco = mscore2_laplacian3(pcs3);  s = refine13(sco); [~,mni] = max(s,[],2); mni = squeeze(mni); mni = neighb_ref(mni); mnis = cat(3,mnis,squeeze(mni));

mni = mnis(:,:,1);
str4 = linkStrengths(mni,numOfParts,4);

strs = [];
for u = 1:5
  strn4 = linkStrengths(mnis(:,:,u),numOfParts,4);
  strn6 = linkStrengths(mnis(:,:,u),numOfParts,6);
  strn = strn4 + strn6;
  strs = cat(3,strs,strn);
end

agg_mni = zeros(size(mni));
agg_mni(:,5) = (1:numOfParts)';

for i = 1:numOfParts
  for j = 1:4
    [mm,mm_ind] = max(strs(i,j,:));
    if mm ~= 0,
      agg_mni(i,j) = mnis(i,j,mm_ind);
    else
      agg_mni(i,j) = mnis(i,j,1);
    end
  end
end
agg_mni = neighb_ref(agg_mni);
str4a = linkStrengths(agg_mni,numOfParts,4);
for i = 1:numOfParts 
  for j = 1:4
    if str4a(i,j) < str4(i,j),
      agg_mni(i,j) = mni(i,j);
    end
  end
end
mni = agg_mni;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mni = mni(:,1:4);

opp = [2 1 4 3];
greedyOrder = [];
done = false(size(mni));

%[~,~,~,~,str4] = findPaths2(mni,numOfParts,4);
%[~,~,~,~,str6] = findPaths2(mni,numOfParts,6);
%[~,~,~,~,str8] = findPaths2(mni,numOfParts,8);
str4 = linkStrengths(mni,numOfParts,4);
str6 = linkStrengths(mni,numOfParts,6);
%str8 = linkStrengths(mni,numOfParts,8);
%str10 = linkStrengths(mni,numOfParts,10);

%[~,~,~,~,str10] = findPaths2(mni,numOfParts,10);
str = str4 + str6;% + str8 + str10;% + str10;
%imagesc(reshape(str(:,3),rows,cols));
%pause
%str4(137,3)
%str6(137,3)
%str(137,3)
%pause


crct = 0;
incr = 0;

it = 1;
ct = 0;
cor_p = [];
cor_r = [];
while true
  it = it + 1;
  Str = str;
  Str(done) = 0;
  [m,ind] = max(Str(:));
  inds = find(Str == m);
  [ux,ur] = ind2sub(size(Str),inds);
  s = sub2ind(size(N),ux,mni(inds),ur);
  [~,mi] = min(N(s));
  [p v r] = ind2sub(size(N),s(mi));
  sum(done(:))
  %if m == 0, break; end
  %if it == 500, break; end
  if m <= 8, break; end
  if ~sum(~done(:)), break; end
  %[p r] = ind2sub(size(mni),ind);
  %if ~isCorrectTrue(p,mni(p,r),r), break; end

  if str4(p,r) >= 4, 
    greedyOrder = [greedyOrder;p r mni(p,r) m];
    ct = ct + 1;
    %if ~isCorrect(p,mni(p,r),r,rows)
    %  ct
    %  [p r mni(p,r)]
    %  displayIncorrect(cor_p,cor_r);
    %  pause
    %else
    %  cor_p = [cor_p p];
    %  cor_r = [cor_r r];
    %end
  end
  %[p r mni(p,r)]

  %if isCorrect(p,mni(p,r),r)
  %  crct = crct + 1;
  %else 
  %  incr = incr + 1;
  %end

  done(p,r) = true;
  done(mni(p,r),opp(r)) = true;

  flg = false;
  for w = 1:2
    if w == 1,  rel = r;      prt = mni(p,r); src_part = p;
    else        rel = opp(r); prt = p;        src_part = mni(p,r);
    end
    cost_ranks((rel-1)*numOfParts*numOfParts + find(cost_ranks(:,:,rel) == prt)) = 0;
    parts = setdiff(find(mni(:,rel) == prt),src_part);
    if ~isempty(parts), flg = true; end
    for t = 1:length(parts)
      part = parts(t);
      nxtBest = find(cost_ranks(part,:,rel));
      if numel(nxtBest) == 0,
        done(part,rel) = true;
        %it
        %pause
      else
        mni(part,rel) = cost_ranks(part,nxtBest(1),rel);
      end
    end
  end

  if flg, 
    %[~,~,~,~,str4] = findPaths2(mni,numOfParts,4);
    %[~,~,~,~,str6] = findPaths2(mni,numOfParts,6);
    %[~,~,~,~,str8] = findPaths2(mni,numOfParts,8);
    str4 = linkStrengths_fast(mni,numOfParts,4);
    str6 = linkStrengths_fast(mni,numOfParts,6);
    %str8 = findPaths2_fast(mni,numOfParts,8);
    %[~,~,~,~,str10] = findPaths2(mni,numOfParts,10);
    str = str4 + str6;% + str8;% + str10;
  end

end




%%[~,~,~,~,str4] = findPaths2(mni,numOfParts,4);
%%[~,~,~,~,str6] = findPaths2(mni,numOfParts,6);
%%[~,~,~,~,str8] = findPaths2(mni,numOfParts,8);
%str4 = findPaths2_cor(mni,done,numOfParts,4);
%str6 = findPaths2_cor(mni,done,numOfParts,6);
%%str8 = findPaths2_cor(mni,done,numOfParts,8);
%
%%[~,~,~,~,str10] = findPaths2(mni,numOfParts,10);
%str = str4 + str6 ;%+ str8;% + str10;
%
%[~,cost_ranks] = sort(s,2,'descend');
%crct = 0;
%incr = 0;
%
%it = 1;
%
%while true
%  it = it + 1;
%  Str = str;
%  Str(done) = 0;
%  [m,ind] = max(Str(:));
%  if m == 0, break; end
%  %if it == 500, break; end
%  %if m <= 4, break; end
%  sum(done(:))
%  if ~sum(~done(:)), break; end
%  [p r] = ind2sub(size(mni),ind);
%  %if ~isCorrectTrue(p,mni(p,r),r), break; end
%  greedyOrder = [greedyOrder;p r mni(p,r) m];
%  %[p r mni(p,r)]
%
%  %if isCorrect(p,mni(p,r),r)
%  %  crct = crct + 1;
%  %else 
%  %  incr = incr + 1;
%  %end
%
%  done(p,r) = true;
%  done(mni(p,r),opp(r)) = true;
%
%  flg = false;
%  for w = 1:2
%    if w == 1,  rel = r;      prt = mni(p,r); src_part = p;
%    else        rel = opp(r); prt = p;        src_part = mni(p,r);
%    end
%    cost_ranks((rel-1)*numOfParts*numOfParts + find(cost_ranks(:,:,rel) == prt)) = 0;
%    parts = setdiff(find(mni(:,rel) == prt),src_part);
%    if ~isempty(parts), flg = true; end
%    for t = 1:length(parts)
%      part = parts(t);
%      nxtBest = find(cost_ranks(part,:,rel));
%      if numel(nxtBest) == 0,
%        done(part,rel) = true;
%      else
%        mni(part,rel) = cost_ranks(part,nxtBest(1),rel);
%      end
%    end
%  end
%
%  if flg, 
%    %[~,~,~,~,str4] = findPaths2(mni,numOfParts,4);
%    %[~,~,~,~,str6] = findPaths2(mni,numOfParts,6);
%    %[~,~,~,~,str8] = findPaths2(mni,numOfParts,8);
%    str4 = findPaths2_cor(mni,done,numOfParts,4);
%    str6 = findPaths2_cor(mni,done,numOfParts,6);
%    %str8 = findPaths2_cor(mni,done,numOfParts,8);
%    %[~,~,~,~,str10] = findPaths2(mni,numOfParts,10);
%    str = str4 + str6;% + str8;% + str10;
%  end
%
%end

time = cputime - start_time;


cor = 0; tot = 0;
mark = 0;
found_mark = false;
for i = 1:size(greedyOrder,1)
  pair = greedyOrder(i,:);
  if isCorrectTrue(pair(1),pair(3),pair(2))
    cor = cor + 1;
  else 
    if ~found_mark,
      mark = i;
      found_mark = true;
    end
  end
  tot = tot + 1;
  %[c r] = ind2sub([cols rows],pair(1)); greedyOrder(i,1) = sub2ind([rows cols],r,c);
  %[c r] = ind2sub([cols rows],pair(3)); greedyOrder(i,3) = sub2ind([rows cols],r,c);
end
%save('greedyOrder','greedyOrder');

[tot cor mark]

function displayIncorrect(incorrect_parts,incorrect_relation) 
  impcs = zeros((partSize+2)*rows,(partSize+4)*cols,3);
  orig_image = imread('../images/3.png');
  for index = 1 : numOfParts

    rowStartIndex = (ceil(index / cols)  - 1) * partSize + 1;
    rowEndIndex   = rowStartIndex + (partSize -  1);
    colStartIndex = mod(index - 1, cols)  * partSize + 1;
    colEndIndex   = colStartIndex + (partSize -  1);

    rowStartIndx  = (ceil(index / cols)  - 1) * (partSize+4) + 3;
    rowEndIndx    = rowStartIndx + (partSize -  1);
    colStartIndx  = mod(index - 1, cols)  * (partSize+4) + 3;
    colEndIndx    = colStartIndx + (partSize -  1);

    impcs(rowStartIndx :rowEndIndx, colStartIndx :colEndIndx, :) = orig_image(rowStartIndex :rowEndIndex, colStartIndex :colEndIndex, :);
  end

  %incorrect_parts = [4,3,28,28];
  %incorrect_relation = [1,2,3,4];

  for index = 1:length(incorrect_parts)
    part = incorrect_parts(index);
    rel  = incorrect_relation(index);
    rowStartIndx  = (ceil(part / cols)  - 1) * (partSize+4) + 3;
    rowEndIndx    = rowStartIndx + partSize-1;
    colStartIndx  = mod(part - 1, cols)  * (partSize+4) + 3;
    colEndIndx    = colStartIndx + partSize-1;
    switch rel
      case 1, impcs(rowStartIndx:rowEndIndx,colStartIndx-2:colStartIndx-1,1) = bitmax;
      case 2, impcs(rowStartIndx:rowEndIndx,colEndIndx+1:colEndIndx+2,2) = bitmax;
      case 3, impcs(rowStartIndx-2:rowStartIndx-1,colStartIndx:colEndIndx,3) = bitmax;
      case 4, 
        impcs(rowEndIndx+1:rowEndIndx+2,colStartIndx:colEndIndx,1) = bitmax;
        impcs(rowEndIndx+1:rowEndIndx+2,colStartIndx:colEndIndx,2) = bitmax;
    end  
  end
  figure,
  imshow(impcs/65536,[0 1]);
  %imwrite(impcs/65536,'b.bmp','bmp');
end


function cor = isCorrectTrue(u,v,r)
  cor = false;
  if (r==1) && (v==u-rows),  cor = true; end
  if (r==2) && (v==u+rows),  cor = true; end
  if (r==3) && (v==u-1), cor = true; end
  if (r==4) && (v==u+1), cor = true; end
end



%function cor = isCorrectTrue(u,v,r)
%  cor = false;
%  if (r==1) && (v==u-1),  cor = true; end
%  if (r==2) && (v==u+1),  cor = true; end
%  if (r==3) && (v==u-cols), cor = true; end
%  if (r==4) && (v==u+cols), cor = true; end
%end

function C = refine13(cost);

  cost = double(cost);
  mask = logical(eye(numOfParts));
  mask = cat(3,mask,mask,mask,mask);
  cost(mask) = 0;  
  %sm = squeeze(sum(cost,2));
  %for i = 1:numOfParts
  %  cost(i,:,1) = cost(i,:,1)/sm(i,1);
  %  cost(i,:,2) = cost(i,:,2)/sm(i,2);
  %  cost(i,:,3) = cost(i,:,3)/sm(i,3);
  %  cost(i,:,4) = cost(i,:,4)/sm(i,4);
  %end
  sm = max(max(max(cost)));
  cost = cost / sm;

  cost(mask) = 1;
  cost = 1 - cost;

  C = cost;
  old_ct = 0;

  %while true,
  for tc = 1:10,

    [mn,mni] = max(C,[],2);
    mn = squeeze(mn);
    mni = squeeze(mni);
    ct = 0;
    for i = 1:numOfParts - 1
      if(mni(mni(mni(mni(i,2),4),1),3) == i), ct = ct + 1; end
      if(mni(mni(mni(mni(i,4),1),3),2) == i), ct = ct + 1; end
      if(mni(mni(mni(mni(i,1),3),2),4) == i), ct = ct + 1; end
      if(mni(mni(mni(mni(i,3),2),4),1) == i), ct = ct + 1; end
    end

    hits = 0;
    %if ct <= old_ct,
    %  C = old_C;
    %  break;
    %end
    %[num] = numCorPairs(true,C);
    %num
    %fprintf('%d,%g%%\n',ct,num);


    cost = C;
    old_C = C;
    old_ct = ct;
  
    L = cost(:,:,1);
    R = cost(:,:,2);
    U = cost(:,:,3);
    D = cost(:,:,4);
    
    %sorted = sort(cost,2,'descend');
    %sorted(:,1:2,1)
    %flag = squeeze(sorted(:,1,:) < 3*sorted(:,2,:));
    
    [mn,mni] = max(cost,[],2);
    mn = squeeze(mn);
    mni = squeeze(mni);
    C = zeros(numOfParts,numOfParts,4);
    %C = cost;
    
    for i = 1:numOfParts
    for j = 1:numOfParts
    if i == j, 
      C(i,j,:) = 0;
      continue;
    end
    
    if L(mni(i,3),mni(j,3)) > L(mni(i,4),mni(j,4))
          C(i,j,1) = cost(i,j,1) + L(mni(i,3),mni(j,3));
    else  C(i,j,1) = cost(i,j,1) + L(mni(i,4),mni(j,4));
    end
    if R(mni(i,3),mni(j,3)) > R(mni(i,4),mni(j,4))
          C(i,j,2) = cost(i,j,2) + R(mni(i,3),mni(j,3));
    else  C(i,j,2) = cost(i,j,2) + R(mni(i,4),mni(j,4));
    end
    if U(mni(i,1),mni(j,1)) > U(mni(i,2),mni(j,2))
          C(i,j,3) = cost(i,j,3) + U(mni(i,1),mni(j,1));
    else  C(i,j,3) = cost(i,j,3) + U(mni(i,2),mni(j,2));
    end
    if D(mni(i,1),mni(j,1)) > D(mni(i,2),mni(j,2))
          C(i,j,4) = cost(i,j,4) + D(mni(i,1),mni(j,1));
    else  C(i,j,4) = cost(i,j,4) + D(mni(i,2),mni(j,2));
    end
    %C(i,j,1) = cost(i,j,1) + (L(mni(i,3),mni(j,3)) + L(mni(i,4),mni(j,4)))/2;
    %C(i,j,2) = cost(i,j,2) + (R(mni(i,3),mni(j,3)) + R(mni(i,4),mni(j,4)))/2;  
    %C(i,j,3) = cost(i,j,3) + (U(mni(i,1),mni(j,1)) + U(mni(i,2),mni(j,2)))/2;
    %C(i,j,4) = cost(i,j,4) + (D(mni(i,1),mni(j,1)) + D(mni(i,2),mni(j,2)))/2;
   
    end
    
    if  ( mni(mni(mni(mni(i,1),3),2),4) == i || ...
          mni(mni(mni(mni(i,1),4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == mni(i,1) || ...
          mni(mni(mni(i,4),1),3) == mni(i,1) || ...
          mni(mni(i,1),2) == i),
        
        C(i,:,1) = cost(i,:,1);
    end
    if  ( mni(mni(mni(mni(i,2),3),1),4) == i || ...
          mni(mni(mni(mni(i,2),4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == mni(i,2) || ...
          mni(mni(mni(i,4),2),3) == mni(i,2) || ...
          mni(mni(i,2),1) == i),
        
        C(i,:,2) = cost(i,:,2);
    end
    if  ( mni(mni(mni(mni(i,3),1),4),2) == i || ...
          mni(mni(mni(mni(i,3),2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == mni(i,3) || ...
          mni(mni(mni(i,2),3),1) == mni(i,3) || ...
          mni(mni(i,3),4) == i),
        
        C(i,:,3) = cost(i,:,3);
    end
    if  ( mni(mni(mni(mni(i,4),1),3),2) == i || ...
          mni(mni(mni(mni(i,4),2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == mni(i,4) || ...
          mni(mni(mni(i,2),4),1) == mni(i,4) || ...
          mni(mni(i,4),3) == i),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    [~,l] = max(C(i,:,1));
    [~,r] = max(C(i,:,2));
    [~,u] = max(C(i,:,3));
    [~,d] = max(C(i,:,4));
  %
    if ~( mni(mni(mni(l,3),2),4) == i || ...
          mni(mni(mni(l,4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == l || ...
          mni(mni(mni(i,4),1),3) == l || ...
          mni(mni(i,1),2) == i),
        
        C(i,:,1) = cost(i,:,1);
    end
    if ~( mni(mni(mni(r,3),1),4) == i || ...
          mni(mni(mni(r,4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == r || ...
          mni(mni(mni(i,4),2),3) == r || ...
          mni(mni(i,2),1) == i),
        
        C(i,:,2) = cost(i,:,2);
    end
    if ~( mni(mni(mni(u,1),4),2) == i || ...
          mni(mni(mni(u,2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == u || ...
          mni(mni(mni(i,2),3),1) == u || ...
          mni(mni(i,3),4) == i),
        
        C(i,:,3) = cost(i,:,3);
    end
    if ~( mni(mni(mni(d,1),3),2) == i || ...
          mni(mni(mni(d,2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == d || ...
          mni(mni(mni(i,2),4),1) == d || ...
          mni(mni(i,4),3) == i),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    end
  end
end

end

function cor = isCorrect(u,v,r,rows)
  cor = false;
  if (r==1) && (v==u-rows),  cor = true; end
  if (r==2) && (v==u+rows),  cor = true; end
  if (r==3) && (v==u-1), cor = true; end
  if (r==4) && (v==u+1), cor = true; end
end


function [perc,incorrect_parts,incorrect_relation] = numCorPairsmni(flag,cost,cols,rows,partsExpMat)
  incorrect_parts = [];
  incorrect_relation = [];
  ht = 0;
  for i = 1:rows
    for j = 1:cols
        if (i > 1),   
          if(wasCompatibilityFunctionCorrectmni(partsExpMat(i,j), partsExpMat(i-1,j), 3, flag, cost)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 3]; 
          end 
        end
        if (i<rows),  
          if(wasCompatibilityFunctionCorrectmni(partsExpMat(i,j), partsExpMat(i+1,j), 4, flag, cost)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 4]; 
          end 
        end
        if (j > 1),   
          if(wasCompatibilityFunctionCorrectmni(partsExpMat(i,j), partsExpMat(i,j-1), 1, flag, cost)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 1]; 
          end 
        end
        if (j<cols),  
          if(wasCompatibilityFunctionCorrectmni(partsExpMat(i,j), partsExpMat(i,j+1), 2, flag, cost)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 2]; 
          end 
        end
    end
  end
  perc = (ht /  (2 * (rows * (cols - 1) + (rows - 1) * cols)))*100;
end

function [hit] = wasCompatibilityFunctionCorrectmni(part, actualNeighbor, relation, flag, cost)
  hit = 0;

  %[part relation actualNeighbor cost(part,relation)]
  %pause
  if (cost(part,relation) == actualNeighbor)
    hit = 1;
  end
end

function normSCO = normalizeScores(partsCompVal)
  t = 0.000000000000001; 
  SCO = partsCompVal;
  normSCO = SCO;
  for ii = 1:1:size(SCO,3) %over each possible arrangement.
      %fprintf('Processing Scores Matrix %d\n',ii);
      %    ii
      [aaa,bbb] = sort(SCO(:,:,ii),2); % sorted over each row.
      rowmins = aaa(:,1:2); %2smallest over each row.
      rowminloc = bbb(:,1); %location of minimum in each row.
      
      [aaa,bbb] = sort(SCO(:,:,ii)); % sorted over each column.
      colmins = aaa(1:2,:); %2smallest over each column.
      colminloc = bbb(1,:); %location of minimum in each column.
      
      for jj = 1:1:size(SCO,1) %over each row.
          values = SCO(jj,:,ii);% the values in the row.
          rowmins(jj,1);
          n1 = values.*0 + rowmins(jj,1); %the minimum for that row...
          n1(rowminloc) = rowmins(jj,2); % the second lowest
          %each position can also be replaced by the smallest nonsame value in  the column...
          n2 = values.*0 + colmins(1,:); %the smallest value in each column.
          % but whereever the row is the same, we use the second lowest value instead
          n2(jj==colminloc)= colmins(2,jj==colminloc);
          nval = (values+t)./(min([n1;n2])+t);
          normSCO(jj,:,ii) = nval;
      end
  end
end
