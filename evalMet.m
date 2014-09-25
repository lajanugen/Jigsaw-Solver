
%function [scores,aggr_score] = evalMet();
%function [l,p,pr,csts,cors,mn,mni] = evalMet();
%function [sort_cost,rankLim,cdif,md] = evalMet(map);
%function [mni] = evalMet(mni,mniC);
function [perc] = evalMet();
%function [partsCompVal1,partsCompVal2,partsCompVal3,partsCompVal4,partsCompVal5,partsCompVal6,partsCompVal7,partsCompVal8] = evalMetricsOR(compatibility)
compatibility = '';
compat_args = [];
oppSide = [2 1 4 3];
normalMode  = 1;
colorScheme = 0;
noImages    = 1;
accr        = zeros(1,20);
write       = 0;

bestNeighborsMGC = {};
bestNeighborsSSD = {};
pA11 = {};

lst = [];
pairCnt = 0;


partSize    = 14;
partSze     = 2*partSize;

n = 6;
m = 8;
prt_sz = 2*partSize+4;
page_im = zeros(1:(prt_sz+3)*m+3, 1:(prt_sz+3)*n+3, 3);
page_cnt = 1;
page_count = 0;

imNum = 1;
imNo = 1;

scores = cell(1,imNum);
aggr_score = zeros(1,20);

perc = [];
for imNo = 1:20
imNo
hits = 0;

% Reading the image and converting to the new color scheme.
img = strcat(int2str(imNo),'.png');
%img = strcat(int2str(imNo),'.jpg');
path = strcat('/home/hevc/Desktop/Lajan_Research/lajan_puzzle/images', '/', img);
%path = strcat('/home/hevc/Desktop/Lajan_Research/lajan_puzzle/datasets/mcgill', '/', img);

image = imread(path);
%image = image(size(image,1)/2:end,1:size(image,1)/2,:);

orig_image = image;
%image = imresize(image,0.5);
imageWithNewColorScheme = convertToUserSelectedColorScheme();
pcs = {};

image = double(image);
imageSize = size(image);

% Calculate a parts data
rows = floor(imageSize(1) / partSize);
cols = floor(imageSize(2) / partSize);
numOfParts = rows * cols;
partsCorrectOrder = 1:numOfParts;
correct = false(numOfParts,numOfParts,4);

% Creating the expected result matrix so we will be able to
% know a parts neighbors.
partsExpMat = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        partsExpMat(i,j) = (i - 1) * cols + j;
    end
end

edgeMag = cell(numOfParts,4);
edgePos = cell(numOfParts,4);
edges = cell(numOfParts,4);

actualNeighbors = zeros(numOfParts,4);
actualNeighbors(:,1) = (1:numOfParts)-1;
actualNeighbors(:,2) = (1:numOfParts)+1;
actualNeighbors(:,3) = (1:numOfParts)-cols;
actualNeighbors(:,4) = (1:numOfParts)+cols;

rgbPartsArray =  zeros(partSize, partSize, 3, numOfParts);

% Splits the image into parts.
cutImageToParts();

pcs3 = [];
for i = 1:numel(pcs)
	pcs3 = cat(3,pcs3,pcs{i});
end
pcs3 = double(pcs3);

r2 = mscore_mgc_orig(pcs3);   perc1 = numCorPairs(false,r2);
s = refine13_2(r2);             perc2 = numCorPairs(false,s);
perc2
%pause
%[mn,mni] = max(s,[],2); 
%mn = squeeze(mn);
%mni = squeeze(mni);
%mni = neighb_ref(mni);        perc3 = numCorPairsmni(false,mni);
%perc4 = multiple_metrics(imNo);
%
%
%perc = [perc;perc1 perc2 perc3 perc4]  
%%[perc,incorrect_parts,incorrect_relation] = numCorPairsmni(false,mni);
%%perc
%%if ~dbg, displayIncorrect(incorrect_parts,incorrect_relation); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%displayIncorrect(incorrect_parts,incorrect_relation);

end

function [perc,incorrect_parts,incorrect_relation] = numCorPairsmni(flag,cost)
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

  if (cost(part,relation) == actualNeighbor)
    hit = 1;
  end
end

function score = bhattacharyya(ap)
  A = cell(4,1); 
  for edge = [1 2 3 4]
    %A{edge}.exists = false(1,numel(ap));
    A{edge}.grad = zeros(numel(ap),768);
  end
  for w = 1:1:numel(ap)
    P1 = ap{w};
    %class(P1)
    for edge = [1 2 3 4]
    	p = getPlanes(P1,edge,2);
      	[c1,~] = imhist(p(:,:,1));
      	[c2,~] = imhist(p(:,:,2));
      	[c3,~] = imhist(p(:,:,3));
      	A{edge}.grad(w,:) = [c1/sum(c1);c2/sum(c2);c3/sum(c3)];
    end
  end

  N =numel(ap);
  score = zeros(N,N,4); 
  oppSide = [2 1 4 3];
  cnt = 0;
  for jj = 1:1:4
    for ii = 1:1:N-1
      score(ii,ii,:) = inf;             
      for kk = ii+1:1:N
        s = oppSide(jj);
        %Dist = sum(abs(A{jj}.grad(ii,:) - A{s}.grad(kk,:)).^2);
        Dist = 0;
        c1 = A{jj}.grad(ii,:);
        c2 = A{s}.grad(kk,:);
        for ct = 0:2
        	u = c1((1:256)+256*ct);
        	v = c2((1:256)+256*ct);          

          Dist = Dist - log(sum(sqrt(u.*v)));          
        end
        score(ii,kk,jj) = Dist;
        score(kk,ii,oppSide(jj)) = Dist;
      end
    end
  end
  score(N,N,:) = inf;
end

function score = KL(ap)
  A = cell(4,1); 
  for edge = [1 2 3 4]
    %A{edge}.exists = false(1,numel(ap));
    A{edge}.grad = zeros(numel(ap),768);
  end
  for w = 1:1:numel(ap)
    P1 = ap{w};
    %class(P1)
    for edge = [1 2 3 4]
      p = getPlanes(P1,edge,2);
      %p = P1;
        %A{edge}.grad(w,:) = [imhist(p(:,1));imhist(p(:,2));imhist(p(:,3))];
        [c1,~] = imhist(p(:,:,1));
        [c2,~] = imhist(p(:,:,2));
        [c3,~] = imhist(p(:,:,3));
        A{edge}.grad(w,:) = [c1/sum(c1);c2/sum(c2);c3/sum(c3)];
        %pause
    end
  end

  N =numel(ap);
  score = zeros(N,N,4); 
  oppSide = [2 1 4 3];
  cnt = 0;
  for jj = 1:1:4
    for ii = 1:1:N-1
      score(ii,ii,:) = inf;            
      for kk = ii+1:1:N
        s = oppSide(jj);
        %Dist = sum(abs(A{jj}.grad(ii,:) - A{s}.grad(kk,:)).^2);
        Dist = 0;
        c1 = A{jj}.grad(ii,:);
        c2 = A{s}.grad(kk,:);
        for ct = 0:2
          u = c1((1:256)+256*ct);
          v = c2((1:256)+256*ct);
          
          v(v==0) = 0.0000001;
          u(u==0) = 0.0000001;
          %v = v(u~=0);
          %u = u(u~=0);
          %pause
          t1 = u.*log(u./v);
          t2 = v.*log(v./u);
          Dist = Dist + sum(t1) + sum(t2);
        end
        score(ii,kk,jj) = Dist;
        score(kk,ii,oppSide(jj)) = Dist;
      end
    end
  end
  score(N,N,:) = inf;
end

function score = score_mgc1(pcs)
  partSize = size(pcs,1);
  N = size(pcs,3)/3;
  oppSide = [2 1 4 3];

  dummyDiffs = [ 0 0 0 ; 1 1 1; -1 -1 -1; 0 0 1; 0 1 0; 1 0 0 ; -1 0 0 ; 0 -1 0; 0 0 -1];
   
  pix = zeros(N,3*partSize,4);
  D_mu = zeros(N,3,4);
  D_cov = zeros(3,3,N,4);
  for edge = [1 2 3 4] 
    pix(:,:,edge) = zeros(N,partSize*3);           
    D_mu(:,:,edge)  = zeros(N,3); 
    D_cov(:,:,:,edge) = zeros(3,3,N);
  end
  for w = 1:1:N
    P1 = pcs(:,:,3*w-2:3*w);
    P1 = single(P1);
    for edge = [1 2 3 4]
      P1Dif = getPlane(P1,edge,1) - getPlane(P1,edge,2);
      P1D_mu =  mean(P1Dif);
      %P1D_cov = cov(double([P1Dif]));

      Diffs1 = getPlane(P1,edge,2) - getPlane(P1,edge,3);
      Diffs2 = getPlane(P1,edge,3) - getPlane(P1,edge,4);

      X = [P1Dif;Diffs1;Diffs2];
      C = ones(56,4);
      C((1:28),1) = (1:28)';
      C((1:28),2) = (1:28)' + 28;
      C(28+(1:28),1) = (1:28)' + 28;
      C(28+(1:28),2) = (1:28)' + 56;

      %A0 = inv(P1D_cov);
      A0 = eye(3);
      A = ItmlAlg(C, X, A0);
      %A
      %pause
      P1D_cov = A;
      
      pixels = getPlane(P1,edge,1);
      pix(w,:,edge) = pixels(:);
      D_mu(w,:,edge) = P1D_mu;
      D_cov(:,:,w,edge) = P1D_cov;
    end
  end


  score = zeros(N,N,4,'single'); 
  %for ii = 1:1:N, pcs = single(pcs); end
  onemat = ones(partSize,1);

  for jj = 1:1:4
    for ii = 1:1:N-1
      P1D_mu =    D_mu(ii,:,jj);
      P1D_cov =   D_cov(:,:,ii,jj);
      p1S =       pix(ii,:,jj);
      score(ii,ii,:) = inf;             

      for kk = ii+1:1:N
        s = oppSide(jj);

        P2D_mu =    D_mu(kk,:,s);
        P2D_cov =   D_cov(:,:,kk,s);
        p2S =       pix(kk,:,s);

        % now, compute the score:
        P12DIF = p1S-p2S;
        P12DIF = reshape(P12DIF,partSize,3);
        P21DIF = -P12DIF;
        
        D12 = (P12DIF-(onemat*P2D_mu))*(P2D_cov); D12 = sum(D12 .* (P12DIF-(onemat*P2D_mu)),2);
        D21 = (P21DIF-(onemat*P1D_mu))*(P1D_cov); D21 = sum(D21 .* (P21DIF-(onemat*P1D_mu)),2);

        Dist = sum(sqrt(D12)+sqrt(D21));           

        score(ii,kk,jj)           = single(Dist);
        score(kk,ii,oppSide(jj))  = single(Dist);
      end
    end
  end
  score(N,N,:) = inf;
end

%function score = score_mgc1(pcs)
%  partSize = size(pcs,1);
%  N = size(pcs,3)/3;
%  oppSide = [2 1 4 3];
%
%  dummyDiffs = [ 0 0 0 ; 1 1 1; -1 -1 -1; 0 0 1; 0 1 0; 1 0 0 ; -1 0 0 ; 0 -1 0; 0 0 -1];
%   
%  pix = zeros(N,3*partSize,4);
%  D_mu = zeros(N,3,4);
%  D_cov = zeros(3,3,N,4);
%  for edge = [1 2 3 4] 
%    pix(:,:,edge) = zeros(N,partSize*3);           
%    D_mu(:,:,edge)  = zeros(N,3); 
%    D_cov(:,:,:,edge) = zeros(3,3,N);
%  end
%  for w = 1:1:N
%    P1 = pcs(:,:,3*w-2:3*w);
%    P1 = single(P1);
%    for edge = [1 2 3 4]
%      P1Dif = getPlane(P1,edge,1) - getPlane(P1,edge,2);
%      P1D_mu =  mean(P1Dif);
%      %P1D_cov = cov(double([P1Dif]));
%
%      Diffs1 = getPlane(P1,edge,2) - getPlane(P1,edge,3);
%      Diffs2 = getPlane(P1,edge,3) - getPlane(P1,edge,4);
%
%      X = [P1Dif;Diffs1;Diffs2];
%      C = ones(56,4);
%      C((1:28),1) = (1:28)';
%      C((1:28),2) = (1:28)' + 28;
%      C(28+(1:28),1) = (1:28)' + 28;
%      C(28+(1:28),2) = (1:28)' + 56;
%
%      %A0 = inv(P1D_cov);
%      A0 = eye(3);
%      A = ItmlAlg(C, X, A0);
%      %A
%      %pause
%      P1D_cov = A;
%      
%      pixels = getPlane(P1,edge,1);
%      pix(w,:,edge) = pixels(:);
%      D_mu(w,:,edge) = P1D_mu;
%      D_cov(:,:,w,edge) = P1D_cov;
%    end
%  end
%
%
%  score = zeros(N,N,4,'single'); 
%  %for ii = 1:1:N, pcs = single(pcs); end
%  onemat = ones(partSize,1);
%
%  for jj = 1:1:4
%    for ii = 1:1:N-1
%      P1D_mu =    D_mu(ii,:,jj);
%      P1D_cov =   D_cov(:,:,ii,jj);
%      p1S =       pix(ii,:,jj);
%      score(ii,ii,:) = inf;             
%
%      for kk = ii+1:1:N
%        s = oppSide(jj);
%
%        P2D_mu =    D_mu(kk,:,s);
%        P2D_cov =   D_cov(:,:,kk,s);
%        p2S =       pix(kk,:,s);
%
%        % now, compute the score:
%        P12DIF = p1S-p2S;
%        P12DIF = reshape(P12DIF,partSize,3);
%        P21DIF = -P12DIF;
%        
%        D12 = (P12DIF-(onemat*P2D_mu))*(P2D_cov); D12 = sum(D12 .* (P12DIF-(onemat*P2D_mu)),2);
%        D21 = (P21DIF-(onemat*P1D_mu))*(P1D_cov); D21 = sum(D21 .* (P21DIF-(onemat*P1D_mu)),2);
%
%        Dist = sum(sqrt(D12)+sqrt(D21));           
%
%        score(ii,kk,jj)           = single(Dist);
%        score(kk,ii,oppSide(jj))  = single(Dist);
%      end
%    end
%  end
%  score(N,N,:) = inf;
%end

function score = score_mgc2(pcs)
  partSize = size(pcs,1);
  N = size(pcs,3)/3;
  oppSide = [2 1 4 3];

  C = ones(56,4);
  C((1:28),1)     = (1:28)';
  C((1:28),2)     = (1:28)' + 28;
  C(28+(1:28),1)  = (1:28)' + 28*2;
  C(28+(1:28),2)  = (1:28)' + 28*3;

  dummyDiffs = [ 0 0 0 ; 1 1 1; -1 -1 -1; 0 0 1; 0 1 0; 1 0 0 ; -1 0 0 ; 0 -1 0; 0 0 -1];
   
  pix = zeros(N,3*partSize,4);
  D_mu = zeros(N,3,4);
  D_cov = zeros(3,3,N,4);
  diffs = zeros(56,3,N,4);
  for edge = [1 2 3 4] 
    pix(:,:,edge) = zeros(N,partSize*3);           
    D_mu(:,:,edge)  = zeros(N,3); 
    D_cov(:,:,:,edge) = zeros(3,3,N);
    diffs(:,:,:,edge) = zeros(56,3,N);
  end
  for w = 1:1:N
    P1 = pcs(:,:,3*w-2:3*w);
    P1 = single(P1);
    for edge = [1 2 3 4]
      P1Dif = getPlane(P1,edge,1) - getPlane(P1,edge,2);
      P1D_mu =  mean(P1Dif);
      Diffs1 = getPlane(P1,edge,2) - getPlane(P1,edge,3);

      diffs(:,:,w,edge) = [P1Dif;Diffs1];

      pixels = getPlane(P1,edge,1);
      pix(w,:,edge) = pixels(:);
      D_mu(w,:,edge) = P1D_mu;
      %D_cov(:,:,w,edge) = P1D_cov;
    end
  end


  score = zeros(N,N,4,'single'); 
  %for ii = 1:1:N, pcs = single(pcs); end
  onemat = ones(partSize,1);

  for jj = 1:1:4
    for ii = 1:1:N-1
      P1D_mu =    D_mu(ii,:,jj);
      P1D_cov =   D_cov(:,:,ii,jj);
      p1S =       pix(ii,:,jj);
      score(ii,ii,:) = inf;             

      for kk = ii+1:1:N
        s = oppSide(jj);

        P2D_mu =    D_mu(kk,:,s);
        P2D_cov =   D_cov(:,:,kk,s);
        p2S =       pix(kk,:,s);

        % now, compute the score:
        P12DIF = p1S-p2S;
        P12DIF = reshape(P12DIF,partSize,3);
        P21DIF = -P12DIF;

        X = [diffs(:,:,ii,jj);-diffs(:,:,kk,s)];
        P1D_cov = ItmlAlg(C, X, eye(3));

        X = [diffs(:,:,kk,s);-diffs(:,:,ii,jj)];
        P2D_cov = ItmlAlg(C, X, eye(3));

        
        D12 = (P12DIF-(onemat*P2D_mu))*(P2D_cov); D12 = sum(D12 .* (P12DIF-(onemat*P2D_mu)),2);
        D21 = (P21DIF-(onemat*P1D_mu))*(P1D_cov); D21 = sum(D21 .* (P21DIF-(onemat*P1D_mu)),2);

        Dist = sum(sqrt(D12)+sqrt(D21));           

        score(ii,kk,jj)           = single(Dist);
        score(kk,ii,oppSide(jj))  = single(Dist);
      end
    end
  end
  score(N,N,:) = inf;
end

function score = GKL(ap)
  A = cell(4,1); 
  for edge = [1 2 3 4]
    %A{edge}.exists = false(1,numel(ap));
    A{edge}.grad = zeros(numel(ap),768);
  end
  for w = 1:1:numel(ap)
    P1 = ap{w};
    %class(P1)
    for edge = [1 2 3 4]
      p = getPlanes(P1,edge,4);
      %p = P1;
        %A{edge}.grad(w,:) = [imhist(p(:,1));imhist(p(:,2));imhist(p(:,3))];
        [c1,~] = imhist(p(:,:,1));
        [c2,~] = imhist(p(:,:,2));
        [c3,~] = imhist(p(:,:,3));
        A{edge}.grad(w,:) = [c1/sum(c1);c2/sum(c2);c3/sum(c3)];
        %pause
    end
  end

  N =numel(ap);
  score = zeros(N,N,4); 
  oppSide = [2 1 4 3];
  cnt = 0;
  for jj = 1:1:4
    for ii = 1:1:N-1
      score(ii,ii,:) = inf;            
      for kk = 1:1:N
        if (ii==kk) continue; end
        s = oppSide(jj);
        %Dist = sum(abs(A{jj}.grad(ii,:) - A{s}.grad(kk,:)).^2);
        Dist = 0;
        c1 = A{jj}.grad(ii,:);
        c2 = A{s}.grad(kk,:);
        for ct = 0:2
          u = c1((1:256)+256*ct);
          v = c2((1:256)+256*ct);
          
          v(v==0) = 0.0000001;
          v = v(u~=0);
          u = u(u~=0);
          %pause
          t1 = u.*log(u./v);
          Dist = Dist + sum(t1);
        end
        score(ii,kk,jj) = Dist;
      end
    end
  end
  score(N,N,:) = inf;
end

function score = bregman(ap)
  A = cell(4,1); 
  for edge = [1 2 3 4]
    %A{edge}.exists = false(1,numel(ap));
    A{edge}.grad = zeros(numel(ap),768);
  end
  for w = 1:1:numel(ap)
    P1 = ap{w};
    %class(P1)
    for edge = [1 2 3 4]
      p = getPlanes(P1,edge,4);
      %p = P1;
        %A{edge}.grad(w,:) = [imhist(p(:,1));imhist(p(:,2));imhist(p(:,3))];
        [c1,~] = imhist(p(:,:,1));
        [c2,~] = imhist(p(:,:,2));
        [c3,~] = imhist(p(:,:,3));
        A{edge}.grad(w,:) = [c1/sum(c1);c2/sum(c2);c3/sum(c3)];
        %pause
    end
  end

  N =numel(ap);
  score = zeros(N,N,4); 
  oppSide = [2 1 4 3];
  cnt = 0;
  for jj = 1:1:4
    for ii = 1:1:N-1
      score(ii,ii,:) = inf;            
      for kk = 1:1:N
        if (ii==kk) continue; end
        s = oppSide(jj);
        %Dist = sum(abs(A{jj}.grad(ii,:) - A{s}.grad(kk,:)).^2);
        Dist = 0;
        c1 = A{jj}.grad(ii,:);
        c2 = A{s}.grad(kk,:);
        for ct = 0:2
          u = c1((1:256)+256*ct);
          v = c2((1:256)+256*ct);
          
          v(v==0) = 0.0000001;
          v = v(u~=0);
          u = u(u~=0);
          %pause
          t1 = u.*log(u./v);
          Dist = Dist + sum(t1);
        end
        score(ii,kk,jj) = Dist;
      end
    end
  end
  score(N,N,:) = inf;
end

function C = refine8(cost);

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
    [~,num] = numCorPairs(true,C);
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
    
    C(i,j,1) = cost(i,j,1) + (L(mni(i,3),mni(j,3)) + L(mni(i,4),mni(j,4)))/2;
    C(i,j,2) = cost(i,j,2) + (R(mni(i,3),mni(j,3)) + R(mni(i,4),mni(j,4)))/2;  
    C(i,j,3) = cost(i,j,3) + (U(mni(i,1),mni(j,1)) + U(mni(i,2),mni(j,2)))/2;
    C(i,j,4) = cost(i,j,4) + (D(mni(i,1),mni(j,1)) + D(mni(i,2),mni(j,2)))/2;
   
    end
    
    if  ( mni(mni(mni(mni(i,1),3),2),4) == i || ...
          mni(mni(mni(mni(i,1),4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == mni(i,1) || ...
          mni(mni(mni(i,4),1),3) == mni(i,1) ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if  ( mni(mni(mni(mni(i,2),3),1),4) == i || ...
          mni(mni(mni(mni(i,2),4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == mni(i,2) || ...
          mni(mni(mni(i,4),2),3) == mni(i,2) ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if  ( mni(mni(mni(mni(i,3),1),4),2) == i || ...
          mni(mni(mni(mni(i,3),2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == mni(i,3) || ...
          mni(mni(mni(i,2),3),1) == mni(i,3) ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if  ( mni(mni(mni(mni(i,4),1),3),2) == i || ...
          mni(mni(mni(mni(i,4),2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == mni(i,4) || ...
          mni(mni(mni(i,2),4),1) == mni(i,4) ),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    [~,l] = max(C(i,:,1));
    [~,r] = max(C(i,:,2));
    [~,u] = max(C(i,:,3));
    [~,d] = max(C(i,:,4));
  
    if ~( mni(mni(mni(l,3),2),4) == i || ...
          mni(mni(mni(l,4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == l || ...
          mni(mni(mni(i,4),1),3) == l ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if ~( mni(mni(mni(r,3),1),4) == i || ...
          mni(mni(mni(r,4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == r || ...
          mni(mni(mni(i,4),2),3) == r ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if ~( mni(mni(mni(u,1),4),2) == i || ...
          mni(mni(mni(u,2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == u || ...
          mni(mni(mni(i,2),3),1) == u ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if ~( mni(mni(mni(d,1),3),2) == i || ...
          mni(mni(mni(d,2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == d || ...
          mni(mni(mni(i,2),4),1) == d ),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    end
  end
end

function C = refine9(cost);

  cost = double(cost);
  mask = logical(eye(numOfParts));
  mask = cat(3,mask,mask,mask,mask);
  cost(mask) = 0;  
  sm = squeeze(sum(cost,2));
  for i = 1:numOfParts
    cost(i,:,1) = (max(cost(i,:,1)) - cost(i,:,1))/sum(max(cost(i,:,1)) - cost(i,:,1));
    cost(i,:,2) = (max(cost(i,:,2)) - cost(i,:,2))/sum(max(cost(i,:,2)) - cost(i,:,2));
    cost(i,:,3) = (max(cost(i,:,3)) - cost(i,:,3))/sum(max(cost(i,:,3)) - cost(i,:,3));
    cost(i,:,4) = (max(cost(i,:,4)) - cost(i,:,4))/sum(max(cost(i,:,4)) - cost(i,:,4));
  end
  %sm = max(max(max(cost)));
  %cost = cost / sm;
  cost(mask) = 0;

  C = cost;
  old_ct = 0;

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
    [ht,perc,incorrect_parts,incorrect_relation] = numCorPairs(true,C);
    %displayIncorrect(incorrect_parts,incorrect_relation);
    %[~,num] = numCorPairs(true,C);
    fprintf('%d,%g%%\n',ct,perc);


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
    
    C(i,j,1) = cost(i,j,1) + (L(mni(i,3),mni(j,3)) + L(mni(i,4),mni(j,4)))/2;
    C(i,j,2) = cost(i,j,2) + (R(mni(i,3),mni(j,3)) + R(mni(i,4),mni(j,4)))/2;  
    C(i,j,3) = cost(i,j,3) + (U(mni(i,1),mni(j,1)) + U(mni(i,2),mni(j,2)))/2;
    C(i,j,4) = cost(i,j,4) + (D(mni(i,1),mni(j,1)) + D(mni(i,2),mni(j,2)))/2;
   
    %prbR = 0;
    %prbL = 0;
    %prbU = 0;
    %prbD = 0;
    %for u = 1:numOfParts
    %  for v = 1:numOfParts
    %    prbR = prbR + U(i,u)*R(u,v)*D(v,j) + D(i,u)*R(u,v)*U(v,j);
    %    prbL = prbL + U(i,u)*L(u,v)*D(v,j) + D(i,u)*L(u,v)*U(v,j);
    %    prbU = prbU + R(i,u)*U(u,v)*L(v,j) + L(i,u)*U(u,v)*R(v,j);
    %    prbD = prbD + R(i,u)*D(u,v)*L(v,j) + L(i,u)*D(u,v)*R(v,j);
    %  end
    %end

    %C(i,j,1) = cost(i,j,1) + prbL;
    %C(i,j,2) = cost(i,j,2) + prbR;
    %C(i,j,3) = cost(i,j,3) + prbU;
    %C(i,j,4) = cost(i,j,4) + prbD;

    end
    
    if  ( mni(mni(mni(mni(i,1),3),2),4) == i || ...
          mni(mni(mni(mni(i,1),4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == mni(i,1) || ...
          mni(mni(mni(i,4),1),3) == mni(i,1) || ...
          mni(mni(i,4),1) == mni(mni(i,1),3) || ...
          mni(mni(i,3),1) == mni(mni(i,1),4) ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if  ( mni(mni(mni(mni(i,2),3),1),4) == i || ...
          mni(mni(mni(mni(i,2),4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == mni(i,2) || ...
          mni(mni(mni(i,4),2),3) == mni(i,2) || ...
          mni(mni(i,4),2) == mni(mni(i,2),3) || ...
          mni(mni(i,3),2) == mni(mni(i,2),4) ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if  ( mni(mni(mni(mni(i,3),1),4),2) == i || ...
          mni(mni(mni(mni(i,3),2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == mni(i,3) || ...
          mni(mni(mni(i,2),3),1) == mni(i,3) || ...
          mni(mni(i,1),3) == mni(mni(i,3),1) || ...
          mni(mni(i,2),3) == mni(mni(i,3),2) ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if  ( mni(mni(mni(mni(i,4),1),3),2) == i || ...
          mni(mni(mni(mni(i,4),2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == mni(i,4) || ...
          mni(mni(mni(i,2),4),1) == mni(i,4) || ...
          mni(mni(i,1),4) == mni(mni(i,4),1) || ...
          mni(mni(i,2),4) == mni(mni(i,4),2) ),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    [~,l] = max(C(i,:,1));
    [~,r] = max(C(i,:,2));
    [~,u] = max(C(i,:,3));
    [~,d] = max(C(i,:,4));
  
    if ~( mni(mni(mni(l,3),2),4) == i || ...
          mni(mni(mni(l,4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == l || ...
          mni(mni(mni(i,4),1),3) == l ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if ~( mni(mni(mni(r,3),1),4) == i || ...
          mni(mni(mni(r,4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == r || ...
          mni(mni(mni(i,4),2),3) == r ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if ~( mni(mni(mni(u,1),4),2) == i || ...
          mni(mni(mni(u,2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == u || ...
          mni(mni(mni(i,2),3),1) == u ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if ~( mni(mni(mni(d,1),3),2) == i || ...
          mni(mni(mni(d,2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == d || ...
          mni(mni(mni(i,2),4),1) == d ),
        
        C(i,:,4) = cost(i,:,4); 
    end

    [~,Cmni] = max(C,[],2);
    Cmni = squeeze(Cmni);
    for h = 1:numOfParts
      if Cmni(h,1) ~= mni(h,1),
        if ( Cmni(Cmni(Cmni(Cmni(h,1),4),2),3) ~= h || ...
            Cmni(Cmni(Cmni(Cmni(h,1),3),2),4) ~= h || ...
            Cmni(Cmni(h,1),3) ~= Cmni(Cmni(h,3),1) || ...
            Cmni(Cmni(h,1),4) ~= Cmni(Cmni(h,4),1) || ...
            Cmni(h,1) ~= Cmni(Cmni(Cmni(h,3),1),4) || ...
            Cmni(h,1) ~= Cmni(Cmni(Cmni(h,4),1),3) ) && ...
        ( Cmni(Cmni(Cmni(mni(h,1),4),2),3) == h || ...
            Cmni(Cmni(Cmni(mni(h,1),3),2),4) == h || ...
            Cmni(mni(h,1),3) == Cmni(Cmni(h,3),1) || ...
            Cmni(mni(h,1),4) == Cmni(Cmni(h,4),1) || ...
            mni(h,1) == Cmni(Cmni(Cmni(h,3),1),4) || ...
            mni(h,1) == Cmni(Cmni(Cmni(h,4),1),3) ),
          C(h,:,1) = cost(h,:,1);
          Cmni(h,1) = mni(h,1);
        end
      end
      if Cmni(h,2) ~= mni(h,2),
        if ( Cmni(Cmni(Cmni(Cmni(h,2),4),1),3) ~= h || ...
            Cmni(Cmni(Cmni(Cmni(h,2),3),1),4) ~= h || ...
            Cmni(Cmni(h,2),3) ~= Cmni(Cmni(h,3),2) || ...
            Cmni(Cmni(h,2),4) ~= Cmni(Cmni(h,4),2) || ...
            Cmni(h,2) ~= Cmni(Cmni(Cmni(h,3),2),4) || ...
            Cmni(h,2) ~= Cmni(Cmni(Cmni(h,4),2),3) ) && ...
        ( Cmni(Cmni(Cmni(mni(h,2),4),1),3) == h || ...
            Cmni(Cmni(Cmni(mni(h,2),3),1),4) == h || ...
            Cmni(mni(h,2),3) == Cmni(Cmni(h,3),2) || ...
            Cmni(mni(h,2),4) == Cmni(Cmni(h,4),2) || ...
            mni(h,2) == Cmni(Cmni(Cmni(h,3),2),4) || ...
            mni(h,2) == Cmni(Cmni(Cmni(h,4),2),3) ),
          C(h,:,2) = cost(h,:,2);
          Cmni(h,2) = mni(h,2);
        end
      end
      if Cmni(h,3) ~= mni(h,3),
        if ( Cmni(Cmni(Cmni(Cmni(h,3),2),4),1) ~= h || ...
            Cmni(Cmni(Cmni(Cmni(h,3),1),4),2) ~= h || ...
            Cmni(Cmni(h,3),1) ~= Cmni(Cmni(h,1),3) || ...
            Cmni(Cmni(h,3),2) ~= Cmni(Cmni(h,2),3) || ...
            Cmni(h,3) ~= Cmni(Cmni(Cmni(h,1),3),2) || ...
            Cmni(h,3) ~= Cmni(Cmni(Cmni(h,2),3),1) ) && ...
        ( Cmni(Cmni(Cmni(mni(h,3),2),4),1) == h || ...
            Cmni(Cmni(Cmni(mni(h,3),1),4),2) == h || ...
            Cmni(mni(h,3),1) == Cmni(Cmni(h,1),3) || ...
            Cmni(mni(h,3),2) == Cmni(Cmni(h,2),3) || ...
            mni(h,3) == Cmni(Cmni(Cmni(h,1),3),2) || ...
            mni(h,3) == Cmni(Cmni(Cmni(h,2),3),1) ),
          C(h,:,3) = cost(h,:,3);
          Cmni(h,3) = mni(h,3);
        end
      end
      if Cmni(h,4) ~= mni(h,4),
        if ( Cmni(Cmni(Cmni(Cmni(h,4),2),3),1) ~= h || ...
            Cmni(Cmni(Cmni(Cmni(h,4),1),3),2) ~= h || ...
            Cmni(Cmni(h,4),1) ~= Cmni(Cmni(h,1),4) || ...
            Cmni(Cmni(h,4),2) ~= Cmni(Cmni(h,2),4) || ...
            Cmni(h,4) ~= Cmni(Cmni(Cmni(h,1),4),2) || ...
            Cmni(h,4) ~= Cmni(Cmni(Cmni(h,2),4),1) ) && ...
        ( Cmni(Cmni(Cmni(mni(h,4),2),3),1) == h || ...
            Cmni(Cmni(Cmni(mni(h,4),1),3),2) == h || ...
            Cmni(mni(h,4),1) == Cmni(Cmni(h,1),4) || ...
            Cmni(mni(h,4),2) == Cmni(Cmni(h,2),4) || ...
            mni(h,4) == Cmni(Cmni(Cmni(h,1),4),2) || ...
            mni(h,4) == Cmni(Cmni(Cmni(h,2),4),1) ),
          C(h,:,4) = cost(h,:,4);
          Cmni(h,4) = mni(h,4);
        end
      end     


      %if Cmni(h,1) ~= mni(h,1),
      %  if( Cmni(Cmni(Cmni(Cmni(h,1),4),2),3) ~= h || ...
      %      Cmni(Cmni(Cmni(Cmni(h,1),3),2),4) ~= h || ...
      %      Cmni(Cmni(h,1),3) ~= Cmni(Cmni(h,3),1) || ...
      %      Cmni(Cmni(h,1),4) ~= Cmni(Cmni(h,4),1) || ...
      %      Cmni(h,1) ~= Cmni(Cmni(Cmni(h,3),1),4) || ...
      %      Cmni(h,1) ~= Cmni(Cmni(Cmni(h,4),1),3) ),
      %    C(h,:,1) = cost(h,:,1);
      %    Cmni(h,1) = mni(h,1);
      %  end
      %end
      %if Cmni(h,2) ~= mni(h,2),
      %  if( Cmni(Cmni(Cmni(Cmni(h,2),4),1),3) ~= h || ...
      %      Cmni(Cmni(Cmni(Cmni(h,2),3),1),4) ~= h || ...
      %      Cmni(Cmni(h,2),3) ~= Cmni(Cmni(h,3),2) || ...
      %      Cmni(Cmni(h,2),4) ~= Cmni(Cmni(h,4),2) || ...
      %      Cmni(h,2) ~= Cmni(Cmni(Cmni(h,3),2),4) || ...
      %      Cmni(h,2) ~= Cmni(Cmni(Cmni(h,4),2),3) ),
      %    C(h,:,2) = cost(h,:,2);
      %    Cmni(h,2) = mni(h,2);
      %  end
      %end
      %if Cmni(h,3) ~= mni(h,3),
      %  if( Cmni(Cmni(Cmni(Cmni(h,3),2),4),1) ~= h || ...
      %      Cmni(Cmni(Cmni(Cmni(h,3),1),4),2) ~= h || ...
      %      Cmni(Cmni(h,3),1) ~= Cmni(Cmni(h,1),3) || ...
      %      Cmni(Cmni(h,3),2) ~= Cmni(Cmni(h,2),3) || ...
      %      Cmni(h,3) ~= Cmni(Cmni(Cmni(h,1),3),2) || ...
      %      Cmni(h,3) ~= Cmni(Cmni(Cmni(h,2),3),1) ),
      %    C(h,:,3) = cost(h,:,3);
      %    Cmni(h,3) = mni(h,3);
      %  end
      %end
      %if Cmni(h,4) ~= mni(h,4),
      %  if( Cmni(Cmni(Cmni(Cmni(h,4),2),3),1) ~= h || ...
      %      Cmni(Cmni(Cmni(Cmni(h,4),1),3),2) ~= h || ...
      %      Cmni(Cmni(h,4),1) ~= Cmni(Cmni(h,1),4) || ...
      %      Cmni(Cmni(h,4),2) ~= Cmni(Cmni(h,2),4) || ...
      %      Cmni(h,4) ~= Cmni(Cmni(Cmni(h,1),4),2) || ...
      %      Cmni(h,4) ~= Cmni(Cmni(Cmni(h,2),4),1) ),
      %    C(h,:,4) = cost(h,:,4);
      %    Cmni(h,4) = mni(h,4);
      %  end
      %end 

    end
    
    end
  end
end

function C = refine10(cost);

  cost = double(cost);
  mask = logical(eye(numOfParts));
  mask = cat(3,mask,mask,mask,mask);
  cost(mask) = 0;  
  sm = squeeze(sum(cost,2));
  for i = 1:numOfParts
    cost(i,:,1) = (max(cost(i,:,1)) - cost(i,:,1))/sum(max(cost(i,:,1)) - cost(i,:,1));
    cost(i,:,2) = (max(cost(i,:,2)) - cost(i,:,2))/sum(max(cost(i,:,2)) - cost(i,:,2));
    cost(i,:,3) = (max(cost(i,:,3)) - cost(i,:,3))/sum(max(cost(i,:,3)) - cost(i,:,3));
    cost(i,:,4) = (max(cost(i,:,4)) - cost(i,:,4))/sum(max(cost(i,:,4)) - cost(i,:,4));
  end
  %sm = max(max(max(cost)));
  %cost = cost / sm;
  cost(mask) = 0;

  C = cost;
  old_ct = 0;

  while true,

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
    if ct <= old_ct,
      C = old_C;
      break;
    end
    [ht,perc,incorrect_parts,incorrect_relation] = numCorPairs(true,C);
    %displayIncorrect(incorrect_parts,incorrect_relation);
    %[~,num] = numCorPairs(true,C);
    fprintf('%d,%g%%\n',ct,perc);


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
    
    C(i,j,1) = cost(i,j,1) + (L(mni(i,3),mni(j,3)) + L(mni(i,4),mni(j,4)))/2;
    C(i,j,2) = cost(i,j,2) + (R(mni(i,3),mni(j,3)) + R(mni(i,4),mni(j,4)))/2;  
    C(i,j,3) = cost(i,j,3) + (U(mni(i,1),mni(j,1)) + U(mni(i,2),mni(j,2)))/2;
    C(i,j,4) = cost(i,j,4) + (D(mni(i,1),mni(j,1)) + D(mni(i,2),mni(j,2)))/2;

    end
    
    if  ( mni(mni(mni(mni(i,1),3),2),4) == i || ...
          mni(mni(mni(mni(i,1),4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == mni(i,1) || ...
          mni(mni(mni(i,4),1),3) == mni(i,1) || ...
          mni(mni(i,4),1) == mni(mni(i,1),3) || ...
          mni(mni(i,3),1) == mni(mni(i,1),4) ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if  ( mni(mni(mni(mni(i,2),3),1),4) == i || ...
          mni(mni(mni(mni(i,2),4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == mni(i,2) || ...
          mni(mni(mni(i,4),2),3) == mni(i,2) || ...
          mni(mni(i,4),2) == mni(mni(i,2),3) || ...
          mni(mni(i,3),2) == mni(mni(i,2),4) ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if  ( mni(mni(mni(mni(i,3),1),4),2) == i || ...
          mni(mni(mni(mni(i,3),2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == mni(i,3) || ...
          mni(mni(mni(i,2),3),1) == mni(i,3) || ...
          mni(mni(i,1),3) == mni(mni(i,3),1) || ...
          mni(mni(i,2),3) == mni(mni(i,3),2) ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if  ( mni(mni(mni(mni(i,4),1),3),2) == i || ...
          mni(mni(mni(mni(i,4),2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == mni(i,4) || ...
          mni(mni(mni(i,2),4),1) == mni(i,4) || ...
          mni(mni(i,1),4) == mni(mni(i,4),1) || ...
          mni(mni(i,2),4) == mni(mni(i,4),2) ),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    [~,l] = max(C(i,:,1));
    [~,r] = max(C(i,:,2));
    [~,u] = max(C(i,:,3));
    [~,d] = max(C(i,:,4));
  
    if ~( mni(mni(mni(l,3),2),4) == i || ...
          mni(mni(mni(l,4),2),3) == i || ...
          mni(mni(mni(i,3),1),4) == l || ...
          mni(mni(mni(i,4),1),3) == l ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if ~( mni(mni(mni(r,3),1),4) == i || ...
          mni(mni(mni(r,4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == r || ...
          mni(mni(mni(i,4),2),3) == r ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if ~( mni(mni(mni(u,1),4),2) == i || ...
          mni(mni(mni(u,2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == u || ...
          mni(mni(mni(i,2),3),1) == u ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if ~( mni(mni(mni(d,1),3),2) == i || ...
          mni(mni(mni(d,2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == d || ...
          mni(mni(mni(i,2),4),1) == d ),
        
        C(i,:,4) = cost(i,:,4); 
    end

    [~,Cmni] = max(C,[],2);
    Cmni = squeeze(Cmni);
    for h = 1:numOfParts
      if Cmni(h,1) ~= mni(h,1),
        if ( Cmni(Cmni(Cmni(mni(h,1),4),2),3) == h || ...
            Cmni(Cmni(Cmni(mni(h,1),3),2),4) == h ),
          C(h,:,1) = cost(h,:,1);
          Cmni(h,1) = mni(h,1);
        end
      end
      if Cmni(h,2) ~= mni(h,2),
        if ( Cmni(Cmni(Cmni(mni(h,2),4),1),3) == h || ...
            Cmni(Cmni(Cmni(mni(h,2),3),1),4) == h ),
          C(h,:,2) = cost(h,:,2);
          Cmni(h,2) = mni(h,2);
        end
      end
      if Cmni(h,3) ~= mni(h,3),
        if ( Cmni(Cmni(Cmni(mni(h,3),2),4),1) == h || ...
            Cmni(Cmni(Cmni(mni(h,3),1),4),2) == h ),
          C(h,:,3) = cost(h,:,3);
          Cmni(h,3) = mni(h,3);
        end
      end
      if Cmni(h,4) ~= mni(h,4),
        if ( Cmni(Cmni(Cmni(mni(h,4),2),3),1) == h || ...
            Cmni(Cmni(Cmni(mni(h,4),1),3),2) == h ),
          C(h,:,4) = cost(h,:,4);
          Cmni(h,4) = mni(h,4);
        end
      end     
    end
    
    end
  end
end

function C = refine11(cost);

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
    [num] = numCorPairs(true,C);
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
          mni(mni(mni(i,4),1),3) == mni(i,1) ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if  ( mni(mni(mni(mni(i,2),3),1),4) == i || ...
          mni(mni(mni(mni(i,2),4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == mni(i,2) || ...
          mni(mni(mni(i,4),2),3) == mni(i,2) ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if  ( mni(mni(mni(mni(i,3),1),4),2) == i || ...
          mni(mni(mni(mni(i,3),2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == mni(i,3) || ...
          mni(mni(mni(i,2),3),1) == mni(i,3) ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if  ( mni(mni(mni(mni(i,4),1),3),2) == i || ...
          mni(mni(mni(mni(i,4),2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == mni(i,4) || ...
          mni(mni(mni(i,2),4),1) == mni(i,4) ),
        
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
          mni(mni(mni(i,4),1),3) == l ),
        
        C(i,:,1) = cost(i,:,1);
    end
    if ~( mni(mni(mni(r,3),1),4) == i || ...
          mni(mni(mni(r,4),1),3) == i || ...
          mni(mni(mni(i,3),2),4) == r || ...
          mni(mni(mni(i,4),2),3) == r ),
        
        C(i,:,2) = cost(i,:,2);
    end
    if ~( mni(mni(mni(u,1),4),2) == i || ...
          mni(mni(mni(u,2),4),1) == i || ...
          mni(mni(mni(i,1),3),2) == u || ...
          mni(mni(mni(i,2),3),1) == u ),
        
        C(i,:,3) = cost(i,:,3);
    end
    if ~( mni(mni(mni(d,1),3),2) == i || ...
          mni(mni(mni(d,2),3),1) == i || ...
          mni(mni(mni(i,1),4),2) == d || ...
          mni(mni(mni(i,2),4),1) == d ),
        
        C(i,:,4) = cost(i,:,4); 
    end
    
    end
  end
end

function C = refine12(cost);

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
    [num] = numCorPairs(true,C);
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

function C = refine13_2(cost);

  cost = double(cost);
  %mask = logical(eye(numOfParts));
  %mask = cat(3,mask,mask,mask,mask);
  %cost(mask) = 0;  
  %sm = squeeze(sum(cost,2));
  %for i = 1:numOfParts
  %  cost(i,:,1) = cost(i,:,1)/sm(i,1);
  %  cost(i,:,2) = cost(i,:,2)/sm(i,2);
  %  cost(i,:,3) = cost(i,:,3)/sm(i,3);
  %  cost(i,:,4) = cost(i,:,4)/sm(i,4);
  %end
  
  %sm = max(max(max(cost)));
  %cost = cost / sm;

  %cost(mask) = 1;
  %cost = 1 - cost;

  C = cost;
  old_ct = 0;

  %while true,
  for tc = 1:10,

    [mn,mni] = min(C,[],2);
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
    [num] = numCorPairs(false,C);
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
    
    [mn,mni] = min(cost,[],2);
    mn = squeeze(mn);
    mni = squeeze(mni);
    C = zeros(numOfParts,numOfParts,4);
    %C = cost;
    
    for i = 1:numOfParts
    for j = 1:numOfParts
    if i == j, 
      C(i,j,:) = inf;
      continue;
    end
    
    if L(mni(i,3),mni(j,3)) < L(mni(i,4),mni(j,4))
          C(i,j,1) = cost(i,j,1) + L(mni(i,3),mni(j,3));
    else  C(i,j,1) = cost(i,j,1) + L(mni(i,4),mni(j,4));
    end
    if R(mni(i,3),mni(j,3)) < R(mni(i,4),mni(j,4))
          C(i,j,2) = cost(i,j,2) + R(mni(i,3),mni(j,3));
    else  C(i,j,2) = cost(i,j,2) + R(mni(i,4),mni(j,4));
    end
    if U(mni(i,1),mni(j,1)) < U(mni(i,2),mni(j,2))
          C(i,j,3) = cost(i,j,3) + U(mni(i,1),mni(j,1));
    else  C(i,j,3) = cost(i,j,3) + U(mni(i,2),mni(j,2));
    end
    if D(mni(i,1),mni(j,1)) < D(mni(i,2),mni(j,2))
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
    
    [~,l] = min(C(i,:,1));
    [~,r] = min(C(i,:,2));
    [~,u] = min(C(i,:,3));
    [~,d] = min(C(i,:,4));
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
    [num] = numCorPairs(true,C);
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

function [l,path,prts,costs,cors,mn,mni] = findCycles(cost)
  [mn,mni] = min(cost,[],2);
  mn = squeeze(mn);
  mni = squeeze(mni);
  l = [];
  path = [];
  prts = [];
  costs = [];
  cors = [];
  for c1 = 1:4
    for c2 = 1:4
      for c3 = 1:4
        for c4 = 1:4
          for c5 = 1:4
            for c6 = 1:4
              if f(c1) + f(c2) + f(c3) + f(c4) + f(c5) + f(c6) ~= 0, continue; end
              if rem(c1+c2+1,4) == 0, continue; end
              if rem(c2+c3+1,4) == 0, continue; end
              if rem(c3+c4+1,4) == 0, continue; end
              if rem(c4+c5+1,4) == 0, continue; end
              if rem(c5+c6+1,4) == 0, continue; end

              for k = 1:size(mni,1)
                if mni(mni(mni(mni(mni(mni(k,c1),c2),c3),c4),c5),c6) == k,
                  pc = [k mni(k,c1) mni(mni(k,c1),c2) mni(mni(mni(k,c1),c2),c3) mni(mni(mni(mni(k,c1),c2),c3),c4) mni(mni(mni(mni(mni(k,c1),c2),c3),c4),c5)];
                  cfg = [c1 c2 c3 c4 c5 c6];
                  cst = [mn(pc(1),cfg(1)) mn(pc(2),cfg(2)) mn(pc(3),cfg(3)) mn(pc(4),cfg(4)) mn(pc(5),cfg(5)) mn(pc(6),cfg(6))];
                  costs = [costs;cst];
                  prts = [prts;pc];
                  l = [l k];
                  path = [path;cfg];
                  cor = isCorrect(pc(1),pc(2),cfg(1)) && ...
                  isCorrect(pc(2),pc(3),cfg(2)) && ...
                  isCorrect(pc(3),pc(4),cfg(3)) && ...
                  isCorrect(pc(4),pc(5),cfg(4)) && ...
                  isCorrect(pc(5),pc(6),cfg(5)) && ...
                  isCorrect(pc(6),pc(1),cfg(6));
                  cors = [cors cor];
                end
              end
            end
          end
        end
      end
    end
  end
end

function k = f(c)

if c < 3, k = 2*c - 3;
else,     k = (2*c - 7)*i;
end

end

function [count,box] = countBox(cost,flag)
  if flag, [mn,mni] = max(cost,[],2);
  else [mn,mni] = min(cost,[],2);
  end
  
  mn = squeeze(mn);
  mni = squeeze(mni);
  count = 0;
  box = 0;
  for i = 1:numOfParts
    if mni(mni(mni(mni(i,1),3),2),4) == i,
      box = box + 1;
      if ~((mni(i,1) == i-1) && (mni(mni(i,1),3) == i-1-cols) && ...
        (mni(mni(mni(i,1),3),2) == i - cols))
        count = count + 1;
          %i
          %[mni(i,1) mni(mni(i,1),3) mni(mni(mni(i,1),3),2) mni(mni(mni(mni(i,1),3),2),4)]
      end
    end
  end
end

function [n,b] = boxNeighbor(cost)
  [mn,mni] = max(cost,[],2);
  mn = squeeze(mn);
  mni = squeeze(mni);
  n = zeros(size(mni));
  b = false(size(mni));
  for i = 1:numOfParts
    if mni(mni(mni(mni(i,1),3),2),4) == i || mni(mni(mni(mni(i,1),4),2),3) == i,
      n(i,1) = mni(i,1);
      b(i,1) = true;
    end
    if mni(mni(mni(mni(i,2),3),1),4) == i || mni(mni(mni(mni(i,2),4),1),3) == i,
      n(i,2) = mni(i,2);
      b(i,2) = true;
    end
    if mni(mni(mni(mni(i,3),1),4),2) == i || mni(mni(mni(mni(i,3),2),4),1) == i,
      n(i,3) = mni(i,3);
      b(i,3) = true;
    end
    if mni(mni(mni(mni(i,4),1),3),2) == i || mni(mni(mni(mni(i,4),2),3),1) == i,
      n(i,4) = mni(i,4);
      b(i,4) = true;
    end
  end
end

function [p,p_,q,q_,mni] = matches(cost)
  [mn,mni] = max(cost,[],2);
  mn = squeeze(mn);
  mni = squeeze(mni);
  p = zeros(size(mni));
  q = zeros(size(mni));
  p_ = false(size(mni));
  q_ = false(size(mni));
  c = 1;
  d = 1;
  cl = [];
  dl = [];
  for i = 1:numOfParts
    if mni(mni(mni(mni(i,1),3),2),4) == i || mni(mni(mni(mni(i,1),4),2),3) == i,
      p(i,1) = mni(i,1);
      p_(i,1) = true;
      if mni(i,1) == i-1,
        c = c + 1;
      else 
        d = d + 1;
      p_(i,1) = false;
      end
    end
    if mni(mni(mni(mni(i,2),3),1),4) == i || mni(mni(mni(mni(i,2),4),1),3) == i,
      p(i,2) = mni(i,2);
      p_(i,2) = true;
      if mni(i,2) == i+1,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,2) = false
      end
    end
    if mni(mni(mni(mni(i,3),1),4),2) == i || mni(mni(mni(mni(i,3),2),4),1) == i,
      p(i,3) = mni(i,3);
      p_(i,3) = true;
      if mni(i,3) == i-24,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,3) = false
      end
    end
    if mni(mni(mni(mni(i,4),1),3),2) == i || mni(mni(mni(mni(i,4),2),3),1) == i,
      p(i,4) = mni(i,4);
      p_(i,4) = true;
      if mni(i,4) == i+24,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,4) = false
      end
    end
  end
  for i = 1:numOfParts
    if mni(mni(i,1),2) == i,
      q(i,1) = mni(i,1);
      q_(i,1) = true;
      if mni(i,1) == i-1,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,1) = false
      end
    end
    if mni(mni(i,2),1) == i,
      q(i,2) = mni(i,2);
      q_(i,2) = true;
      if mni(i,2) == i+1,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,2) = false
      end
    end
    if mni(mni(i,3),4) == i,
      q(i,3) = mni(i,3);
      q_(i,3) = true;
      if mni(i,3) == i-24,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,3) = false
      end
    end
    if mni(mni(i,4),3) == i,
      q(i,4) = mni(i,4);
      q_(i,4) = true;
      if mni(i,4) == i+24,
        c = c + 1;
      else 
        d = d + 1;
        p_(i,4) = false
      end
    end
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

function normSCO = normalizeScores1(partsCompVal)
  t = 0.000000000000001; 
  SCO = partsCompVal;
  normSCO = SCO;
  for ii = 1:1:size(SCO,3) %over each possible arrangement.
      [aaa,bbb] = sort(SCO(:,:,ii),2); % sorted over each row.
      rowmins = aaa(:,1:2); %2smallest over each row.
      rowminloc = bbb(:,1); %location of minimum in each row.
      
      [aaa,bbb] = sort(SCO(:,:,ii)); % sorted over each column.
      colmins = aaa(1:2,:); %2smallest over each column.
      colminloc = bbb(1,:); %location of minimum in each column.
    for kk = 1:1:size(SCO,2)
      
      for jj = 1:1:size(SCO,1) %over each row.
          


          nval = (SCO(jj,kk,ii)+t)./(max([rowmins(jj,2);colmins(2,kk)])+t);
          normSCO(jj,kk,ii) = nval;
      end
    end
  end
end

function normSCO = normalizeScores2(partsCompVal)
  t = 0.000000000000001; 
  SCO = partsCompVal;
  normSCO = SCO;
  for ii = 1:1:size(SCO,3) %over each possible arrangement.
      [aaa,bbb] = sort(SCO(:,:,ii),2,'descend'); % sorted over each row.
      rowmins = aaa(:,1:2); %2smallest over each row.
      rowminloc = bbb(:,1); %location of minimum in each row.
      
      [aaa,bbb] = sort(SCO(:,:,ii),'descend'); % sorted over each column.
      colmins = aaa(1:2,:); %2smallest over each column.
      colminloc = bbb(1,:); %location of minimum in each column.
    for kk = 1:1:size(SCO,2)
      
      for jj = 1:1:size(SCO,1) %over each row. 


          nval = (SCO(jj,kk,ii)+t)./(min([rowmins(jj,2);colmins(2,kk)])+t);
          normSCO(jj,kk,ii) = nval;
      end
    end
  end
end

function sco = greedyScore(cost,flag)
  num = 0;
  hits = 0;
  while true
    if flag,  [a,b] = max(cost(:));
    else,     [a,b] = min(cost(:));
    end
    if isnan(a), break; end
    [u,v,r] = ind2sub(size(cost),b);
    hits = hits + isCorrect(u,v,r);
    num = num + 1;
    cost(u,:,r) = NaN;
    cost(:,v,r) = NaN;
  end
  sco = 100*double(hits)/double(num);
end

function cor = isCorrect(u,v,r)
  cor = false;
  if (r==1) && (v==u-1),  cor = true; end
  if (r==2) && (v==u+1),  cor = true; end
  if (r==3) && (v==u-24), cor = true; end
  if (r==4) && (v==u+24), cor = true; end

  if (r==1) && (rem(u,24)==1), cor = true; end
  if (r==2) && (rem(u,24)==0), cor = true; end
  if (r==3) && (u<=24), cor = true; end
  if (r==4) && (u>432-24), cor = true; end
end

function [perc,incorrect_parts,incorrect_relation] = numCorPairs(flag,cost)
  incorrect_parts = [];
  incorrect_relation = [];
  ht = 0;
  for i = 1:rows
    for j = 1:cols
        if (i > 1),
        	if(wasCompatibilityFunctionCorrect(partsExpMat(i,j), partsExpMat(i-1,j), 3, flag, cost)), 
            
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 3]; 
        	end 
        end
        if (i<rows),  
        	if(wasCompatibilityFunctionCorrect(partsExpMat(i,j), partsExpMat(i+1,j), 4, flag, cost)), 
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 4]; 
        	end 
        end
        if (j > 1),   
        	if(wasCompatibilityFunctionCorrect(partsExpMat(i,j), partsExpMat(i,j-1), 1, flag, cost)), 
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 1]; 
        	end 
        end
        if (j<cols),  
        	if(wasCompatibilityFunctionCorrect(partsExpMat(i,j), partsExpMat(i,j+1), 2, flag, cost)), 
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

function [ht,perc,incorrect_parts,incorrect_relation] = numCorPairsOR(flag,varargin)
  incorrect_parts = [];
  incorrect_relation = [];
  ht = 0;
  for i = 1:rows
    for j = 1:cols
        if (i > 1),   
        	if(wasCompatibilityFunctionCorrectOR(partsExpMat(i,j), partsExpMat(i-1,j), 3, flag, varargin)), 
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 3]; 
        	end 
        end
        if (i<rows),  
        	if(wasCompatibilityFunctionCorrectOR(partsExpMat(i,j), partsExpMat(i+1,j), 4, flag, varargin)), 
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 4]; 
        	end 
        end
        if (j > 1),   
        	if(wasCompatibilityFunctionCorrectOR(partsExpMat(i,j), partsExpMat(i,j-1), 1, flag, varargin)), 
        		ht = ht + 1; 
        	else, 
        		incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
        		incorrect_relation = [incorrect_relation 1]; 
        	end 
        end
        if (j<cols),  
        	if(wasCompatibilityFunctionCorrectOR(partsExpMat(i,j), partsExpMat(i,j+1), 2, flag, varargin)), 
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

function [hit] = wasCompatibilityFunctionCorrect(part, actualNeighbor, relation, flag, cost)
  hit = 0;
  partsVec = cost(part,:,relation);

  if flag, minNdxVec = find(partsVec==max(partsVec));
  else minNdxVec = find(partsVec==min(partsVec));
  end

  if ((length(minNdxVec) == 1) && minNdxVec == actualNeighbor)
    hit = 1;
  end 
end

function [hit] = wasCompatibilityFunctionCorrectOR(part, actualNeighbor, relation, flag, args)
  hit = 0;

  for nv = 1:length(args)
  	partsVec = args{nv}(part,:,relation);

  	if flag, minNdxVec = find(partsVec==max(partsVec));
  	else minNdxVec = find(partsVec==min(partsVec));
  	end

  	if ((length(minNdxVec) == 1) && minNdxVec == actualNeighbor)
  	  hit = 1;
  	  break;
  	end
  end
end

function [perc,incorrect_parts,incorrect_relation] = numCorPairsMostProb(flag,varargin)
  incorrect_parts = [];
  incorrect_relation = [];
  ht = 0;
  for i = 1:rows
    for j = 1:cols
        if (i > 1),   
          if(wasCompatibilityFunctionCorrectMostProb(partsExpMat(i,j), partsExpMat(i-1,j), 3, flag, varargin)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 3]; 
          end 
        end
        if (i<rows),  
          if(wasCompatibilityFunctionCorrectMostProb(partsExpMat(i,j), partsExpMat(i+1,j), 4, flag, varargin)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 4]; 
          end 
        end
        if (j > 1),   
          if(wasCompatibilityFunctionCorrectMostProb(partsExpMat(i,j), partsExpMat(i,j-1), 1, flag, varargin)), 
            ht = ht + 1; 
          else, 
            incorrect_parts = [incorrect_parts partsExpMat(i,j)]; 
            incorrect_relation = [incorrect_relation 1]; 
          end 
        end
        if (j<cols),  
          if(wasCompatibilityFunctionCorrectMostProb(partsExpMat(i,j), partsExpMat(i,j+1), 2, flag, varargin)), 
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

function sco = mostProbScoPopn(args)
  sco = zeros([numOfParts numOfParts 4]);
  for part = 1:numOfParts
    for relation = 1:4
      neighb = 0;
      bestVec = [];
      mx = 0;
      for nv = 1:5
        partsVec = args{nv}(part,:,relation);
        if args{nv+10}(part,relation),
          neighb = args{nv+5}(part,relation);
          bestVec = partsVec;
        end
      end
      if neighb == 0,
        for nv = 1:5
          partsVec = args{nv}(part,:,relation);      
          if max(partsVec) > mx, 
            mx = max(partsVec);
            bestVec = partsVec;
          end
        end
      end
      sco(part,:,relation) = bestVec;
    end
  end
end


function [hit] = wasCompatibilityFunctionCorrectMostProb(part, actualNeighbor, relation, flag, args)
  hit = 0;

  mx = 0;

  neighb = 0;
  for nv = 1:5
    partsVec = args{nv}(part,:,relation);
    if args{nv+10}(part,relation),
      neighb = args{nv+5}(part,relation);
    end
  end
  if neighb ~= 0, mxVec = neighb;
  else
    for nv = 1:5
      partsVec = args{nv}(part,:,relation);      
      if max(partsVec) > mx, 
        mx = max(partsVec);
        mxVec = find(partsVec == max(partsVec));
      end
    end
  end

  if ((length(mxVec) == 1) && mxVec == actualNeighbor)
    hit = 1;
    %break;
  end
end

function displayIncorrect(incorrect_parts,incorrect_relation) 
  impcs = zeros((partSize+2)*rows,(partSize+4)*cols,3);
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

function compVal = computeCompatibility(first, second, relation)
  p1 = pcs{first}; 
  p2 = pcs{second};

  b1 = getPlane(p1,relation,1);
  b2 = getPlane(p2,oppSide(relation),1);

  switch compatibility
    case 'diff'
      switch compat_args(1)
        case 1, 
                [p q] = compat_args(2:3);
                compVal = sum(sum((abs(b1 - b2)).^p)).^(q/p);
        case 2, 
                compVal = sum(sum((b1 - b2).^2));
        case 3, 
                compVal = sum(sum(sqrt(abs(b1 - b2))));
        case 4, 
                compVal = sum(sum((abs(b1 - b2)).^p)).^(1/p);
      end
    case 'diss'
      b1 = zeros(3,partSize,3);
      b2 = zeros(3,partSize,3);

      switch relation
          case 1 % 'left'
              b1 = transpose3D(p1,2,1);
              b2 = transpose3D(p2,partSize - 1,partSize);
          case 2 % 'right'
              b1 = transpose3D(p1,partSize - 1,partSize);
              b2 = transpose3D(p2,2,1);
          case 3 % 'up'
              b1(1:2,:,:) = p1(2:-1:1,:,:);
              b2(1:2,:,:) = p2(partSize - 1:partSize,:,:);
          case 4 % 'down'
              b1(1:2,:,:) = p1(partSize - 1:partSize,:,:);
              b2(1:2,:,:) = p2(2:-1:1,:,:);
      end

      switch compat_args(1)
          case 1        
              b1(3,:,:) = b1(2,:,:) + (b1(2,:,:) - b1(1,:,:));
              b2(3,:,:) = b2(2,:,:) + (b2(2,:,:) - b2(1,:,:));
              compVal = sum(sum(sum((b1(3,:,:) - b2(2,:,:)).^2))) + sum(sum(sum((b1(2,:,:) - b2(3,:,:)).^2)));
          case 2
              p = 0.3;
              q = 1/32;
              compVal = sum(sum(sum((abs(b1(3,:,:) - b2(2,:,:))).^p))).^(q/p) + sum(sum(sum((abs(b1(2,:,:) - b2(3,:,:))).^p))).^(q/p);
          case 3
              b1(3,:,:) = abs(b1(2,:,:) - b1(1,:,:));
              b2(3,:,:) = abs(b2(2,:,:) - b2(1,:,:));
              compVal = sum(sum(sum(((b1(3,:,:) + b2(3,:,:)) / 2 - abs(b1(2,:,:) - b2(2,:,:))).^2)));
          case 4
              b1(3,:,:) = (b2(2,:,:) - b1(1,:,:)) / 2 + b1(2,:,:);
              b2(3,:,:) = (b2(1,:,:) - b1(2,:,:)) / 2 + b2(2,:,:);
              compVal = sum(sum(sum((b1(3,:,:) - b2(3,:,:)).^2)));
          case 5
              b1(3,:,:) = (b2(2,:,:) - b1(1,:,:)) / 2 + b1(2,:,:);
              b2(3,:,:) = (b2(1,:,:) - b1(2,:,:)) / 2 + b2(2,:,:);
              compVal = sum(sum(sum(sqrt(abs(b1(3,:,:) - b2(3,:,:))))));
      end
  end
end

function [partsCompVal] = initializePartsCompatibility()
  if length(compatibility) == 4 && (all(compatibility == 'diff') || all(compatibility == 'diss'))
    for i = 1:numOfParts
        for j = 1:i
            for l = 1:4
                if (i == j),    partsCompVal(i,j,l) = bitmax;
                else
                  partsCompVal(i,j,l) = computeCompatibility(i, j, l);
                  if (l < 3) partsCompVal(j,i,3-l) = partsCompVal(i,j,l); 
                  else       partsCompVal(j,i,7-l) = partsCompVal(i,j,l); 
                  end
                end
            end
        end
    end
  else
    switch compatibility
    	case 'mgc'
        	partsCompVal = score_mgc(pcs3);		
    end
  end  
end

function plane = getPlane(P1,edge,n)
  if      (edge==3), plane = squeeze(P1(n,:,:));
  elseif  (edge==2), plane = squeeze(P1(:,end+1-n,:));
  elseif  (edge==4), plane = squeeze(P1(end+1-n,:,:));
  else, plane = squeeze(P1(:,n,:));
  end
end

function plane = getPlanes(P1,edge,n)
  if      (edge==3), plane = squeeze(P1(1:n,:,:));
  elseif  (edge==2), plane = squeeze(P1(:,end:-1:end+1-n,:));
  elseif  (edge==4), plane = squeeze(P1(end:-1:end+1-n,:,:));
  elseif  (edge==1), plane = squeeze(P1(:,1:n,:));
  end
end

function [convertedImage] = convertToUserSelectedColorScheme()
  if (colorScheme == 0)
    convertedImage = image;
  elseif (colorScheme == 1)
    cformRgbToLab = makecform('srgb2lab');
    convertedImage = double(applycform(image, cformRgbToLab));
  elseif (colorScheme == 2)
    convertedImage = rgb2hsv(image);
  elseif (colorScheme == 3)
    convertedImage = rgb2gray(image);
  else
    assert(false, strcat('Color scheme:',num2str(colorScheme),' is not supported.'));
  end
end
% Cuts the images into parts.
function cutImageToParts()
  for index = 1 : numOfParts
    rowStartIndex = (ceil(index / cols)  - 1) * partSize + 1;
    rowEndIndex = rowStartIndex + (partSize -  1);
    colStartIndex = mod(index - 1, cols)  * partSize + 1;
    colEndIndex = colStartIndex + (partSize -  1);
    pcs{index} = imageWithNewColorScheme(rowStartIndex :rowEndIndex, colStartIndex :colEndIndex, :);
    rgbPartsArray(:,:,:, index) = image(rowStartIndex :rowEndIndex, colStartIndex :colEndIndex, :);
    %save('rgbParts.mat','rgbPartsArray','-mat')
  end
end

end

