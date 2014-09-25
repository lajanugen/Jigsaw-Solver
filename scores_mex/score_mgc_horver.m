function score = score_mgc_horver(pcs)

  partSize = size(pcs,1);
  N = size(pcs,3)/3;
  oppSide = [2 1 4 3];
  dummyDiffs = [ 0 0 0 0 0 0; 1 1 1 1 1 1; -1 -1 -1 -1 -1 -1; 0 0 1 0 0 0; 0 0 0 0 0 1; 0 1 0 0 0 0; 0 0 0 0 1 0; 1 0 0 0 0 0; 0 0 0 1 0 0; -1 0 0 0 0 0; 0 0 0 -1 0 0; 0 -1 0 0 0 0; 0 0 0 0 -1 0; 0 0 -1 0 0 0; 0 0 0 0 0 -1];

  pix = zeros(N,3*partSize,4);
  D_mu = zeros(N,6,4);
  D_cov = zeros(6,6,N,4);
  for edge = [1 2 3 4] 
    pix(:,:,edge) = zeros(N,partSize*3);           
    D_mu(:,:,edge)  = zeros(N,6); 
    D_cov(:,:,:,edge) = zeros(6,6,N);
  end

  for w = 1:1:N
    P1 = pcs(:,:,3*w-2:3*w);
    P1 = single(P1);
    for edge = [1 2 3 4]
      S = getPlane(P1,edge,1);
      P1Dif1 = getPlane(P1,edge,1) - getPlane(P1,edge,2);
      P1Dif2 = S(1:end-1,:) - S(2:end,:);
      P1Dif = [P1Dif1(1:end-1,:) P1Dif2];
      P1D_mu =  mean(P1Dif);
      P1D_cov = cov(double([P1Dif;dummyDiffs]));

      pixels = getPlane(P1,edge,1);
      pix(w,:,edge) = pixels(:);
      D_mu(w,:,edge) = P1D_mu;
      D_cov(:,:,w,edge) = P1D_cov;
    end
  end

  score = zeros(N,N,4,'single'); 
  %for ii = 1:1:N, P{ii} = single(P{ii}); end
  onemat = ones(partSize-1,1);

  for jj = 1:1:4
    for ii = 1:1:N-1

      P1D_mu =    D_mu(ii,:,jj);
      P1D_cov =   D_cov(:,:,ii,jj);
      p1S =       pix(ii,:,jj);
      p1S = reshape(p1S,partSize,3);
      score(ii,ii,:) = inf;             

      for kk = ii+1:1:N
        s = oppSide(jj);

        P2D_mu =    D_mu(kk,:,s);
        P2D_cov =   D_cov(:,:,kk,s);
        p2S =       pix(kk,:,s);

        % now, compute the score:
        p2S = reshape(p2S,partSize,3);

        hor = p1S - p2S;
        ver1 = p1S(1:end-1,:) - p1S(2:end,:);
        P12DIF = [hor(1:end-1,:) ver1];
        
        hor = -hor;
        ver2 = p2S(1:end-1,:) - p2S(2:end,:);
        P21DIF = [hor(1:end-1,:) ver2];
        
        D12 = (P12DIF-(onemat*P2D_mu))/(P2D_cov); 
        D12 = sum(D12 .* (P12DIF-(onemat*P2D_mu)),2);
        D21 = (P21DIF-(onemat*P1D_mu))/(P1D_cov); 
        D21 = sum(D21 .* (P21DIF-(onemat*P1D_mu)),2);

        Dist = sum(sqrt(D12)+sqrt(D21));           

        score(ii,kk,jj)           = single(Dist);
        score(kk,ii,oppSide(jj))  = single(Dist);
      end
    end
  end
  score(N,N,:) = inf;
end

function plane = getPlane(P1,edge,n)
  if      (edge==3), plane = squeeze(P1(n,:,:));
  elseif  (edge==2), plane = squeeze(P1(:,end+1-n,:));
  elseif  (edge==4), plane = squeeze(P1(end+1-n,:,:));
  else, plane = squeeze(P1(:,n,:));
  end
end