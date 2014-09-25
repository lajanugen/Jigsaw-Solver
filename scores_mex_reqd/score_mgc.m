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
      C = ones(2*partSize,4);
      C((1:partSize),1) = (1:partSize)';
      C((1:partSize),2) = (1:partSize)' + partSize;
      C(partSize+(1:partSize),1) = (1:partSize)' + partSize;
      C(partSize+(1:partSize),2) = (1:partSize)' + 2*partSize;

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

function plane = getPlane(P1,edge,n)
  if      (edge==3), plane = squeeze(P1(n,:,:));
  elseif  (edge==2), plane = squeeze(P1(:,end+1-n,:));
  elseif  (edge==4), plane = squeeze(P1(end+1-n,:,:));
  else, plane = squeeze(P1(:,n,:));
  end
end