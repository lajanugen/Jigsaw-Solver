function score = score_ssd(pcs)
  partSize = size(pcs,1);
  N = size(pcs,3)/3;
  oppSide = [2 1 4 3];
  score = zeros(N,N,4,'single'); 
  %for ii = 1:1:N, P{ii} = single(P{ii}); end

  for jj = 1:1:4
    for ii = 1:1:N-1
      p1S =       getPlane(pcs(:,:,3*ii-2:3*ii),jj,1);
      score(ii,ii,:) = inf;             

      for kk = ii+1:1:N
        s = oppSide(jj);
        p2S =       getPlane(pcs(:,:,3*kk-2:3*kk),s,1);

        Dist = sum(sum((p1S-p2S).^2));

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