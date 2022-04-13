function top88_plusAM_qp_primal(nelx,nely,volfrac,penal,rmin,ft)
if nargin==0
 clc
 clear
 close all

 fact = 1;
 nelx= 60*fact;
 nely= 30*fact;
 nele = nelx*nely;
 volfrac = 0.5;
 penal = 3;
 rmin = 1.5;
 MOVE_LIMIT = 0.8;
 epsilon = 1e-5;
 ft = 2;
end
 MOVE_LIMIT = 0.8;
 epsilon = 1e-5;
nele = nelx*nely;
options = optimoptions('quadprog','Display','off');

E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%  F = sparse(2*(nelx/2*(nely+1)+(nely+1)),1,-30,2*(nely+1)*(nelx+1),1);
% %F = sparse([2*(nelx/2*(nely+1)+(nely+1))-1 2*(nelx/2*(nely+1)+(nely+1))]...,[1 1],100*[20/sqrt(2),-30-60/sqrt(2)],2*(nely+1)*(nelx+1),1);
% U = zeros(2*(nely+1)*(nelx+1),1);
%  fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% %fixeddofs = union(1:2*(nely+1):2*(nelx+1)*(nely+1),2:2*(nely+1):2*(nelx+1)*(nely+1));
% alldofs = 1:2*(nely+1)*(nelx+1);
% freedofs = setdiff(alldofs,fixeddofs);
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
 change = 1;

 nConstraints = (nely-1)*(nelx-2)+1+2*(nely-1);
 lambda{1} = ones(1,nConstraints);
 %% START ITERATION
 while loop < 200
 loop = loop + 1;
 if loop < 30
 penal = 1;
 elseif penal < 3 % after 30 iterations, before penal=3
 penal = penal*1.01;
 else
 penal = 3;
 end

 %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
  
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  
  %BUILDING A QUADRATIC SUBPROBLEM
 f= dc(:)'; % Define the term in the objective function that is multiplied ...
            %by s
            
            g=zeros(nConstraints,1);
 dg = sparse(nConstraints,nele);
 g(1) = sum(xPhys(:)) / (volfrac*nele) - 1;
 dg(1,:) = dv(:)' / (volfrac*nele); % the total volume constraint is g = p * ...x-1≤0

 
 %% Overhang constraints
 i=2;
 Size = [nely nelx];

 for ely=1:nely-1
 for elx = 1:nelx

 Xunder = [elx-1, elx, elx+1];
 Yunder = ely+1; 
 if elx > 1
 indexA = sub2ind(Size,Yunder,Xunder(1));
 AA = xPhys(indexA);
 else
 AA = 0;
 indexA = -1;
 end

 indexB = sub2ind(Size,Yunder,Xunder(2));
 BB = xPhys(indexB);
 
 if(elx<nelx)
 indexC = sub2ind(Size,Yunder,Xunder(3));
 CC = xPhys(indexC);
 else
 CC = 0;
 indexC = -1;
 end

 indexX = sub2ind(Size,ely,elx); 
 %x<a+b+c
 g(i) = xPhys(ely,elx) - AA - BB - CC;
 dg(i,indexX) = 1;
 if indexA > 0
 dg(i,indexA) = -1;
 end
 dg(i,indexB) = -1;
 if indexC > 0
 dg(i,indexC) = -1;
 end

 i=i+1;
 end
 end

 %% Objective function and c20
 c20 = 2./xPhys(:) .* abs(dc(:));
 c20 = max(1e-5,c20); % To guarantee convexity in the subproblem (Etman 2012)


 c2 = 2./xPhys(:)' .* abs(dg); 
 c2 = max(0,c2);
 Qdata = zeros(nele,1); % Var to build up Q


for i = 1:nele
    
Qdata(i) = c20(i) + lambda{loop}(:)' * c2(:,i);
end

Q=sparse(1:nele,1:nele,Qdata); % Create Q
A=dg; % Inequality constraint matrix
b=-(g); % Right side of inequality constraints
 %own implementation of 0.2 maximum move limit:
 lb = max(-xPhys(:)+epsilon,-MOVE_LIMIT * ones(nele,1));
 ub = min(1-xPhys(:), MOVE_LIMIT * ones(nele,1));

 %% SOLVING THE PRIMAL QUADRATIC SUBPROBLEM
 [s,fval,exitflag,output,lam] = quadprog(Q,f,A,b,[],[],lb,ub,[],options); % ...
          %Solve the subproblem
 lambda{loop+1} = lam.ineqlin; % Save the Lagrange multiplier of the next ...
%iteration
 xnew = xPhys + reshape(s,size(xPhys)); % Obtain the next design
 
 
 if ft == 1
 xPhys = xnew;
 elseif ft == 2
 xPhys(:) = (H*xnew(:))./Hs;
 end
 change = max(abs(xnew(:)-xPhys(:)));
 x = xnew;
 %% PRINT RESULTS
 fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
 mean(xPhys(:)),change);
 %% PLOT DENSITIES
 xDouble = [flip(xPhys,2), xPhys];

 colormap(gray); imagesc(1-xDouble); caxis([0 1]); axis equal; axis off; drawnow;
 %colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
A=1-xDouble;%1-xPhys
save A.mat A%
 end