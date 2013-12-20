
function [] = pltMatNoScale(H)
	imagesc(logNAN10(abs(H)));
	daspect([1 1 1]);
	set(gca, "XTickLabel", "");
	set(gca, "YTickLabel", "");
end#function

function [out] = logNAN10(in)
	out = log10(in);
	out(isinf(out)) = nan;
	out(out < -16) = nan;
end#function



% so we get the same images every time
rand("seed", 324324);
mkdir("plots");

% image of original matrix
N = 10; I = eye(N);
r = rand([N,N]);
M = r * diag(rand([N, 1])*10) * r^(-1);

% image of hessenberg decomposition
pltMatNoScale(M); print -deps -tight plots/original.eps;
X = hess(M);      pltMatNoScale(X); print -deps -tight plots/hess.eps;

% two images of steps of the QR algorithm
pltMatNoScale(X); print -deps -tight plots/qr1.eps;

for i = 1:13
	[Q, R] = qr(X); X = R*Q;
end#for

qr_intermediate = X;

pltMatNoScale(X); print -deps -tight plots/qr2.eps;

% image of completed QR
X = schur(M);
pltMatNoScale(X); print -deps -tight plots/qrdone.eps;




% beyond here: faking it for the win...


% one bulge
X=hess(M); X(3,1)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge1.eps;
X=hess(M); X(6,4)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge2.eps;
X=hess(M); X(9,7)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge3.eps;

% done version
pltMatNoScale(qr_intermediate); print -deps -tight plots/bulgedone.eps;




% two bulges
X=hess(M); X(9,7)=X(1,1); X(4,2)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulges1.eps;
X=hess(M); X(7,4)=(X(1,1)+X(end,end))/2; X(7,5)=X(1,1); X(6,4)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulges2.eps;

[U, S] = schur(X(4:7, 4:7));

X(4:7, :) = U' * X(4:7, :);
X(:, 4:7) = X(:, 4:7) * U;
pltMatNoScale(X); print -deps -tight plots/bulgesspiked.eps;

% spiked version
X = qr_intermediate;
%X(1:5, 6:10) = 0; % this sort of makes it unclear
X(6:10, 1:5) = 0;
pltMatNoScale(X); print -deps -tight plots/bulgesdeflated.eps;
