source("trainBust.m");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = pltMatNoScale(H)
	pltMatNoScaleOffset(H, 1:length(H));
end#function

function [] = pltMatNoScaleOffset(H, offset)
	imagesc(offset, offset, logNAN10(abs(H)));
	daspect([1 1 1]);
	set(gca, "XTickLabel", ""); set(gca, "XTick", [])
	set(gca, "YTickLabel", ""); set(gca, "YTick", [])
end#function

function [out] = logNAN10(in)
	out = log10(in);
	out(isinf(out)) = nan;
	out(out < -16) = nan;
end#function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output from our code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand("seed", 324324); % MAGIC
rand("seed", 324); % MAGIC
eVals = 1 - 2*rand(20,1);
Q = rand(20);

A = hess(inv(Q) * diag(eVals) * Q);
[H,INFO] = trainBust(A,eVals(1:4),true,true,'m', 'plots2');

X = hess(H);
split1 = 12;
X(split1+1, split1) = 0;
%X(1:split1, split1+1:end) = 0; % zero upper right corner
pltMatNoScale(X); print -deps -tight plots2/final.eps;

% recurse top left
A = X(1:split1, 1:split1);
[H,INFO] = trainBust(A,eVals(1:4),true,true,'m', 'plots3');

A = hess(H);
p = 7;
A(p+1, p) = 0;
%A(1:p, p+1:end) = 0; % zero upper right corner
pltMatNoScale(A); print -deps -tight plots3/final.eps;

% recurse bottom right
A = X((split1+1):end, (split1+1):end);
[H,INFO] = trainBust(A,eVals(1:3),true,true,'m', 'plots4');

A = hess(H);
p = 4;
A(p+1, p) = 0;
%A(1:p, p+1:end) = 0; % zero upper right corner
pltMatNoScale(A); print -deps -tight plots4/final.eps;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of figures for "Existing techniques" section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% so we get the same images every time
rand("seed", 324326);
mkdir("plots");

% image of original matrix
N = 10; I = eye(N);
%evals = rand([N, 1])*10;
evals = [2 3 4 5 6 7 8 9 10 11];
r = rand([N,N]);
M = r * diag(evals) * r^(-1);

% image of hessenberg decomposition
pltMatNoScale(M); print -deps -tight plots/original.eps;
X = hess(M);      pltMatNoScale(X); print -deps -tight plots/hess.eps;

% two images of steps of the QR algorithm
[Q, R] = qr(X); X = R*Q;
qr_first_iteration_done = X;
pltMatNoScale(X); print -deps -tight plots/qr1.eps;

for i = 1:150
	shiftm = eye(size(X)) * X(end,end-1);
	[Q, R] = qr(X - shiftm); X = R*Q + shiftm;
end#for

pltMatNoScale(X); print -deps -tight plots/qr2.eps;

% image of completed QR
X = schur(M);
pltMatNoScale(X); print -deps -tight plots/qrdone.eps;



% beyond here: faking it for the win...


% one bulge
X=hess(M); X(3,1)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge1.eps;
X=hess(M); X(5,3)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge2.eps;
X=hess(M); X(7,5)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulge3.eps;

% done version
pltMatNoScale(qr_first_iteration_done); print -deps -tight plots/bulgedone.eps;




% two bulges
X=hess(M); X(9,7)=X(1,1); X(4,2)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulges1.eps;
X=hess(M); X(7,4)=(X(1,1)+X(end,end))/2; X(7,5)=X(1,1); X(6,4)=X(end,end); pltMatNoScale(X); print -deps -tight plots/bulges2.eps;

[U, S] = schur(X(4:7, 4:7));

X(4:7, :) = U' * X(4:7, :);
X(:, 4:7) = X(:, 4:7) * U;
pltMatNoScale(X);
print -deps -tight plots/bulgesspiked.eps;

% spiked version
X = qr_first_iteration_done;
%X(1:5, 6:10) = 0; % this sort of makes it unclear
X(6:10, 1:5) = 0;
pltMatNoScale(X); print -deps -tight plots/bulgesdeflated.eps;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCHUR DECOMPOSITION IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_poster_matrix_figure(A, filename)
	pltMatNoScale(A);
	create_poster_figure_from_gcf(filename);
end#function

function create_poster_figure_from_gcf(filename)
	print(['plots/' filename '.eps'], '-deps', '-tight');
end#function

rand('seed', 324);
r = rand([15 15]);
A = r * diag(rand([1 15])) * r^(-1);
Ahess = hess(A);

sub = 6:10;
[smallQ, ~] = schur(Ahess(sub, sub));

# embed in identity matrix
Q = eye(length(A));
Q(sub,sub) = smallQ;

Az = zeros(size(A));
Qz = Az;

Az(sub, sub) = Ahess(sub, sub);
Qz(sub, sub) = smallQ;

create_poster_matrix_figure(A, 'schur_A');
create_poster_matrix_figure(Ahess, 'schur_Ahess');
create_poster_matrix_figure(Az, 'schur_Az');
create_poster_matrix_figure(Qz, 'schur_Qz');
create_poster_matrix_figure(Q, 'schur_Q');
create_poster_matrix_figure(Q', 'schur_Qt');
create_poster_matrix_figure(Q' * Ahess * Q, 'schur_final');

